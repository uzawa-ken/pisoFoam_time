// pEqn_timed.H  (PISO + MRF あり版)
// 前提：これを include する .C 側で
//   #include "cpuTime.H"
//   double accPEqnAsm=0.0, accPEqnSol=0.0;
// が定義済みであること。

{
    // ===== Assemble: rAU, HbyA, phiHbyA, adjustPhi, constrainPressure =====
    cpuTime __tAsm0;

    volScalarField rAU(1.0/UEqn.A());

    volVectorField HbyA
    (
        constrainHbyA(rAU*UEqn.H(), U, p)
    );

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi))
    );

    MRF.makeRelative(phiHbyA);

    adjustPhi(phiHbyA, U, p);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAU, MRF);

    accPEqnAsm += __tAsm0.elapsedCpuTime();

    // ================= Non-orthogonal pressure corrector loop ==============
    while (piso.correctNonOrthogonal())
    {
        // assemble matrix
        cpuTime __tAsm1;

        ///------ A_p * p = b_p
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
        );
        pEqn.setReference(pRefCell, pRefValue);
        ///------

        accPEqnAsm += __tAsm1.elapsedCpuTime();

/*
    	// 行列サイズログ
        const label nCellsMat = pEqn.diag().size();   // = 行数
        const label nUpper    = pEqn.upper().size();  // = 内部面の数
        const label nnz       = nCellsMat + 2*nUpper;

        Info<< "p-matrix: nCellsMat=" << nCellsMat
            << "  nUpper=" << nUpper
            << "  nnz="    << nnz
            << nl;
*/
        // ローカル
        const label nCellsLocal = pEqn.diag().size();
        const label nUpperLocal = pEqn.upper().size();
        const label nnzLocal    = nCellsLocal + 2*nUpperLocal;

        // グローバル合計（sum）
        const label nCellsGlobal = returnReduce(nCellsLocal, sumOp<label>());
        const label nUpperGlobal = returnReduce(nUpperLocal, sumOp<label>());
        const label nnzGlobal    = returnReduce(nnzLocal,   sumOp<label>());

        // master だけに全体を出す
        if (Pstream::master())
        {
            Info<< "p-matrix global: nCells=" << nCellsGlobal
                << "  nUpper=" << nUpperGlobal
                << "  nnz="    << nnzGlobal
                << nl;
        }

        // 各ランクごとの値はさっきの Pout で
        Pout<< "p-matrix(proc " << Pstream::myProcNo() << "): "
            << "nCellsMat=" << nCellsLocal
            << "  nUpper="  << nUpperLocal
            << "  nnz="     << nnzLocal
            << nl;
    	
        ///============
        const bool doGnnWrite =
            writePressureSystem
         && piso.finalNonOrthogonalIter()
         && (runTime.timeIndex() % gnnWriteInterval == 0);

        // ===== ここで A, b, 座標などをスナップショットとして出力する =====
        // ===== ここで A, b, 座標などをスナップショットとして出力する =====
        if (doGnnWrite)
        {
            // 1. A（lduMatrix 構造）を取得
            const lduMatrix& A = pEqn; // ← pEqn は fvScalarMatrix なので、そのまま lduMatrix として見てよい

            const scalarField& diag  = A.diag();  // 対角成分 A_ii
            const scalarField& lower = A.lower(); // 下三角（内部フェイスごと）
            const scalarField& upper = A.upper(); // 上三角（内部フェイスごと）

            const lduAddressing& addr     = A.lduAddr();
            const labelUList&    lowerAddr = addr.lowerAddr();
            const labelUList&    upperAddr = addr.upperAddr();

            // 2. 右辺ベクトル b
            const scalarField& b = pEqn.source(); 

            // 3. セル中心座標（グラフノード特徴用に）
            const vectorField& C = mesh.C();

            // 4. 出力ファイル名（例：<case>/gnn/pEqn_0.5.dat など）
            fileName outDir(runTime.path()/"gnn");
            mkDir(outDir);
            fileName outName
            (
                outDir
              / (
                    "pEqn_"
                  + runTime.timeName()
                  + "_rank" + Foam::name(Pstream::myProcNo())
                  + ".dat"
                )
            );

            // メッシュサイズ
            const label nCells = mesh.nCells();
            const label nFaces = lowerAddr.size();

            // 境界メッシュ（CoCell 計算と WALL_FACES 出力用）
            const fvBoundaryMesh&  fvb = mesh.boundary();
            const polyBoundaryMesh& pbm = mesh.boundaryMesh();

            // ---------------------------------------------------------
            // 1) cellQuality を使って「セル skewness / non-orthogonality」
            // ---------------------------------------------------------
            cellQuality cq(mesh);

            tmp<scalarField> tCellNonOrtho = cq.nonOrthogonality();
            const scalarField& cellNonOrtho = tCellNonOrtho();    // [deg]

            tmp<scalarField> tCellSkew = cq.skewness();
            const scalarField& cellSkew = tCellSkew();            // 無次元

            // ---------------------------------------------------------
            // 2) とりあえずの cell aspect ratio（バウンディングボックス比）
            // ---------------------------------------------------------
            const cellList&  cells = mesh.cells();
            const faceList&  faces = mesh.faces();
            const pointField& pts  = mesh.points();

            scalarField cellAspect(nCells, 1.0);

            for (label celli=0; celli<nCells; ++celli)
            {
                const labelList& cFaces = cells[celli];

                scalar minX = GREAT, maxX = -GREAT;
                scalar minY = GREAT, maxY = -GREAT;
                scalar minZ = GREAT, maxZ = -GREAT;

                forAll(cFaces, fi)
                {
                    const label facei = cFaces[fi];
                    const face& f = faces[facei];

                    forAll(f, pi)
                    {
                        const point& p = pts[f[pi]];
                        minX = Foam::min(minX, p.x());
                        maxX = Foam::max(maxX, p.x());
                        minY = Foam::min(minY, p.y());
                        maxY = Foam::max(maxY, p.y());
                        minZ = Foam::min(minZ, p.z());
                        maxZ = Foam::max(maxZ, p.z());
                    }
                }

                const scalar dx = maxX - minX;
                const scalar dy = maxY - minY;
                const scalar dz = maxZ - minZ;

                const scalar minDim = Foam::max(SMALL, Foam::min(dx, Foam::min(dy, dz)));
                const scalar maxDim = Foam::max(dx, Foam::max(dy, dz));

                cellAspect[celli] = maxDim/minDim;   // 簡易 aspect ratio
            }

            // ---------------------------------------------------------
            // 3) diagContrast: 近傍セルの対角成分の max/min 比
            // ---------------------------------------------------------
            scalarField diagMax(nCells, 0.0);
            scalarField diagMin(nCells, GREAT);

            // まず自セルで初期化
            for (label i = 0; i < nCells; ++i)
            {
                const scalar d = Foam::mag(diag[i]);
                diagMax[i] = d;
                diagMin[i] = d;
            }

            // 近傍セル（lower/upper アドレス）を見て更新
            for (label f = 0; f < nFaces; ++f)
            {
                const label own = lowerAddr[f];
                const label nei = upperAddr[f];

                const scalar dOwn = Foam::mag(diag[own]);
                const scalar dNei = Foam::mag(diag[nei]);

                diagMax[own] = Foam::max(diagMax[own], dNei);
                diagMin[own] = Foam::min(diagMin[own], dNei);

                diagMax[nei] = Foam::max(diagMax[nei], dOwn);
                diagMin[nei] = Foam::min(diagMin[nei], dOwn);
            }

            scalarField diagContrast(nCells, 1.0);
            for (label i = 0; i < nCells; ++i)
            {
                // 均一なら ~1, スケールが不均一なところほど >1
                diagContrast[i] = diagMax[i]/(diagMin[i] + SMALL);
            }

            // ---------------------------------------------------------
            // 4) セル体積に基づく代表長さとサイズジャンプ
            // ---------------------------------------------------------
            const scalarField& V = mesh.V();
            scalarField cellSize(nCells, 0.0);

            for (label celli = 0; celli < nCells; ++celli)
            {
                // 代表長さ h_i ≒ V_i^(1/3)
                cellSize[celli] = Foam::pow(V[celli], 1.0/3.0);
            }

            // 隣接セルとのサイズジャンプ:
            //   sizeJump_i = max_j max(h_i/h_j, h_j/h_i)
            scalarField sizeJump(nCells, 1.0);

            for (label f = 0; f < nFaces; ++f)
            {
                const label own = lowerAddr[f];
                const label nei = upperAddr[f];

                const scalar hi = cellSize[own];
                const scalar hj = cellSize[nei];

                scalar ratio = hi/(hj + VSMALL);
                if (ratio < 1.0) ratio = 1.0/ratio;

                sizeJump[own] = Foam::max(sizeJump[own], ratio);
                sizeJump[nei] = Foam::max(sizeJump[nei], ratio);
            }

            // ---------------------------------------------------------
            // 4’) セルごとのクーラン数 CoCell を計算
            //      CoCell_i ≒ 0.5 * Δt / V_i * Σ_{faces∈cell i} |phiHbyA_face|
            // ---------------------------------------------------------
            scalarField CoCell(nCells, 0.0);

            const scalarField& phiIn = phiHbyA.internalField();
            const auto&        phiBf = phiHbyA.boundaryField();

            const label  nInternalFaces = mesh.nInternalFaces();
            const scalar deltaT         = runTime.deltaTValue();

            for (label celli = 0; celli < nCells; ++celli)
            {
                const labelList& cFaces = cells[celli];
                scalar sumMagPhi = 0.0;

                forAll(cFaces, fi)
                {
                    const label facei = cFaces[fi];

                    scalar phiVal = 0.0;
                    if (facei < nInternalFaces)
                    {
                        // 内部フェイス
                        phiVal = phiIn[facei];
                    }
                    else
                    {
                        // 境界フェイス
                        const label patchi     = pbm.whichPatch(facei);
                        const label localFacei = facei - pbm[patchi].start();
                        phiVal = phiBf[patchi][localFacei];
                    }

                    sumMagPhi += Foam::mag(phiVal);
                }

                CoCell[celli] = 0.5*sumMagPhi*deltaT/(V[celli] + VSMALL);
            }

            // ---------------------------------------------------------
            // 5) ファイル出力
            // ---------------------------------------------------------
            OFstream os(outName);

            // --- ヘッダ：セル数と内部フェイス数
            os  << "nCells " << nCells << '\n'
                << "nFaces " << nFaces << '\n';

            // --- ノード（セル）情報：
            // id, x, y, z, diag, b, skew, nonOrtho, aspectRatio,
            // diagContrast, cellSize, sizeJump, CoCell
            os << "CELLS\n";
            for (label i = 0; i < nCells; ++i)
            {
                const vector& ci = C[i];
                os  << i << ' '
                    << ci.x() << ' ' << ci.y() << ' ' << ci.z() << ' '
                    << diag[i]        << ' '
                    << b[i]           << ' '
                    << cellSkew[i]    << ' '
                    << cellNonOrtho[i]<< ' '
                    << cellAspect[i]  << ' '
                    << diagContrast[i]<< ' '
                    << V[i]           << ' '
                    << cellSize[i]    << ' '
                    << sizeJump[i]    << ' '
                    << CoCell[i]
                    << '\n';
            }

            // --- エッジ情報：faceId, lowerCell, upperCell, lowerCoeff, upperCoeff
            os << "EDGES\n";
            for (label f = 0; f < nFaces; ++f)
            {
                os  << f << ' '
                    << lowerAddr[f] << ' ' << upperAddr[f] << ' '
                    << lower[f] << ' ' << upper[f] << '\n';
            }

            // ---------------------------------------------------------
            // 壁境界フェイス情報の出力
            //   セクション名: WALL_FACES
            //   フォーマット:
            //     patchIndex  patchName  localFaceId  globalFaceId  ownerCell
            //     cx  cy  cz  Sfx  Sfy  Sfz
            // ---------------------------------------------------------
            os << "WALL_FACES\n";

            const vectorField& Cf = mesh.Cf();
            const vectorField& Sf = mesh.Sf();

            forAll(fvb, patchi)
            {
                const fvPatch& fvp = fvb[patchi];

                if (fvp.type() == "wall")
                {
                    const polyPatch& pp        = pbm[patchi];
                    const label      start     = pp.start();
                    const label      size      = pp.size();
                    const labelList& faceCells = pp.faceCells();
                    const word&      patchName = fvp.name();

                    for (label i = 0; i < size; ++i)
                    {
                        const label faceI  = start + i;   // グローバル face index
                        const label cellI  = faceCells[i];
                        const vector& cf   = Cf[faceI];   // face center
                        const vector& sf   = Sf[faceI];   // face area vector (n*|S|)

                        os << patchi      << ' '   // patch index
                           << patchName   << ' '   // patch name
                           << i           << ' '   // local face index on the patch
                           << faceI       << ' '   // global face index
                           << cellI       << ' '   // owner cell
                           << cf.x() << ' ' << cf.y() << ' ' << cf.z() << ' '
                           << sf.x() << ' ' << sf.y() << ' ' << sf.z()
                           << '\n';
                    }
                }
            }

            // ---------------------------------------------------------
            // 6) 係数行列 A を CSR 形式で出力（各プロセスローカル）
            //     - nRows = nCells
            //     - nnz   = nCells (diag) + 2*nFaces (上下三角)
            //     - 行ごとに diag + 近傍セルを詰める
            // ---------------------------------------------------------
            // ---------------------------------------------------------
            // 6) 係数行列 A を CSR 形式で出力（各プロセスローカル）
            //     - nRows = nCells
            //     - nnz   = nCells (diag) + 2*nFaces (上下三角)
            //     - 行ごとに diag + 近傍セルを詰める
            // ---------------------------------------------------------
            {
                const label nRows = nCells;
                const label nnz   = nCells + 2*nFaces;

                // ★ 修正: A_csr_<time>_rank<proc>.dat に出力する
                fileName csrName
                (
                    outDir
                  / (
                        "A_csr_"
                      + runTime.timeName()
                      + "_rank" + Foam::name(Pstream::myProcNo())
                      + ".dat"
                    )
                );

                OFstream osCsr(csrName);

                osCsr << "nRows " << nRows << '\n'
                      << "nCols " << nRows << '\n'
                      << "nnz "   << nnz   << '\n';

                // 行ごとの非ゼロ数を数える（対角1 + 近傍）
                List<label> rowNNZ(nRows, 1);   // まず対角分を1としてセット
                for (label f = 0; f < nFaces; ++f)
                {
                    const label own = lowerAddr[f];
                    const label nei = upperAddr[f];
                    ++rowNNZ[own];
                    ++rowNNZ[nei];
                }

                // rowPtr を prefix-sum で構成
                List<label> rowPtr(nRows + 1, 0);
                for (label i = 0; i < nRows; ++i)
                {
                    rowPtr[i+1] = rowPtr[i] + rowNNZ[i];
                }

                osCsr << "ROW_PTR\n";
                for (label i = 0; i <= nRows; ++i)
                {
                    osCsr << rowPtr[i] << (i < nRows ? ' ' : '\n');
                }

                // colInd / val を埋める
                List<label>  colInd(nnz);
                List<scalar> val(nnz);

                // 行ごとの現在書き込み位置（rowPtr のコピー）
                List<label> rowPos(rowPtr);

                // まず対角成分を追加
                for (label i = 0; i < nRows; ++i)
                {
                    const label idx = rowPos[i]++;
                    colInd[idx] = i;
                    val[idx]    = diag[i];
                }

                // 次にオフ対角（上下三角）を追加
                for (label f = 0; f < nFaces; ++f)
                {
                    const label own = lowerAddr[f];
                    const label nei = upperAddr[f];

                    // A(own, nei) = upper[f]
                    {
                        const label idx = rowPos[own]++;
                        colInd[idx] = nei;
                        val[idx]    = upper[f];
                    }

                    // A(nei, own) = lower[f]
                    {
                        const label idx = rowPos[nei]++;
                        colInd[idx] = own;
                        val[idx]    = lower[f];
                    }
                }

                // colInd を出力
                osCsr << "COL_IND\n";
                for (label k = 0; k < nnz; ++k)
                {
                    osCsr << colInd[k] << (k+1 < nnz ? ' ' : '\n');
                }

                // 値を出力
                osCsr << "VALUES\n";
                for (label k = 0; k < nnz; ++k)
                {
                    osCsr << val[k] << (k+1 < nnz ? ' ' : '\n');
                }
            }
        } // if (doGnnWrite)
    	
        // solve
        cpuTime __tSol;
        pEqn.solve(p.select(piso.finalInnerIter()));
        accPEqnSol += __tSol.elapsedCpuTime();

        // 解ベクトル x と「真の残差 r = A_full x - b」を出力
        if (doGnnWrite)
        {
            const scalarField& x = p.internalField();
            const scalarField& b = pEqn.source();
            const label nCells = mesh.nCells();

            // ---- x を出力
            {
                fileName outDir(runTime.path()/"gnn");
                mkDir(outDir);
                fileName outName
                (
                    outDir
                  / (
                        "x_"
                      + runTime.timeName()
                      + "_rank" + Foam::name(Pstream::myProcNo())
                      + ".dat"
                    )
                );

            	OFstream osX(outName);

                for (label i=0; i<nCells; ++i)
                {
                    osX << i << ' ' << x[i] << '\n';
                }
            }

            // ---- fvScalarMatrix::residual() を使って「インターフェイス込み」の r = A x - b を取得
            // v2412 では residual() の戻り値は tmp<Field<scalar>> = tmp<scalarField>
            tmp<scalarField> tR = pEqn.residual();
            const scalarField& r = tR();   // セルごとの残差 r_i（長さ nCells）

            scalar sumR2  = 0.0;
            scalar sumB2  = 0.0;
            scalar sumAx2 = 0.0;
            scalar maxAbsR = 0.0;

            for (label i=0; i<nCells; ++i)
            {
                const scalar ri = r[i];
                const scalar bi = b[i];
                const scalar Axi = ri + bi;  // r = Ax - b → Ax = r + b

                sumR2  += ri*ri;
                sumB2  += bi*bi;
                sumAx2 += Axi*Axi;
                maxAbsR = Foam::max(maxAbsR, Foam::mag(ri));
            }

            const scalar normR  = Foam::sqrt(sumR2);
            const scalar normB  = Foam::sqrt(sumB2) + SMALL;
            const scalar normAx = Foam::sqrt(sumAx2) + SMALL;

            Info<< "pEqn true(Ax-b) via residual(): "
                << "||r||_2="        << normR
                << ", max|r_i|="     << maxAbsR
                << ", ||r||/||b||="  << normR/normB
                << ", ||r||/||Ax||=" << normR/normAx
                << nl;

            // ---- セルごとの r をファイルに出力
            {
                fileName outDir(runTime.path()/"gnn");
                mkDir(outDir);
            	fileName rName
            	(
            	    outDir
            	  / (
            	        "rTrue_"
            	      + runTime.timeName()
            	      + "_rank" + Foam::name(Pstream::myProcNo())
            	      + ".dat"
            	    )
            	);

            	OFstream osR(rName);

                for (label i=0; i<nCells; ++i)
                {
                    osR << i << ' ' << r[i] << '\n';
                }
            }

            // ---------------------------------------------------------
            // 条件数評価用データの出力
            //
            // 理論的背景:
            //   |e|/|x_true| <= kappa(A) * |r|/|b|
            // ここで e = x_pred - x_true, r = A*x_pred - b
            //
            // 出力内容:
            //   1. ゲルシュゴリンの円定理に基づく固有値境界
            //   2. べき乗法による最大固有値推定
            //   3. 必要なノルム
            // ---------------------------------------------------------
            {
                const lduMatrix& A = pEqn;
                const scalarField& diag  = A.diag();
                const scalarField& lower = A.lower();
                const scalarField& upper = A.upper();
                const lduAddressing& addr = A.lduAddr();
                const labelUList& lowerAddr = addr.lowerAddr();
                const labelUList& upperAddr = addr.upperAddr();
                const label nFacesLocal = lowerAddr.size();

                // ----- ゲルシュゴリンの円定理に基づく固有値境界 -----
                // 各行の非対角成分の絶対値和を計算
                scalarField rowOffDiagSum(nCells, 0.0);

                for (label f = 0; f < nFacesLocal; ++f)
                {
                    const label own = lowerAddr[f];
                    const label nei = upperAddr[f];

                    rowOffDiagSum[own] += Foam::mag(upper[f]);
                    rowOffDiagSum[nei] += Foam::mag(lower[f]);
                }

                // ゲルシュゴリン境界: λ ∈ [a_ii - R_i, a_ii + R_i]
                // 対称正定値行列なので固有値は正
                scalar lambdaMinGershgorin = GREAT;
                scalar lambdaMaxGershgorin = -GREAT;

                for (label i = 0; i < nCells; ++i)
                {
                    const scalar d = diag[i];  // 対角成分（負値、ラプラシアン）
                    const scalar R = rowOffDiagSum[i];

                    // ポアソン行列の場合、対角成分は負（-Σ係数）
                    // 固有値の範囲: [|d| - R, |d| + R]
                    const scalar absD = Foam::mag(d);
                    const scalar lambdaMin_i = Foam::max(SMALL, absD - R);
                    const scalar lambdaMax_i = absD + R;

                    lambdaMinGershgorin = Foam::min(lambdaMinGershgorin, lambdaMin_i);
                    lambdaMaxGershgorin = Foam::max(lambdaMaxGershgorin, lambdaMax_i);
                }

                // 条件数の上界（ゲルシュゴリン推定）
                const scalar kappaGershgorin = lambdaMaxGershgorin / (lambdaMinGershgorin + SMALL);

                // ----- べき乗法による最大固有値推定 -----
                // A*v を計算する関数（lduMatrix構造を使用）
                // 注: OpenFOAMのラプラシアン行列は -div(grad) なので負定値
                //     固有値評価のため絶対値を使う

                const label powerIterMax = 30;
                const scalar powerIterTol = 1e-4;

                scalarField v(nCells, 1.0/Foam::sqrt(scalar(nCells)));
                scalarField Av(nCells, 0.0);
                scalar lambdaMaxPower = 0.0;

                for (label iter = 0; iter < powerIterMax; ++iter)
                {
                    // Av = A * v （行列ベクトル積）
                    Av = 0.0;

                    // 対角成分
                    for (label i = 0; i < nCells; ++i)
                    {
                        Av[i] = diag[i] * v[i];
                    }

                    // 非対角成分
                    for (label f = 0; f < nFacesLocal; ++f)
                    {
                        const label own = lowerAddr[f];
                        const label nei = upperAddr[f];

                        Av[own] += upper[f] * v[nei];
                        Av[nei] += lower[f] * v[own];
                    }

                    // レイリー商: λ ≈ (v^T A v) / (v^T v)
                    // ただし |Av| を使って正の固有値として扱う
                    scalar vTAv = 0.0;
                    scalar vTv = 0.0;
                    for (label i = 0; i < nCells; ++i)
                    {
                        vTAv += v[i] * Av[i];
                        vTv += v[i] * v[i];
                    }
                    const scalar lambdaNew = Foam::mag(vTAv) / (vTv + SMALL);

                    // 収束判定
                    if (iter > 0 && Foam::mag(lambdaNew - lambdaMaxPower) / (lambdaMaxPower + SMALL) < powerIterTol)
                    {
                        lambdaMaxPower = lambdaNew;
                        break;
                    }
                    lambdaMaxPower = lambdaNew;

                    // ベクトルの正規化: v = Av / |Av|
                    scalar normAv = 0.0;
                    for (label i = 0; i < nCells; ++i)
                    {
                        normAv += Av[i] * Av[i];
                    }
                    normAv = Foam::sqrt(normAv) + SMALL;

                    for (label i = 0; i < nCells; ++i)
                    {
                        v[i] = Av[i] / normAv;
                    }
                }

                // ----- 対角スケーリング後の条件数推定 -----
                // D^(-1/2) A D^(-1/2) の条件数（対角前処理後）
                scalar diagMin = GREAT;
                scalar diagMax = -GREAT;
                for (label i = 0; i < nCells; ++i)
                {
                    const scalar absD = Foam::mag(diag[i]);
                    diagMin = Foam::min(diagMin, absD);
                    diagMax = Foam::max(diagMax, absD);
                }
                const scalar diagRatio = diagMax / (diagMin + SMALL);

                // ----- ノルム計算 -----
                // |x_true|_2
                scalar sumX2 = 0.0;
                for (label i = 0; i < nCells; ++i)
                {
                    sumX2 += x[i] * x[i];
                }
                const scalar normXtrue = Foam::sqrt(sumX2) + SMALL;

                // ----- 出力 -----
                fileName outDir(runTime.path()/"gnn");
                mkDir(outDir);
                fileName condName
                (
                    outDir
                  / (
                        "conditionEval_"
                      + runTime.timeName()
                      + "_rank" + Foam::name(Pstream::myProcNo())
                      + ".dat"
                    )
                );

                OFstream osCond(condName);

                osCond << "# Condition number evaluation data\n";
                osCond << "# |e|/|x_true| <= kappa(A) * |r|/|b|\n";
                osCond << "#\n";
                osCond << "nCells " << nCells << '\n';
                osCond << '\n';

                osCond << "# Gershgorin bounds (eigenvalue estimates)\n";
                osCond << "lambda_min_gershgorin " << lambdaMinGershgorin << '\n';
                osCond << "lambda_max_gershgorin " << lambdaMaxGershgorin << '\n';
                osCond << "kappa_gershgorin " << kappaGershgorin << '\n';
                osCond << '\n';

                osCond << "# Power iteration estimate (max eigenvalue)\n";
                osCond << "lambda_max_power " << lambdaMaxPower << '\n';
                osCond << '\n';

                osCond << "# Diagonal scaling info\n";
                osCond << "diag_min " << diagMin << '\n';
                osCond << "diag_max " << diagMax << '\n';
                osCond << "diag_ratio " << diagRatio << '\n';
                osCond << '\n';

                osCond << "# Norms for error bound evaluation\n";
                osCond << "norm_b " << normB << '\n';
                osCond << "norm_x_true " << normXtrue << '\n';
                osCond << "norm_r " << normR << '\n';
                osCond << "max_abs_r " << maxAbsR << '\n';
                osCond << '\n';

                osCond << "# Relative quantities\n";
                osCond << "rel_residual_norm_r_over_b " << normR/normB << '\n';
                osCond << "rel_residual_norm_r_over_Ax " << normR/normAx << '\n';
                osCond << '\n';

                osCond << "# Error bound (upper bound on relative error)\n";
                osCond << "# |e|/|x_true| <= kappa * |r|/|b|\n";
                osCond << "error_bound_gershgorin " << kappaGershgorin * normR/normB << '\n';

                Info<< "Condition number evaluation: "
                    << "kappa_gershgorin=" << kappaGershgorin
                    << ", lambda_max_power=" << lambdaMaxPower
                    << ", |r|/|b|=" << normR/normB
                    << ", error_bound=" << kappaGershgorin * normR/normB
                    << nl;
            }
        }
        ///============

/// div flux
        if (piso.finalNonOrthogonalIter())
        {
            cpuTime __tAsm2;
            phi = phiHbyA - pEqn.flux();
            accPEqnAsm += __tAsm2.elapsedCpuTime();

            // ===== ここから：セル間フラックス発散を出力 =====
            if (doGnnWrite)
            {
                // phiHbyA: 予測フラックス（u* に対応）
                // phi    : 圧力修正後フラックス（最終 U に対応）

                // OpenFOAM の離散スキームそのままの fvc::div を使う
                tmp<volScalarField> tDivPhiHbyA = fvc::div(phiHbyA);
                const scalarField& divPhiHbyA =
                    tDivPhiHbyA().internalField();   // 各セルの ∇·(u*)

                tmp<volScalarField> tDivPhi = fvc::div(phi);
                const scalarField& divPhi =
                    tDivPhi().internalField();       // 各セルの ∇·u^(n+1)

                const scalarField& V = mesh.V();
                const label nCells = mesh.nCells();

                fileName outDir(runTime.path()/"gnn");
                mkDir(outDir);

                // 例: divPhi_10_rank7.dat
                fileName divName
                (
                    outDir
                  / (
                        "divPhi_"
                      + runTime.timeName()
                      + "_rank" + Foam::name(Pstream::myProcNo())
                      + ".dat"
                    )
                );

                OFstream osDiv(divName);

                // フォーマット:
                //   cellId  divPhiHbyA  divPhi  volume
                //   0       ...         ...     ...
                //   1       ...         ...     ...
                //   ...
                for (label i=0; i<nCells; ++i)
                {
                    osDiv
                        << i << ' '
                        << divPhiHbyA[i] << ' '
                        << divPhi[i]     << ' '
                        << V[i]
                        << '\n';
                }
            }
            // ===== ここまで：セル間フラックス発散を出力 =====
        }
///
    	
    }

    // ===== continuityErrs + 速度補正（assemble 側に加算） ===================
    cpuTime __tAsm3;
    #include "continuityErrs.H"

    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
    accPEqnAsm += __tAsm3.elapsedCpuTime();
}

