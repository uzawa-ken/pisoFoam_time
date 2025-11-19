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

//
        const label nCells = pEqn.diag().size();   // = 行数
        const label nUpper = pEqn.upper().size();  // = 内部面の数 (下三角と同じ)

        const label nnz = nCells + 2*nUpper;

        Info<< "p-matrix: nCells="   << nCells
            << "  nUpper="           << nUpper
            << "  nnz="              << nnz
            << nl;
//    	
    	
        ///============
        const bool doGnnWrite =
            writePressureSystem
         && piso.finalNonOrthogonalIter()
         && (runTime.timeIndex() % gnnWriteInterval == 0);

        // ===== ここで A, b, 座標などをスナップショットとして出力する =====
        if (doGnnWrite)
        {
            // 1. A（lduMatrix 構造）を取得
            const lduMatrix& A = pEqn; // ← pEqn は fvScalarMatrix なので、そのまま lduMatrix として見てよい

            const scalarField& diag  = A.diag();  // 対角成分 A_ii
            const scalarField& lower = A.lower(); // 下三角（内部フェイスごと）
            const scalarField& upper = A.upper(); // 上三角（内部フェイスごと）

            const lduAddressing& addr = A.lduAddr();
            const labelUList& lowerAddr = addr.lowerAddr();
            const labelUList& upperAddr = addr.upperAddr();

            // 2. 右辺ベクトル b
            const scalarField& b = pEqn.source(); 

            // 3. セル中心座標（グラフノード特徴用に）
            const vectorField& C = mesh.C();

            // 4. 出力ファイル名（例：<case>/gnn/pEqn_0.5.dat など）
            fileName outDir(runTime.path()/"gnn");
            mkDir(outDir);
            fileName outName(outDir/("pEqn_" + runTime.timeName() + ".dat"));

            OFstream os(outName);

            label nCells = mesh.nCells();
            label nFaces = lowerAddr.size();

            // --- ヘッダ：セル数と内部フェイス数
            os  << "nCells " << nCells << '\n'
                << "nFaces " << nFaces << '\n';

            // --- ノード（セル）情報：id, x, y, z, diag, b
            os << "CELLS\n";
            for (label i=0; i<nCells; ++i)
            {
                const vector& ci = C[i];
                os  << i << ' '
                    << ci.x() << ' ' << ci.y() << ' ' << ci.z() << ' '
                    << diag[i] << ' '
                    << b[i] << '\n';
            }

            // --- エッジ情報：faceId, lowerCell, upperCell, lowerCoeff, upperCoeff
            os << "EDGES\n";
            for (label f=0; f<nFaces; ++f)
            {
                os  << f << ' '
                    << lowerAddr[f] << ' ' << upperAddr[f] << ' '
                    << lower[f] << ' ' << upper[f] << '\n';
            }
        }

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
                fileName outName(outDir/("x_" + runTime.timeName() + ".dat"));
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
                fileName rName(outDir/("rTrue_" + runTime.timeName() + ".dat"));
                OFstream osR(rName);

                for (label i=0; i<nCells; ++i)
                {
                    osR << i << ' ' << r[i] << '\n';
                }
            }
        }
        ///============

        if (piso.finalNonOrthogonalIter())
        {
            cpuTime __tAsm2;
            phi = phiHbyA - pEqn.flux();
            accPEqnAsm += __tAsm2.elapsedCpuTime();
        }
    }

    // ===== continuityErrs + 速度補正（assemble 側に加算） ===================
    cpuTime __tAsm3;
    #include "continuityErrs.H"

    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
    accPEqnAsm += __tAsm3.elapsedCpuTime();
}

