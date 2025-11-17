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

        // 解ベクトル x を出力
        if (doGnnWrite)
        {
            const scalarField& x = p.internalField();

            fileName outDir(runTime.path()/"gnn");
            mkDir(outDir);
            fileName outName(outDir/("x_" + runTime.timeName() + ".dat"));
            OFstream osX(outName);

            label nCells = mesh.nCells();
            for (label i=0; i<nCells; ++i)
            {
                osX << i << ' ' << x[i] << '\n';
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

