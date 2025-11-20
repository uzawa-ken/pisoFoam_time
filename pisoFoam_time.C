/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pisoFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, turbulent flow, using the PISO
    algorithm.

    \heading Solver details
    The solver uses the PISO algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}} + \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
    \endvartable

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "cpuTime.H"
#include <sstream>
inline std::string fmt8(double x){ std::ostringstream os; os.setf(std::ios::fixed); os.precision(8); os<<x; return os.str(); }
#include "cellQuality.H"   // skewness, nonOrthogonality 用


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow,"
        " using the PISO algorithm."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

///
    // createFields.H の後、時間ループの前あたり
    const Switch writePressureSystem
    (
        mesh.solutionDict().lookupOrDefault<Switch>("writePressureSystem", false)
    );

    // 何ステップおきに書くか（デフォルト 1 = 毎ステップ）
    const label gnnWriteInterval
    (
        mesh.solutionDict().lookupOrDefault<label>("gnnWriteInterval", 10)
    );
///

    while (runTime.loop())
    {
        cpuTime tStep; // このタイムステップ総時間
        double accUEqn=0.0, accPEqnAsm=0.0, accPEqnSol=0.0, accTurb=0.0, accWrite=0.0;
        auto pct = [&](double v)->std::string{
            const double tot=tStep.elapsedCpuTime();
            std::ostringstream os; os.setf(std::ios::fixed); os.precision(2);
            os << (tot>0.0 ? v*100.0/tot : 0.0); return os.str();
        };
    
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Update settings from the control dictionary
        piso.read();

        // Pressure-velocity PISO corrector
        {
            // UEqn を同一スコープに展開しつつ時間計測
            cpuTime __tUE;
            #include "UEqn.H"   // ← 行頭に単独で置く。前に何も書かない
            accUEqn += __tUE.elapsedCpuTime();

            // --- PISO loop
            while (piso.correct())
            {
                // assemble/solve を分解計測する版
                #include "pEqn_timed.H"
            }
        }

        laminarTransport.correct();
        { cpuTime __tTb; turbulence->correct(); accTurb += __tTb.elapsedCpuTime(); }

        { cpuTime __tW; runTime.write(); accWrite += __tW.elapsedCpuTime(); }

        const double stepTotal = tStep.elapsedCpuTime();
        Info<< "=== Time-step breakdown (pisoFoam) ===" << nl
            << "  UEqn assemble+solve                 : \"" << fmt8(accUEqn)               << "\" [s]  (\"" << pct(accUEqn)               << "\" %)" << nl
            << "  pEqn total                         : \"" << fmt8(accPEqnAsm+accPEqnSol) << "\" [s]  (\"" << pct(accPEqnAsm+accPEqnSol) << "\" %)" << nl
            << "     pEqn(assemble)                  : \"" << fmt8(accPEqnAsm)            << "\" [s]  (\"" << pct(accPEqnAsm)            << "\" %)" << nl
            << "     pEqn(solve)                     : \"" << fmt8(accPEqnSol)            << "\" [s]  (\"" << pct(accPEqnSol)            << "\" %)" << nl
            << "  turbulence.correct                 : \"" << fmt8(accTurb)                << "\" [s]  (\"" << pct(accTurb)                << "\" %)" << nl
            << "  write()                            : \"" << fmt8(accWrite)               << "\" [s]  (\"" << pct(accWrite)               << "\" %)" << nl
            << "Time-step total: \"" << fmt8(stepTotal) << "\" [s]" << nl;
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
