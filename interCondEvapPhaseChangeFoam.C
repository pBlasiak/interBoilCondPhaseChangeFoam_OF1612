/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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
    interCondEvapPhaseChangeFoam

Group
    grpMultiphaseSolvers

Description
    Solver for 2 incompressible, immiscible fluids with phase-change
    (condensation, evaporation).  Uses a VOF (volume of fluid) phase-fraction based
    interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a
    single momentum equation is solved.

    The set of phase-change models provided are designed to simulate evaporation
    and condensation.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "subCycle.H"
#include "smoothInterfaceProperties.H"
#include "phaseChangeTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

//    #include "readGravitationalAcceleration.H"
//    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createMRF.H"
    #include "createFvOptions.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "createTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
            surfaceScalarField rhoPhiCp
            (
                IOobject
                (
                    "rhoPhiCp",
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("0", dimMass/dimTime*dimSpecificHeatCapacity, 0)
            );

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaControls.H"

			dimensionedScalar zeroMassFlux("0", dimensionSet(1,2,-3,-1,0,0,0), 0.0);
			rhoPhiCp = zeroMassFlux;
            //surfaceScalarField rhoPhiCp
            //(
            //    IOobject
            //    (
            //        "rhoPhiCp",
            //        runTime.timeName(),
            //        mesh
            //    ),
            //    mesh,
            //    dimensionedScalar("0", dimMass/dimTime*dimSpecificHeatCapacity, 0)
            //);
            surfaceScalarField rhoPhi
            (
                IOobject
                (
                    "rhoPhi",
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("0", dimMass/dimTime, 0)
            );
//			dimensionedScalar zeroMassFlux("0", dimMass/dimTime, 0);
//			rhoPhi = zeroMassFlux;


            #include "alphaEqnSubCycle.H"
            interface.correct();

            #include "UEqn.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }
        //#include "TEqn.H" // w Nima Sam tutaj jest TEqn
            //mixture->correct(); // w Nima Sam tutaj jest correct()

        mixture->correct();

        runTime.write();
    volScalarField limitedAlpha1 = min(max(alpha1, scalar(0)), scalar(1));

		// TODO zrob save mv i mc zeby nie trzeba tego liczyc ???
        Info<< "****Condensation rate: "
            << gSum(mixture->mDotAlphal()[0]*mesh.V()*(1.0 - limitedAlpha1))*hEvap.value() << " W" << endl;
        Info<< "****Evaporation rate: "
            << gSum(mixture->mDotAlphal()[1]*mesh.V()*limitedAlpha1)*hEvap.value()  << " W" << endl;
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
