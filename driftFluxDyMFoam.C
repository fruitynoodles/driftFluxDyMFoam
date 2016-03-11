/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    driftFluxDyMFoam

Description
    Solver for 2 incompressible fluids using the mixture approach with the
    drift-flux approximation for relative motion of the phases with optional dynamic mesh
    stuff cobbled on.

    Used for simulating the settling of the dispersed phase when a mixer is present,
    and other similar separation problems.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "subCycle.H"
#include "incompressibleTwoPhaseInteractingMixture.H"
#include "relativeVelocityModel.H"
#include "turbulenceModel.H"
#include "CompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "gaussLaplacianScheme.H"
#include "uncorrectedSnGrad.H"
#include "CorrectPhi.H"
//#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"

    pimpleControl pimple(mesh);

    #include "createControls.H"
    #include "createFields.H"
    #include "createMRF.H"
    #include "createFvOptions.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    volScalarField rAU
    (
         IOobject
         (
             "rAU",
             runTime.timeName(),
             mesh,
             IOobject::READ_IF_PRESENT,
             IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1.0)
    );

    #include "correctPhi.H"
    #include "createUf.H"
    #include "readControls.H"
    #include "CourantNo.H"
    #include "setDeltaT.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();
            mesh.update();
            if(mesh.changing())
            {
                Info<< "Executiontime for mesh.update() = "
                    << runTime.elapsedCpuTime()-timeBeforeMeshUpdate
                    << " s" <<endl;
                gh = (g & mesh.C())-ghRef;
                ghf = (g & mesh.Cf())-ghRef;

            }
            if(mesh.changing() && correctPhi)
            {
                phi = mesh.Sf() & Uf;
                #include "correctPhi.H"
                fvc::makeRelative(phi,U);
                mixture.correct();
                #include "alphaControls.H"
                UdmModel.correct();
            }
            if(mesh.changing() && checkMeshCourantNo)
            {
                #include "meshCourantNo.H"
            }

            #include "alphaControls.H"

            UdmModel.correct();

            #include "alphaEqnSubCycle.H"


            #include "UEqn.H"

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

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
