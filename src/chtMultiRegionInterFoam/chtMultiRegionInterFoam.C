/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    chtMultiRegionInterFoam

Group
    grpHeatTransferSolvers

Description
    Transient solver for two compressible, non-isothermal immiscible fluids
    using a VOF (volume of fluid) phase-fraction based interface capturing
    approach, with solid heat conduction and conjugate heat transfer between
    the solid and fluid regions.

    The fluid region uses compressibleInterFoam physics (algebraic MULES);
    the solid regions use chtMultiRegionFoam conduction. Regions are coupled
    through compressible::turbulentTemperatureRadCoupledMixed boundaries.

    Optional adaptive outer-corrector reduction (PIMPLE adaptiveOuterTol,
    adaptiveOuterSteps; off by default): after a run of quasi-steady steps the
    solver drops to a single tight outer corrector.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"


#include "compressibleInterPhaseTransportModel.H"
#include "compressibleCourantNo.H"
#include "pimpleControl.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "fvcSmooth.H"
#include "coordinateSystem.H"
#include "loopControl.H"

#include "turbulentFluidThermoModel.H"
#include "pressureControl.H"

#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for two compressible, non-isothermal immiscible fluids "
        "using VOF phase-fraction based interface capturing, with "
        "solid heat conduction with conjugate heat transfer "
        "between solid and fluid regions."
    );

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    // #include "setRootCase.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"

    #include "createMeshes.H"

    #include "createFields.H"
    
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"

    // Local time stepping (pseudo-transient) is not implemented; the solver
    // is fully transient. LTS is kept as a compile-time switch because the
    // MULES alpha equation selects its ddt scheme from it.
    const bool LTS = false;

    if (!LTS)
    {
        #include "alphaMultiRegionCourantNo.H"
        #include "compressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setInitialMultiRegionDeltaT.H"
    }

    Info<< "\nStarting time loop\n" << endl;

    // Adaptive outer correction: was the previous step fully quasi-steady
    // (flow + thermal)? If so and the option is enabled, the next step uses a
    // single tight outer corrector. Starts false so the run always begins fully coupled.
    // A hysteresis counter requires several CONSECUTIVE quasi-steady steps
    // before reducing correctors, so an isolated low-change step during an
    // active transient cannot trigger an under-resolved (diverging) step.
    bool prevStepQuasiSteady = false;
    label steadyCount = 0;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"
        #include "readPIMPLEControls.H"

        {
            #include "alphaMultiRegionCourantNo.H"
            #include "compressibleMultiRegionCourantNo.H"
            #include "solidRegionDiffusionNo.H"
            #include "setMultiRegionDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Adaptive outer-corrector reduction ---
        // On a step that the previous step flagged as fully quasi-steady (both
        // flow and thermal change below adaptiveOuterTol), run a SINGLE outer
        // corrector. That corrector is finalIter, so the tight Final solver
        // settings apply (p_rghFinal in fluid/pEqn.H, TFinal, hFinal) and the
        // step is not left under-resolved. It skips the second corrector and
        // the energy sub-loop, halving the pressure solves on easy steps.
        // adaptiveOuterTol <= 0 (default) disables this and reproduces the
        // fixed-corrector driver exactly.
        const scalar adaptiveOuterTol =
            pimple.lookupOrDefault<scalar>("adaptiveOuterTol", 0.0);
        const label adaptiveOuterSteps =
            pimple.lookupOrDefault<label>("adaptiveOuterSteps", 5);
        const bool adaptiveOuter = (adaptiveOuterTol > SMALL);

        const int nOuterThis =
            (adaptiveOuter && prevStepQuasiSteady) ? 1 : nOuterCorr;

        if (adaptiveOuter && nOuterThis < nOuterCorr)
        {
            Info<< "OUTER-EXIT: quasi-steady, using " << nOuterThis
                << " outer corrector (tight)" << endl;
        }

        // --- PIMPLE loop
        for (int oCorr=0; oCorr<nOuterThis; ++oCorr)
        {
            const bool finalIter = (oCorr == nOuterThis-1);

            forAll(fluidRegions, i)
            {
                Info<< "\nSolving for fluid region "
                    << fluidRegions[i].name() << endl;
                #include "setRegionFluidFields.H"
                #include "readFluidMultiRegionPIMPLEControls.H"
                #include "solveFluid.H"
            }

            forAll(solidRegions, i)
            {
                Info<< "\nSolving for solid region "
                    << solidRegions[i].name() << endl;
                #include "setRegionSolidFields.H"
                #include "readSolidMultiRegionPIMPLEControls.H"
                #include "solveSolid.H"
            }

            // Additional loops for energy solution only
            if (!oCorr && nOuterThis > 1)
            {
                loopControl looping(runTime, pimple, "energyCoupling");

                while (looping.loop())
                {
                    Info<< nl << looping << nl;

                    forAll(fluidRegions, i)
                    {
                        Info<< "\nSolving for fluid region "
                            << fluidRegions[i].name() << endl;
                       #include "setRegionFluidFields.H"
                       #include "readFluidMultiRegionPIMPLEControls.H"
                       frozenFlow = true;
                       #include "solveFluid.H"
                    }

                    forAll(solidRegions, i)
                    {
                        Info<< "\nSolving for solid region "
                            << solidRegions[i].name() << endl;
                        #include "setRegionSolidFields.H"
                        #include "readSolidMultiRegionPIMPLEControls.H"
                        #include "solveSolid.H"
                    }
                }
            }
        }

        // Adaptive outer correction: measure this step's flow and thermal
        // change (relative to the start-of-step oldTime values) to decide
        // whether the NEXT step is quasi-steady. No extra field storage is needed.
        if (adaptiveOuter)
        {
            scalar flowChange = 0;
            forAll(fluidRegions, i)
            {
                const volVectorField& U = UFluid[i];
                const scalar denom =
                    gMax(mag(U.oldTime().primitiveField())) + SMALL;
                flowChange = max
                (
                    flowChange,
                    gMax(mag(U.primitiveField() - U.oldTime().primitiveField()))
                  / denom
                );
            }

            scalar thermChange = 0;
            forAll(solidRegions, i)
            {
                const volScalarField& h = thermos[i].he();
                const scalar denom =
                    gMax(mag(h.oldTime().primitiveField())) + SMALL;
                thermChange = max
                (
                    thermChange,
                    gMax(mag(h.primitiveField() - h.oldTime().primitiveField()))
                  / denom
                );
            }

            if (flowChange < adaptiveOuterTol && thermChange < adaptiveOuterTol)
            {
                ++steadyCount;
            }
            else
            {
                steadyCount = 0;
            }

            // Only treat the next step as quasi-steady after enough consecutive
            // settled steps (hysteresis), and drop straight back to full
            // coupling on the first non-settled step.
            prevStepQuasiSteady = (steadyCount >= adaptiveOuterSteps);

            Info<< "adaptiveOuter: flowChange=" << flowChange
                << " thermChange=" << thermChange
                << " steadyCount=" << steadyCount
                << " nextCorrectors=" << (prevStepQuasiSteady ? 1 : nOuterCorr)
                << endl;
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
