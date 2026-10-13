// from compressibleInterFoam (unchanged)


if (finalIter)
{
    mesh.data().setFinalIteration(true);
}


if (frozenFlow)
{
    // Energy-only update: flow fields are held fixed (used by the
    // energy-coupling sub-loop to tighten conjugate heat transfer)
    #include "contErr.H"
    volScalarField divUp("divUp", fvc::div(fvc::absolute(phi, U), p));
    #include "TEqn.H"
}
else
{
    #include "alphaControls.H"
    #include "compressibleAlphaEqnSubCycle.H"

    turbulence.correctPhasePhi();

    #include "UEqn.H"
    volScalarField divUp("divUp", fvc::div(fvc::absolute(phi, U), p));
    // Temperature form of the energy equation (no radiation source),
    // as in compressibleInterFoam
    #include "TEqn.H"

    // --- Pressure corrector loop
    // --- PISO loop
    for (int corr=0; corr<nCorr; corr++)
    {
        #include "pEqn.H"
    }

    if((turbOnFinalIterOnly && finalIter) || !turbOnFinalIterOnly)
    {
        turbulence.correct();
    }
}

if (finalIter)
{
    mesh.data().setFinalIteration(false);
}
