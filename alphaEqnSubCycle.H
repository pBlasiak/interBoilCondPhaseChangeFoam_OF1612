{
    // Standard face-flux compression coefficient
    surfaceScalarField phic(interface.cAlpha()*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (interface.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    volScalarField divU(fvc::div(phi));

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        surfaceScalarField rhoPhiSum("rhoPhiSum", rhoPhi);
        surfaceScalarField rhoPhiCpSum("rhoPhiCpSum", rhoPhiCp);

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqn.H"
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
            rhoPhiCpSum += (runTime.deltaT()/totalDeltaT)*rhoPhiCp;
        }

        rhoPhi = rhoPhiSum;
        rhoPhiCp = rhoPhiCpSum;
    }
    else
    {
        #include "alphaEqn.H"
    }

    rho == alpha1*rho1 + alpha2*rho2;
    rhoCp == alpha1*rho1*cp1 + alpha2*rho2*cp2;
}
