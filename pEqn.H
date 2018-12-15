{
	if (!NimaSamPEqn) 
	{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );
    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p_rgh);

    surfaceScalarField phig
    (
        (
            interface.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
//    constrainPressure(p_rgh, U, phiHbyA, rAUf);
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    Pair<tmp<volScalarField>> vDotP = mixture->vDotP();
    const volScalarField& vDotcP = vDotP[0]();
    const volScalarField& vDotvP = vDotP[1]();

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvc::div(phiHbyA) - fvm::laplacian(rAUf, p_rgh)
		==
			fvm::Sp(vDotcP-vDotvP,p_rgh) - (vDotcP-vDotvP)*(pSat-rho*gh) 
        );

        p_rghEqn.setReference(pRefCell, pRefValue);

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + p_rghEqn.flux();

            //p_rgh.relax(); // from interFoam

            U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
	}
	else // NIMA SAM pEqn
	{
	Info << "Nima Sam pEqn is solved." << endl;
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("Dp", fvc::interpolate(rAU));

   volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
      // + fvc::ddtPhiCorr(rAU, rho, U, phi)
    );
    adjustPhi(phiHbyA, U, p_rgh);
    phi = phiHbyA;

    surfaceScalarField phig
    (
        (
           fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          //  fvc::interpolate(interface.sigmaK())*fvc::snGrad(interface.smoothAlpha())
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;
/*
    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
      (
          p_rgh.boundaryField(),
          (
              phiHbyA.boundaryField()
            - (mesh.Sf().boundaryField() & U.boundaryField())
          )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
      );
*/
    Pair<tmp<volScalarField> > vDotP = mixture->vDotP();
    const volScalarField& vDotcP = vDotP[0]();
    const volScalarField& vDotvP = vDotP[1]();

//#include "HelmholtzEqn.H"



    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
           fvc::div(phiHbyA) - fvm::laplacian(rAUf, p_rgh)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));


        		  solve
        		  (
        				  p_rghEqn == fvm::Sp(vDotcP-vDotvP,p_rgh) - (vDotcP-vDotvP)*(pSat-rho*gh) ,
        				  //p_rghEqn == vDotcP-vDotvP
        				 // p_rghEqn == fvm::Sp(vDotPs,p_rgh) - (vDotPs)*(pSat-rho*gh),
        				  mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        		  );

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + p_rghEqn.flux();

            U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    //#include "continuityErrs.H"

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

	}
    Info<< "min/max(U) = " << min(U).value() << ", "
        << max(U).value() <<endl;
}
