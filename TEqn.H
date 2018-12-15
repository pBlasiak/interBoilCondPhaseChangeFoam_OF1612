{
    surfaceScalarField alpha1f =
        min(max(fvc::interpolate(alpha1), scalar(0)), scalar(1));

    surfaceScalarField kf =
            alpha1f*k1 + (scalar(1) - alpha1f)*k2;
    surfaceScalarField cpf =
            alpha1f*cp1 + (scalar(1) - alpha1f)*cp2;
	cp = alpha1*cp1 + (1.0 - alpha1)*cp2;
	alphaEff = alpha1*k1 + (1.0 - alpha1)*k2;

	Pair<tmp<volScalarField> > vDotT = mixture->vDotT();
		const volScalarField& vDotcT = vDotT[0]();
		const volScalarField& vDotvT = vDotT[1]();
	    const volScalarField vDotcmvT(vDotcT - vDotvT);

    volScalarField IF
    (
        IOobject
        (
            "IF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		mesh,
        dimensionedScalar("IF", dimensionSet(1,-1,-3,-1,0,0,0), scalar(0))
    );

	//	TODO it is not used now!!!!!!!!!!
//	IF = twoPhaseProperties->IF();
	forAll(alpha1, celli)
	{
		if (alpha1[celli] < 1.0 - cutoff.value() && alpha1[celli] > cutoff.value()) IF[celli] = 1.0;
	}


	TSat = mixture->TSat();

	if (NimaSamTEqn)
	{
		Info<< "Nima Sam TEqn is solved." << endl;
	//	for Nima Sam TEqn
	volScalarField k= k1*alpha1 + (1.0 - alpha1)*k2;
	volScalarField rhoC = rho1*cp1*alpha1 + (1.0 - alpha1)*rho2*cp2;
	    fvScalarMatrix TEqn
	    (
			  fvm::ddt(T)
			+ fvm::div(phi, T)
		    - fvm::Sp(fvc::div(phi),T)
			- fvm::laplacian( k/rhoC , T,"laplacian(alphaEff,T)")
		);
    TEqn.solve();
	}
	else
	{
    fvScalarMatrix TEqn
    (
//	 // Nima Sam TEqn
//		  fvm::ddt(T)
//		+ fvm::div(phi, T)
//	    - fvm::Sp(fvc::div(phi),T)
//		- fvm::laplacian( k/rhoC , T,"laplacian(alphaEff,T)")
//	 // End of Nima Sam TEqn

        fvm::ddt(rhoCp, T)
      + fvm::div(rhoPhiCp, T)
      - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoPhiCp), T)
      - fvm::laplacian(kf, T)
     ==
		fvm::Sp(vDotcmvT, T) - vDotcmvT*TSat	// TODO update TSat in correct
    );
    TEqn.solve();
	}	
//	// OLD EQUATION:
//	// TODO try TEqn by Nima Sam?
//    fvScalarMatrix TEqn
//    (
//        fvm::ddt(rho, T)
//      + fvm::div(rhoPhi, T)
//      - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), T)
//      - fvm::laplacian(kf/cpf, T)
//     ==
//	    hEvap/cp*(fvm::Sp(mDotTv/TSat, T) + mDotTc - mDotTv)  // tak jest jak m = mc + mv
////	    hEvap/cp*(mDotTc - mDotTv*(T - TSat)/TSat)  // tak jest jak m = mc - mv
////	  - 1./cp*fvm::Sp(pow(10,10)*IF, T) + 1./cp*pow(10,10)*TSat*IF
//    );


//    TEqn.relax();
//    TEqn.solve();

    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;
}
