{
	Pair<tmp<volScalarField> > vDotT = mixture->vDotT();
		const volScalarField& vDotcT = vDotT[0]();
		const volScalarField& vDotvT = vDotT[1]();
	    const volScalarField vDotcmvT(vDotcT - vDotvT);

	TSat = mixture->TSat();

	volScalarField k = mixture->k();
	    fvScalarMatrix TEqn
	    (
			  fvm::ddt(T)
			+ fvm::div(phi, T)
	//	    - fvm::Sp(fvc::div(phi),T)
			- fvm::laplacian( k/rhoCp , T,"laplacian(alphaEff,T)")
		);
    TEqn.relax();
	solve
	(
		TEqn
     ==
		//fvm::Sp(vDotcmvT/rhoC, T) - vDotcmvT/rhoC*TSat	// TODO update TSat in correct
		fvm::Sp(vDotcmvT, T) - vDotcmvT*TSat	// TODO update TSat in correct
	);

	

    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;
}
