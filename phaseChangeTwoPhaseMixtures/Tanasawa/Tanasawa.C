/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "Tanasawa.H"
#include "fvc.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(Tanasawa, 0);
    addToRunTimeSelectionTable(phaseChangeTwoPhaseMixture, Tanasawa, components);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtures::Tanasawa::Tanasawa
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    phaseChangeTwoPhaseMixture(typeName, U, phi),

    cond_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("condensation")),
    evap_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("evaporation")),
    gamma_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("gamma")),
   	mCoeff_(2.0*gamma_/(2.0 - gamma_)/sqrt(2.0*M_PI*R_)*hEvap_*rho2())
{
	Info<< "Tanasawa model settings:  " << endl;
	Info<< "Condensation is " << cond_	<< endl;
	Info<< "Evaporation is "  << evap_  << endl;
	Info<< "gamma = "		  << gamma_ << endl;
	Info<< "2.0*gamma_/(2.0 - gamma_)/sqrt(2.0*M_PI*R_)*hEvap_*rho2() = " << mCoeff_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam ::volScalarField Foam::phaseChangeTwoPhaseMixtures::Tanasawa::calcGradAlphal() const
{
	volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));
	return mag(fvc::grad(limitedAlpha1));
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::mDotAlphal()
//Foam::phaseChangeTwoPhaseMixtures::Tanasawa::mDotAlphal() const
{
    const dimensionedScalar T0("0", dimTemperature, 0.0);

	if (cond_ && evap_)
	{
		return Pair<tmp<volScalarField> >
		(
			-mCoeff_*min(T_ - TSat_ ,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0)),
			-mCoeff_*max(T_ - TSat_ ,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0)) 
		);
	}
	else if (cond_)
	{
		return Pair<tmp<volScalarField> >
		(
			-mCoeff_*min(T_ - TSat_ ,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0)),
			mEvapDotAlphal_*scalar(0)
		//	-mCoeff_*max(T_ - TSat_ ,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0))*scalar(0)
		);
	}
	else if (evap_)
	{
		return Pair<tmp<volScalarField> >
		(
			//-mCoeff_*min(T_ - TSat_ ,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0))*scalar(0),
			mCondDotAlphal_*scalar(0),
			-mCoeff_*max(T_ - TSat_ ,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0))
		);
	}
	else 
	{
		return Pair<tmp<volScalarField> >
		(
			mCondDotAlphal_*scalar(0),
			mEvapDotAlphal_*scalar(0)
			//-mCoeff_*min(T_ - TSat_ ,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0))*scalar(0),
			//-mCoeff_*max(T_ - TSat_ ,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0))*scalar(0)
		);
	}
	//if (cond_)
	//{
	//	mCondDotAlphal_ = -mCoeff_*min(T_ - TSat_ ,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0));
	//}
	//else 
	//{
    //    const dimensionedScalar mCondDotAlphal0
	//							(
	//								"mCondDotAlphal0", 
	//								dimensionSet(1, -3, -1, 0, 0, 0, 0),
	//								0.0
	//							);
	//	mCondDotAlphal_ = mCondDotAlphal0;
	//}

	//if (evap_)
	//{
	//	mEvapDotAlphal_ = -mCoeff_*max(T_ - TSat_ ,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0));
	//}
	//else 
	//{
    //    const dimensionedScalar mEvapDotAlphal0("mEvapDotAlphal0", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0);
	//	mEvapDotAlphal_ = mEvapDotAlphal0;
	//}

	//tmp<volScalarField> tmCondDotAlphal
	//(
	//	new volScalarField
	//	(
	//		IOobject
    //    	(
    //    	    U_.time().timeName(),
    //    	    U_.db(),
	//			IOobject::NO_READ,
	//			IOobject::NO_WRITE
    //    	),
	//		mCondDotAlphal_ 
	//	)
	//);

	//tmp<volScalarField> tmEvapDotAlphal
	//(
	//	new volScalarField
	//	(
	//		IOobject
    //    	(
    //    	    U_.time().timeName(),
    //    	    U_.db(),
	//			IOobject::NO_READ,
	//			IOobject::NO_WRITE
    //    	),
	//		mEvapDotAlphal_
	//	)
	//);

	//return Pair<tmp<volScalarField> >
	//(
	//	tmCondDotAlphal,
	//	tmEvapDotAlphal
	//);

    //return Pair<tmp<volScalarField> >
    //(
    //	-rc_*Cm1_*min(T_ - TSat_ ,T0)*AbyV()/sqrt(pow(TSat_,3.0)),

    //    -rv_*Cm1_*max(T_ - TSat_ ,T0)*AbyV()/sqrt(pow(TSat_,3.0))
    //);

}
Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::mDotP()
//Foam::phaseChangeTwoPhaseMixtures::Tanasawa::mDotP() const
{
    volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));
    const dimensionedScalar T0("0", dimTemperature, 0.0);

	if (cond_ && evap_)
	{
		return Pair<tmp<volScalarField> >
		(
			-mCoeff_*min(T_ - TSat_,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0))
						*pos(p_-pSat_)/max(p_-pSat_,1E-6*pSat_)*(1.0-limitedAlpha1),
			-mCoeff_*max(T_ - TSat_,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0))
						*neg(p_-pSat_)/max(pSat_-p_,1E-05*pSat_)*limitedAlpha1 
		);
	}
	else if (cond_)
	{
		return Pair<tmp<volScalarField> >
		(
			-mCoeff_*min(T_ - TSat_,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0))
						*pos(p_-pSat_)/max(p_-pSat_,1E-6*pSat_)*(1.0-limitedAlpha1),
			mEvapDotP_*scalar(0)
		//	-mCoeff_*max(T_ - TSat_,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0))
		//				*neg(p_-pSat_)/max(pSat_-p_,1E-05*pSat_)*limitedAlpha1*scalar(0)
		);
	}
	else if (evap_)
	{
		return Pair<tmp<volScalarField> >
		(
			//-mCoeff_*min(T_ - TSat_,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0))
			//			*pos(p_-pSat_)/max(p_-pSat_,1E-6*pSat_)*(1.0-limitedAlpha1)*scalar(0),
			mCondDotP_*scalar(0),
			-mCoeff_*max(T_ - TSat_,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0))
						*neg(p_-pSat_)/max(pSat_-p_,1E-05*pSat_)*limitedAlpha1
		);
	}
	else 
	{
		return Pair<tmp<volScalarField> >
		(
			mCondDotP_*scalar(0),
			mEvapDotP_*scalar(0)
		//	-mCoeff_*min(T_ - TSat_,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0))
		//				*pos(p_-pSat_)/max(p_-pSat_,1E-6*pSat_)*(1.0-limitedAlpha1)*scalar(0),
		//	-mCoeff_*max(T_ - TSat_,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0))
		//				*neg(p_-pSat_)/max(pSat_-p_,1E-05*pSat_)*limitedAlpha1*scalar(0)
		);
	}


	//if (cond_)
	//{
	//	mCondDotP_ = -mCoeff_*min(T_ - TSat_,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0))
	//					*pos(p_-pSat_)/max(p_-pSat_,1E-6*pSat_)*(1.0-limitedAlpha1);
	//}
	//else 
	//{
    //    const dimensionedScalar mCondDotP0("mCondDotP0", dimensionSet(0, -2, 1, 0, 0, 0, 0), 0.0);
	//	mCondDotP_ = mCondDotP0; 
	//}

	//if (evap_)
	//{
	//	mEvapDotP_ = -mCoeff_*max(T_ - TSat_,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0))
	//					*neg(p_-pSat_)/max(pSat_-p_,1E-05*pSat_)*limitedAlpha1;
	//}
	//else 
	//{
    //    const dimensionedScalar mEvapDotP0("mEvapDotP0", dimensionSet(0, -2, 1, 0, 0, 0, 0), 0.0);
	//	mEvapDotP_ = mEvapDotP0;
	//}

	//tmp<volScalarField> tmCondDotP
	//(
	//	new volScalarField
	//	(
	//		IOobject
    //    	(
    //    	    U_.time().timeName(),
    //    	    U_.db(),
	//			IOobject::NO_READ,
	//			IOobject::NO_WRITE
    //    	),
	//		mCondDotP_
	//	)
	//);

	//tmp<volScalarField> tmEvapDotP
	//(
	//	new volScalarField
	//	(
	//		IOobject
    //    	(
    //    	    U_.time().timeName(),
    //    	    U_.db(),
	//			IOobject::NO_READ,
	//			IOobject::NO_WRITE
    //    	),
	//		mEvapDotP_
	//	)
	//);
	//return Pair<tmp<volScalarField> >
	//(
	//	tmCondDotP, 
	//	tmEvapDotP
	//);
    //return Pair<tmp<volScalarField> >
    //(
    //		-rc_*Cm1_*min(T_ - TSat_,T0)*AbyV()/sqrt(pow(TSat_,3.0))
    //		*pos(p_-pSat_)/max(p_-pSat_,1E-6*pSat_)
    //		*(1.0-limitedAlpha1),

    //		-rv_*Cm1_*max(T_ - TSat_,T0)*AbyV()/sqrt(pow(TSat_,3.0))
    //		*neg(p_-pSat_)/max(pSat_-p_,1E-05*pSat_)
    //		*limitedAlpha1
    //);
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::mDotT() 
//Foam::phaseChangeTwoPhaseMixtures::Tanasawa::mDotT() const
{
    volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));

	if (cond_ && evap_)
	{
		return Pair<tmp<volScalarField> >
		(
			-mCoeff_*calcGradAlphal()/sqrt(pow(TSat_,3.0))
						*neg(T_ - TSat_)*(1.0-limitedAlpha1),
			mCoeff_*calcGradAlphal()/sqrt(pow(TSat_,3.0))
						*limitedAlpha1*pos(T_ - TSat_) 
		);
	}
	else if (cond_)
	{
		return Pair<tmp<volScalarField> >
		(
			-mCoeff_*calcGradAlphal()/sqrt(pow(TSat_,3.0))
						*neg(T_ - TSat_)*(1.0-limitedAlpha1),
			mEvapDotT_*scalar(0)
			//mCoeff_*calcGradAlphal()/sqrt(pow(TSat_,3.0))
			//			*limitedAlpha1*pos(T_ - TSat_)*scalar(0)
		);
	}
	else if (evap_)
	{
		return Pair<tmp<volScalarField> >
		(
		//	-mCoeff_*calcGradAlphal()/sqrt(pow(TSat_,3.0))
		//				*neg(T_ - TSat_)*(1.0-limitedAlpha1)*scalar(0),
			mCondDotT_*scalar(0),
			mCoeff_*calcGradAlphal()/sqrt(pow(TSat_,3.0))
						*limitedAlpha1*pos(T_ - TSat_)
		);
	}
	else
	{
		return Pair<tmp<volScalarField> >
		(
			mCondDotT_*scalar(0),
			mEvapDotT_*scalar(0)
		//	-mCoeff_*calcGradAlphal()/sqrt(pow(TSat_,3.0))
		//				*neg(T_ - TSat_)*(1.0-limitedAlpha1)*scalar(0),
		//	mCoeff_*calcGradAlphal()/sqrt(pow(TSat_,3.0))
		//				*limitedAlpha1*pos(T_ - TSat_)*scalar(0)
		);
	}
	//if (cond_)
	//{
	//	mCondDotT_ = -mCoeff_*calcGradAlphal()/sqrt(pow(TSat_,3.0))
	//					*neg(T_ - TSat_)*(1.0-limitedAlpha1);
	//}
	//else 
	//{
    //    const dimensionedScalar mCondDotT0("mCondDotT0", dimensionSet(1, -3, -1, -1, 0, 0, 0), 0.0);
	//	mCondDotT_ = mCondDotT0;
	//}

	//if (evap_)
	//{
	//	mEvapDotT_ = mCoeff_*calcGradAlphal()/sqrt(pow(TSat_,3.0))
	//					*limitedAlpha1*pos(T_ - TSat_);
	//}
	//else 
	//{
    //    const dimensionedScalar mEvapDotT0("mEvapDotT0", dimensionSet(1, -3, -1, -1, 0, 0, 0), 0.0);
	//	mEvapDotT_ = mEvapDotT0;
	//}

	//tmp<volScalarField> tmCondDotT
	//(
	//	new volScalarField
	//	(
	//		IOobject
    //    	(
    //    	    U_.time().timeName(),
    //    	    U_.db(),
	//			IOobject::NO_READ,
	//			IOobject::NO_WRITE
    //    	),
	//		mCondDotT_
	//	)
	//);

	//tmp<volScalarField> tmEvapDotT
	//(
	//	new volScalarField
	//	(
	//		IOobject
    //    	(
    //    	    U_.time().timeName(),
    //    	    U_.db(),
	//			IOobject::NO_READ,
	//			IOobject::NO_WRITE
    //    	),
	//		mEvapDotT_
	//	)
	//);
	//return Pair<tmp<volScalarField> >
	//(
	//	tmCondDotT, 
	//	tmEvapDotT
	//);
    //return Pair<tmp<volScalarField> >
    //(
    //       -rc_*Cm1_*AbyV()/sqrt(pow(TSat_,3.0))
    //        *neg(T_ - TSat_)
    //        *(1.0-limitedAlpha1),

    //        rv_*Cm1_*AbyV()/sqrt(pow(TSat_,3.0))
    //       *limitedAlpha1
    //       *pos(T_ - TSat_)
    //);
}

bool Foam::phaseChangeTwoPhaseMixtures::Tanasawa::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("condensation") >> cond_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("evaporation") >> evap_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("gamma") >> gamma_;

		mCoeff_ = 2.0*gamma_/(2.0 - gamma_)/sqrt(2.0*M_PI*R_)*hEvap_*rho2();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
