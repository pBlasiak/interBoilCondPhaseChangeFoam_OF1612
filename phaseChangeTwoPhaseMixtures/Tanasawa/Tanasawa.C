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
	// czy nie powinno byc dla wrzenia mnozone przez rho1?
	// czyli czy nie powinno byc mvCoeff i mlCoeff?
   	mCoeff_(2.0*gamma_/(2.0 - gamma_)/sqrt(2.0*M_PI*R_)*hEvap_*rho2())
{
	Info<< "Tanasawa model settings:  " << endl;
	Info<< "Condensation is " << cond_	<< endl;
	Info<< "Evaporation is "  << evap_  << endl;
	Info<< "gamma = "		  << gamma_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam ::volScalarField Foam::phaseChangeTwoPhaseMixtures::Tanasawa::calcGradAlphal() const
{
	volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));
	return mag(fvc::grad(limitedAlpha1));
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::mDotAlphal() const
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
			mEvapAlphal_*scalar(0)
		);
	}
	else if (evap_)
	{
		return Pair<tmp<volScalarField> >
		(
			mCondAlphal_*scalar(0),
			-mCoeff_*max(T_ - TSat_ ,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0))
		);
	}
	else 
	{
		return Pair<tmp<volScalarField> >
		(
			mCondAlphal_*scalar(0),
			mEvapAlphal_*scalar(0)
		);
	}
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::mDotP() const
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
			mEvapP_*scalar(0)
		);
	}
	else if (evap_)
	{
		return Pair<tmp<volScalarField> >
		(
			mCondP_*scalar(0),
			-mCoeff_*max(T_ - TSat_,T0)*calcGradAlphal()/sqrt(pow(TSat_,3.0))
						*neg(p_-pSat_)/max(pSat_-p_,1E-05*pSat_)*limitedAlpha1
		);
	}
	else 
	{
		return Pair<tmp<volScalarField> >
		(
			mCondP_*scalar(0),
			mEvapP_*scalar(0)
		);
	}
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Tanasawa::mDotT() const
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
			mEvapT_*scalar(0)
		);
	}
	else if (evap_)
	{
		return Pair<tmp<volScalarField> >
		(
			mCondT_*scalar(0),
			mCoeff_*calcGradAlphal()/sqrt(pow(TSat_,3.0))
						*limitedAlpha1*pos(T_ - TSat_)
		);
	}
	else
	{
		return Pair<tmp<volScalarField> >
		(
			mCondT_*scalar(0),
			mEvapT_*scalar(0)
		);
	}
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
