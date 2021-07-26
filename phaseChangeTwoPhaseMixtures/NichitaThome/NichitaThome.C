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

#include "NichitaThome.H"
#include "fvc.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(NichitaThome, 0);
    addToRunTimeSelectionTable(phaseChangeTwoPhaseMixture, NichitaThome, components);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtures::NichitaThome::NichitaThome
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    phaseChangeTwoPhaseMixture(typeName, U, phi),

    cond_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("condensation")),
    evap_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("evaporation"))
{
	Info<< "NichitaThome model settings:  " << endl;
	Info<< "Condensation is " << cond_	<< endl;
	Info<< "Evaporation is "  << evap_  << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam ::volVectorField Foam::phaseChangeTwoPhaseMixtures::NichitaThome::calcGradAlphal() const
{
	volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));
	return fvc::grad(limitedAlpha1);
}

Foam ::volVectorField Foam::phaseChangeTwoPhaseMixtures::NichitaThome::calcGradT() const
{
	return fvc::grad(T_);
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::NichitaThome::mDotAlphal() const
{
    volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));
    const dimensionedScalar T0("0", dimTemperature, 0.0);
	const volScalarField kEff = this->k();
	const volScalarField gradAlphaGradT = calcGradAlphal() & calcGradT();

	if (cond_ && evap_)
	{
		return Pair<tmp<volScalarField> >
		(
			neg(T_ - TSat_)*kEff*gradAlphaGradT/hEvap_/max(1.0 - limitedAlpha1, 1E-6),
		   -pos(T_ - TSat_)*kEff*gradAlphaGradT/hEvap_/max(limitedAlpha1, 1E-6)
		);
	}
	else if (cond_)
	{
		return Pair<tmp<volScalarField> >
		(
			neg(T_ - TSat_)*kEff*gradAlphaGradT/hEvap_,
			mEvapAlphal_*scalar(0)
		);
	}
	else if (evap_)
	{
		return Pair<tmp<volScalarField> >
		(
			mCondAlphal_*scalar(0),
	       -pos(T_ - TSat_)*kEff*gradAlphaGradT/hEvap_
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
Foam::phaseChangeTwoPhaseMixtures::NichitaThome::mDotP() const
{
    volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));
    const dimensionedScalar T0("0", dimTemperature, 0.0);
	const volScalarField kEff = this->k();
	const volScalarField gradAlphaGradT = calcGradAlphal() & calcGradT();

	if (cond_ && evap_)
	{
		return Pair<tmp<volScalarField> >
		(
			neg(T_ - TSat_)*kEff*gradAlphaGradT/hEvap_
						*pos(p_-pSat_)/max(p_-pSat_,1E-6*pSat_),
		   -pos(T_ - TSat_)*kEff*gradAlphaGradT/hEvap_
						*neg(p_-pSat_)/max(pSat_-p_,1E-6*pSat_)
		);
	}
	else if (cond_)
	{
		return Pair<tmp<volScalarField> >
		(
			neg(T_ - TSat_)*kEff*gradAlphaGradT/hEvap_
						*pos(p_-pSat_)/max(p_-pSat_,1E-6*pSat_),
			mEvapP_*scalar(0)
		);
	}
	else if (evap_)
	{
		return Pair<tmp<volScalarField> >
		(
			mCondP_*scalar(0),
		   -pos(T_ - TSat_)*kEff*gradAlphaGradT/hEvap_
						*neg(p_-pSat_)/max(pSat_-p_,1E-6*pSat_)
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
Foam::phaseChangeTwoPhaseMixtures::NichitaThome::mDotT() const
{
    volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));
	const volScalarField kEff = this->k();
	const volScalarField gradAlphaGradT = calcGradAlphal() & calcGradT();

	if (cond_ && evap_)
	{
		return Pair<tmp<volScalarField> >
		(
		   -neg(T_ - TSat_)*kEff*gradAlphaGradT/hEvap_/max(TSat_ - T_, 1E-6*TSat_),
		    pos(T_ - TSat_)*kEff*gradAlphaGradT/hEvap_/max(T_ - TSat_, 1E-6*TSat_)
		);
	}
	else if (cond_)
	{
		return Pair<tmp<volScalarField> >
		(
		   -neg(T_ - TSat_)*kEff*gradAlphaGradT/hEvap_/max(TSat_ - T_, 1E-6*TSat_),
			mEvapT_*scalar(0)
		);
	}
	else if (evap_)
	{
		return Pair<tmp<volScalarField> >
		(
			mCondT_*scalar(0),
		    pos(T_ - TSat_)*kEff*gradAlphaGradT/hEvap_/max(T_ - TSat_, 1E-6*TSat_)
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

bool Foam::phaseChangeTwoPhaseMixtures::NichitaThome::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("condensation") >> cond_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("evaporation") >> evap_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
