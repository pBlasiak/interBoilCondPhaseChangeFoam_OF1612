/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeTwoPhaseMixture, 0);
    defineRunTimeSelectionTable(phaseChangeTwoPhaseMixture, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixture::phaseChangeTwoPhaseMixture
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    thermalIncompressibleTwoPhaseMixture(U, phi),
    phaseChangeTwoPhaseMixtureCoeffs_(subDict(type + "Coeffs")),
    mCondDotAlphal_
    (
        IOobject
        (
            "mCondDotAlphal",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mCondDotAlphal", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
    ),
    mEvapDotAlphal_
    (
        IOobject
        (
            "mEvapDotAlphal",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mEvapDotAlphal", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
    ),
    mCondDotP_
    (
        IOobject
        (
            "mCondDotP",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mCondDotP", dimensionSet(0, -2, 1, 0, 0, 0, 0), 0.0)
    ),
    mEvapDotP_
    (
        IOobject
        (
            "mEvapDotP",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mEvapDotP", dimensionSet(0, -2, 1, 0, 0, 0, 0), 0.0)
    ),
    mCondDotT_
    (
        IOobject
        (
            "mCondDotT",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mCondDotT", dimensionSet(1, -3, -1, -1, 0, 0, 0), 0.0)
    ),
    mEvapDotT_
    (
        IOobject
        (
            "mEvapDotT",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mEvapDotT", dimensionSet(1, -3, -1, -1, 0, 0, 0), 0.0)
    ),
	p_(U_.db().lookupObject<volScalarField>("p")),
	T_(U_.db().lookupObject<volScalarField>("T")),
    TSatG_("TSatGlobal", dimTemperature, phaseChangeTwoPhaseMixtureCoeffs_.lookup("TSatGlobal")),
    TSat_
    (
        IOobject
        (
            "TSat",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
		TSatG_
    ),
    pSat_("pSat", dimPressure, phaseChangeTwoPhaseMixtureCoeffs_.lookup("pSat")),
    hEvap_("hEvap", dimEnergy/dimMass, phaseChangeTwoPhaseMixtureCoeffs_.lookup("hEvap")),
    R_("R", dimGasConstant, phaseChangeTwoPhaseMixtureCoeffs_.lookup("R")),
    TSatLocalPressure_(readBool(phaseChangeTwoPhaseMixtureCoeffs_.lookup("TSatLocalPressure"))),
    cutoff_("cutoff", dimless, phaseChangeTwoPhaseMixtureCoeffs_.lookup("cutoff")),
    NimaSamTEqn_("NimaSamTEqn", phaseChangeTwoPhaseMixtureCoeffs_.lookup("NimaSamTEqn")),
    NimaSamAlphaEqn_("NimaSamAlphaEqn", phaseChangeTwoPhaseMixtureCoeffs_.lookup("NimaSamAlphaEqn"))
{
	Info<< "TSatGlobal = "				<< TSatG_ << endl;
	Info<< "pSat = "		  			<< pSat_ << endl;
	Info<< "hEvap = "		  			<< hEvap_ << endl;
	Info<< "R = "			  			<< R_ << endl;
	Info<< "TSatLocalPressure = "       << TSatLocalPressure_ << endl;
	Info<< "cutoff = "					<< cutoff_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseChangeTwoPhaseMixture::calcTSatLocal() 
{
	if (TSatLocalPressure_)
	{
		Info <<"TSat is calculated based on local pressure field." << endl;
	    TSat_ = 1.0/(1.0/TSatG_ - R_/hEvap_*log(max(p_/pSat_,1E-08)));
	}
	else
	{
		Info <<"TSat is constant, TSat = " << TSatG_ << endl;
	}
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixture::vDotAlphal()
//Foam::phaseChangeTwoPhaseMixture::vDotAlphal() const
{
    volScalarField alphalCoeff(1.0/rho1() - alpha1_*(1.0/rho1() - 1.0/rho2()));
    Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal();

    return Pair<tmp<volScalarField>>
    (
        alphalCoeff*mDotAlphal[0],
        alphalCoeff*mDotAlphal[1]
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixture::vDotT()
//Foam::phaseChangeTwoPhaseMixture::vDotT() const
{
	   if (NimaSamTEqn_)
	   {
	   volScalarField rhoCp =  rho1()*cp1()*alpha1_ + rho2()*cp2()*(1.0-alpha1_);
	   volScalarField TCoeff = hEvap_/rhoCp;
	   Pair<tmp<volScalarField> > mDotT = this->mDotT();

	    return Pair<tmp<volScalarField> >
	    (
	    		TCoeff*mDotT[0],
	    		TCoeff*mDotT[1]
	    );
	   }
	   else
	   {
	   Pair<tmp<volScalarField> > mDotT = this->mDotT();

	    return Pair<tmp<volScalarField> >
	    (
	    		hEvap_*mDotT[0],
	    		hEvap_*mDotT[1]
	    );
	   }
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeTwoPhaseMixture::vDotP()
//Foam::phaseChangeTwoPhaseMixture::vDotP() const
{
    dimensionedScalar pCoeff(1.0/rho1() - 1.0/rho2());
    Pair<tmp<volScalarField>> mDotP = this->mDotP();

    return Pair<tmp<volScalarField>>(pCoeff*mDotP[0], pCoeff*mDotP[1]);
}

void Foam::phaseChangeTwoPhaseMixture::correct()
{
	calcTSatLocal();
}

bool Foam::phaseChangeTwoPhaseMixture::read()
{
    if (incompressibleTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("TSatGlobal") >> TSatG_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("TSatLocalPressure") >> TSatLocalPressure_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("pSat") >> pSat_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("hEvap") >> hEvap_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("R") >> R_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("cutoff") >> cutoff_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
