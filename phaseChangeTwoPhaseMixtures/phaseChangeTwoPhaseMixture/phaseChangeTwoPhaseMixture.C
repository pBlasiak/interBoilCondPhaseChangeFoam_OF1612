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
    mCond_
    (
        IOobject
        (
            "mCond",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mCond", dimensionSet(1, 0, -1, 0, 0, 0, 0), 0.0)
    ),
    mEvap_
    (
        IOobject
        (
            "mEvap",
            U_.time().timeName(),
            U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar("mEvap", dimensionSet(1, 0, -1, 0, 0, 0, 0), 0.0)
    ),
	mCondAlphal_
	(
	    IOobject
	    (
	        "mCondAlphal",
	        U_.time().timeName(),
	        U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U_.mesh(),
	    dimensionedScalar("mCondAlphal", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
	),
	mEvapAlphal_
	(
	    IOobject
	    (
	        "mEvapAlphal",
	        U_.time().timeName(),
	        U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U_.mesh(),
	    dimensionedScalar("mEvapAlphal", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0)
	),
	mCondP_
	(
	    IOobject
	    (
	        "mCondP",
	        U_.time().timeName(),
	        U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U_.mesh(),
	    dimensionedScalar("mCondP", dimensionSet(0, -2, 1, 0, 0, 0, 0), 0.0)
	),
	mEvapP_
	(
	    IOobject
	    (
	        "mEvapP",
	        U_.time().timeName(),
	        U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U_.mesh(),
	    dimensionedScalar("mEvapP", dimensionSet(0, -2, 1, 0, 0, 0, 0), 0.0)
	),
	mCondT_
	(
	    IOobject
	    (
	        "mCondT",
	        U_.time().timeName(),
	        U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U_.mesh(),
	    dimensionedScalar("mCondT", dimensionSet(1, -3, -1, -1, 0, 0, 0), 0.0)
	),
	mEvapT_
	(
	    IOobject
	    (
	        "mEvapT",
	        U_.time().timeName(),
	        U_.db(),
			IOobject::NO_READ,
			IOobject::NO_WRITE
	    ),
	    U_.mesh(),
	    dimensionedScalar("mEvapT", dimensionSet(1, -3, -1, -1, 0, 0, 0), 0.0)
	),
	p_(U_.db().lookupObject<volScalarField>("p")),
	T_(U_.db().lookupObject<volScalarField>("T")),
    TSatG_("TSatGlobal", dimTemperature, lookup("TSatGlobal")),
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
    pSat_("pSat", dimPressure, lookup("pSat")),
    hEvap_("hEvap", dimEnergy/dimMass, lookup("hEvap")),
    R_("R", dimGasConstant, lookup("R")),
    TSatLocalPressure_(readBool(lookup("TSatLocalPressure"))),
	printPhaseChange_(readBool(lookup("printPhaseChange")))
{
	Info<< "TSatGlobal = "				<< TSatG_ << endl;
	Info<< "pSat = "		  			<< pSat_ << endl;
	Info<< "hEvap = "		  			<< hEvap_ << endl;
	Info<< "R = "			  			<< R_ << endl;
	Info<< "TSatLocalPressure = "       << TSatLocalPressure_ << endl;
	Info<< "printPhaseChange = "        << printPhaseChange_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseChangeTwoPhaseMixture::calcTSatLocal() 
{
	if (TSatLocalPressure_)
	{
		//Info <<"TSat is calculated based on local pressure field." << endl;
	    TSat_ = 1.0/(1.0/TSatG_ - R_/hEvap_*log(max(p_/pSat_,1E-08)));
	}
	//else
	//{
	//	Info <<"TSat is constant, TSat = " << TSatG_ << endl;
	//}
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixture::vDotAlphal() const
{
    volScalarField alphalCoeff(1.0/rho1() - alpha1_*(1.0/rho1() - 1.0/rho2()));
    Pair<tmp<volScalarField> > mDotAlphal = this->mDotAlphal();

    return Pair<tmp<volScalarField> >
    (
        alphalCoeff*mDotAlphal[0],
        alphalCoeff*mDotAlphal[1]
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixture::vDotT() const
{
	   Pair<tmp<volScalarField> > mDotT = this->mDotT();

	   return Pair<tmp<volScalarField> >
	   (
	   		hEvap_*mDotT[0],
	   		hEvap_*mDotT[1]
	   );
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixture::vDotP() const
{
    dimensionedScalar pCoeff(1.0/rho1() - 1.0/rho2());
    Pair<tmp<volScalarField> > mDotP = this->mDotP();

    return Pair<tmp<volScalarField> >(pCoeff*mDotP[0], pCoeff*mDotP[1]);
}

void Foam::phaseChangeTwoPhaseMixture::correct()
{
	calcTSatLocal();

    volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));
    const fvMesh& mesh = alpha1_.mesh();
	// prawdopodobnie jest zle liczone mCond_ i mEvap
	// raczej powinno byc policzone bez internalField 
	// ewentualnie dodac boundaryField
    mCond_.ref() = mDotAlphal()[0]().internalField()*(1.0 - limitedAlpha1.internalField())*mesh.V();
    mEvap_.ref() = mDotAlphal()[1]().internalField()*limitedAlpha1.internalField()*mesh.V();

	if (printPhaseChange_)
	{
    	Info<< "****Condensation rate: "
    	    << gSum(mCond())*hEvap_.value() << " W" << endl;
    	Info<< "****Evaporation rate: "
    	    << gSum(mEvap())*hEvap_.value() << " W" << endl;
	}
}

bool Foam::phaseChangeTwoPhaseMixture::read()
{
    if (incompressibleTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");
        lookup("TSatGlobal") >> TSatG_;
        lookup("TSatLocalPressure") >> TSatLocalPressure_;
        lookup("pSat") >> pSat_;
        lookup("hEvap") >> hEvap_;
        lookup("R") >> R_;
        lookup("printPhaseChange") >> printPhaseChange_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
