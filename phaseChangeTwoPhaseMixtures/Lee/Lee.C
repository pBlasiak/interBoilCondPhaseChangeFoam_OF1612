/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "Lee.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{
    defineTypeNameAndDebug(Lee, 0);
    addToRunTimeSelectionTable(phaseChangeTwoPhaseMixture, Lee, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeTwoPhaseMixtures::Lee::Lee
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    phaseChangeTwoPhaseMixture(typeName, U, phi),

    Cc_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cc")),
    Cv_(phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cv")),

    mcCoeff_(Cc_*rho2()),
    mvCoeff_(Cv_*rho1())
{
	//TODO: sprawdz czy model jest w solverze korektowany
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Lee::mDotAlphal() 
//Foam::phaseChangeTwoPhaseMixtures::Lee::mDotAlphal() const
{
    //const volScalarField& T = alpha1_.db().lookupObject<volScalarField>("T");
//	TSatLocal();
//    volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));

    const dimensionedScalar T0("0", dimTemperature, 0.0);

    return Pair<tmp<volScalarField> >
    (
		// minus sign "-" to provide mc > 0  and mv < 0
		-mcCoeff_*min(T_ - TSat_, T0)/TSat_,

		-mvCoeff_*max(T_ - TSat_, T0)/TSat_
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Lee::mDotP() 
//Foam::phaseChangeTwoPhaseMixtures::Lee::mDotP() const
{
    const volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));
    const dimensionedScalar T0("0", dimTemperature, 0.0);

    return Pair<tmp<volScalarField> >
    (
		// minus sign "-" to provide mc > 0  and mv < 0
        -mcCoeff_*(1.0 - limitedAlpha1)*min(T_ - TSat_, T0)/TSat_
		*pos(p_ - pSat_)/max(p_ - pSat_, 1E-8*pSat_),

        -mvCoeff_*limitedAlpha1*max(T_ - TSat_, T0)/TSat_
		*neg(p_ - pSat_)/min(p_ - pSat_, 1E-8*pSat_)
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::phaseChangeTwoPhaseMixtures::Lee::mDotT()
//Foam::phaseChangeTwoPhaseMixtures::Lee::mDotT() const
{
    volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));

    return Pair<tmp<volScalarField> >
    (
		// minus sign "-" to provide mc > 0  
		// now also mv > 0 but in TEqn, term (mc - mv) is applied
        -mcCoeff_*(1.0 - limitedAlpha1)*neg(T_ - TSat_)/TSat_,

        mvCoeff_*limitedAlpha1*pos(T_ - TSat_)/TSat_
    );
}

//void Foam::phaseChangeTwoPhaseMixtures::Lee::correct()
//{}


bool Foam::phaseChangeTwoPhaseMixtures::Lee::read()
{
    if (phaseChangeTwoPhaseMixture::read())
    {
        phaseChangeTwoPhaseMixtureCoeffs_ = subDict(type() + "Coeffs");

        phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cc") >> Cc_;
        phaseChangeTwoPhaseMixtureCoeffs_.lookup("Cv") >> Cv_;

        mcCoeff_ = Cc_*rho2();
        mvCoeff_ = Cv_*rho1();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
