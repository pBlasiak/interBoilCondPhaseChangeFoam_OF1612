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
    under the terms of the GNU Generac License as published by the
    Free Software Foundation; either  2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the ho it will be useful, but WITHOUT
    ANY WARRANTY; without even the imarranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.he GNU General Public License
    for more details.

    You should have received a copy oNU General Public License
    along with OpenFOAM; if not, write Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floorn, MA 02110-1301 USA

Class
    Foam::phaseChangeTwoPhaseMixtures::Tanasawa

Description
    Tanasawa phase change model.
	The model supports condensation and evaporation mass transfer.

    Reference:
    \verbatim
        Tanasawa, I.
        "Advances in condensation heat transfer",
        Advances in Heat Transfer, vol. 21, 1991, pp. 55-139.
    \endverbatim


SourceFiles
    Tanasawa.C

\*--------------------------------------------------------------------*/

#ifndef Tanasawa_H
#define Tanasawa_H

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{

/*--------------------------------------------------------------------*\
                              Class Tanasawa
\*--------------------------------------------------------------------*/

class Tanasawa
:
    public phaseChangeTwoPhaseMixture
{
    // Private data

    	Switch cond_;

    	Switch evap_;

    	dimensionedScalar gamma_;

        dimensionedScalar mCoeff_;
		
    // Private member functions

        volScalarField calcGradAlphal() const;

public:

    //- Runtime type information
    TypeName("Tanasawa");


    // Constructors

        //- construct from components
        Tanasawa
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~Tanasawa()
        {}


    // Member Functions


        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
		//  in alphaEqn
        //virtual Pair<tmp<volScalarField> > mDotAlphal() const;
        virtual Pair<tmp<volScalarField> > mDotAlphal();

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (p - pSat) in pEqn
        //virtual Pair<tmp<volScalarField> > mDotP() const;
        virtual Pair<tmp<volScalarField> > mDotP();

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (T - TSat) in TEqn
        //virtual Pair<tmp<volScalarField> > mDotT() const;
        virtual Pair<tmp<volScalarField> > mDotT();

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseChangeTwoPhaseMixtures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
