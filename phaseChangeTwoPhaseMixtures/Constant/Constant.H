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
    Foam::phaseChangeTwoPhaseMixtures::Constant

Description
    Constant phase change model.
	The model supports condensation and evaporation mass transfer.
	The mass flux due to phase change is set in the dictionary
	as constant value.



SourceFiles
    Constant.C

\*--------------------------------------------------------------------*/

#ifndef Constant_H
#define Constant_H

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{

/*--------------------------------------------------------------------*\
                              Class Constant
\*--------------------------------------------------------------------*/

class Constant
:
    public phaseChangeTwoPhaseMixture
{
    // Private data

	    // in kg/(m2 s)
    	dimensionedScalar mCondFlux_;
    	dimensionedScalar mEvapFlux_;

    // Private member functions
	
        // moze przeniesc ta funkcje do phaseChangeTwoPhaseMixture
		// i dodac zmienna gdzie jest przechowywany grad(alpha)?
		// zeby nie liczyc w kazdej funkcji tego od nowa
        volScalarField calcGradAlphal() const;

public:

    //- Runtime type information
    TypeName("Constant");


    // Constructors

        //- construct from components
        Constant
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~Constant()
        {}


    // Member Functions


        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply by (1 - alphal) for the condensation rate
        //  and a coefficient to multiply by alphal for the vaporisation rate
		//  in alphaEqn
        virtual Pair<tmp<volScalarField> > mDotAlphal() const;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply by (p - pSat) in pEqn
        virtual Pair<tmp<volScalarField> > mDotP() const;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply by (T - TSat) in TEqn
        virtual Pair<tmp<volScalarField> > mDotT() const;

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseChangeTwoPhaseMixtures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
