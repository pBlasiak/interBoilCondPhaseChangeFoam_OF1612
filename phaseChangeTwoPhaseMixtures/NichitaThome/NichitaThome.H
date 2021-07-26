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
    Foam::phaseChangeTwoPhaseMixtures::NichitaThome

Description
    NichitaThome phase change model.
	The model supports condensation and evaporation mass transfer.
	It is simplified sharp model.

    Reference:
    \verbatim
    B.A. Nichita, J.R. Thome, A level set method and a heat transfer model
    implemented into FLUENT for modeling of microscale two phase flows, in:
    AVT-178 Specialists’ Meeting on System Level Thermal Management for
    Enhanced Platform Efficiency, Bucharest, Romania, 2010.    
	\endverbatim


SourceFiles
    NichitaThome.C

\*--------------------------------------------------------------------*/

#ifndef NichitaThome_H
#define NichitaThome_H

#include "phaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeTwoPhaseMixtures
{

/*--------------------------------------------------------------------*\
                              Class NichitaThome
\*--------------------------------------------------------------------*/

class NichitaThome
:
    public phaseChangeTwoPhaseMixture
{
    // Private data

    	Switch cond_;

    	Switch evap_;

    // Private member functions

        volVectorField calcGradAlphal() const;
        volVectorField calcGradT() const;

public:

    //- Runtime type information
    TypeName("NichitaThome");


    // Constructors

        //- construct from components
        NichitaThome
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~NichitaThome()
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
