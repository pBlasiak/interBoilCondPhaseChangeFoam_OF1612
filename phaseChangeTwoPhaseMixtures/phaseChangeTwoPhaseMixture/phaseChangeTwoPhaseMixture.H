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

Class
    Foam::phaseChangeTwoPhaseMixture

Description

SourceFiles
    phaseChangeTwoPhaseMixture.C
    newPhaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeTwoPhaseMixture_H
#define phaseChangeTwoPhaseMixture_H

#include "thermalIncompressibleTwoPhaseMixture.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class phaseChangeTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeTwoPhaseMixture
:
    public thermalIncompressibleTwoPhaseMixture
{

protected:

    // Protected data

        dictionary phaseChangeTwoPhaseMixtureCoeffs_;

		//- Condensation mass flux, kg/s (+)
		volScalarField mCond_;
		
		//- Vaporisation mass flux, kg/s (-)
		volScalarField mEvap_;
		
		//- Used in mDotAlphal()
		volScalarField mCondAlphal_;
		volScalarField mEvapAlphal_;

		//- Used in mDotP()
		volScalarField mCondP_;
		volScalarField mEvapP_;

		//- Used in mDotT()
		volScalarField mCondT_;
		volScalarField mEvapT_;

        //- Reference to pressure field
        const volScalarField& p_;

        //- Reference to temperature field
        const volScalarField& T_;

        //- Global saturation temperature
        dimensionedScalar TSatG_;
		
        //- Local saturation temperature
        volScalarField TSat_;

        //- Saturation vapour pressure
        dimensionedScalar pSat_;

        //- Latent heat of vaporization
        dimensionedScalar hEvap_;

        //- Specific gas constant
        dimensionedScalar R_;

		//-True if saturation vapour temperature is based on pressure
		Switch TSatLocalPressure_;

		Switch printPhaseChange_;
		

    // Private Member Functions

        //- Disallow copy construct
        phaseChangeTwoPhaseMixture(const phaseChangeTwoPhaseMixture&);

        //- Disallow default bitwise assignment
        void operator=(const phaseChangeTwoPhaseMixture&);

//protected:

		//- Calculates local saturation temperature field based on pressure
		//  using simplified Clausius-Clapeyron equation 
		void calcTSatLocal();

public:

    //- Runtime type information
    TypeName("phaseChangeTwoPhaseMixture");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseChangeTwoPhaseMixture,
            components,
            (
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (U, phi)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<phaseChangeTwoPhaseMixture> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Constructors

        //- Construct from components
        phaseChangeTwoPhaseMixture
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~phaseChangeTwoPhaseMixture()
    {}


    // Member Functions

        //- Return const-access to the condensation mass flux
        const volScalarField& mCond() const
        {
            return mCond_;
        }
		
        //- Return const-access to the evaporation mass flux
        const volScalarField& mEvap() const
        {
            return mEvap_;
        }

        //- Return const-access to the global saturation temperature
        const dimensionedScalar& TSatG() const
        {
            return TSatG_;
        }
		
        //- Return const-access to the local saturation temperature
        const volScalarField& TSat() const
        {
            return TSat_;
        }

        //- Return const-access to the global saturation vapour pressure
        const dimensionedScalar& pSat() const
        {
            return pSat_;
        }

        //- Return const-access to the latent heat of vaporization
        const dimensionedScalar& hEvap() const
        {
            return hEvap_;
        }

        //- Return const-access to the specific gas constant
        const dimensionedScalar& R() const
        {
            return R_;
        }

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply by (1 - alphal) for the condensation rate
        //  and a coefficient to multiply by alphal for the vaporisation rate
		//  in alphaEqn 
		//  units are kg/(m3 s)
        virtual Pair<tmp<volScalarField> > mDotAlphal() const = 0;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply by (p - pSat) in pEqn
		//  units are kg/(m3 s Pa)
        virtual Pair<tmp<volScalarField> > mDotP() const = 0;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply by (T - TSat) in TEqn
		//  units are kg/(m3 s K)
        virtual Pair<tmp<volScalarField> > mDotT() const = 0;

        //- Return the volumetric condensation and vaporisation rates as a
        //  coefficient to multiply by (1 - alphal) for the condensation rate
        //  and a coefficient to multiply by alphal for the vaporisation rate
		//  in alphaEqn
        Pair<tmp<volScalarField> > vDotAlphal() const;

        //- Return the volumetric condensation and vaporisation rates as
        //  coefficients to multiply by (p - pSat)
        Pair<tmp<volScalarField> > vDotP() const;
		
        //- Return the volumetric condensation and vaporisation rates as
        //  coefficients to multiply by (T - TSat)
        Pair<tmp<volScalarField> > vDotT() const;

        //- Correct the phaseChange model
        virtual void correct();

        //- Read the transportProperties dictionary and update
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
