/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

//#include "fvCFD.H"
#include "thermalIncompressibleTwoPhaseMixture.H"
#include "fvc.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalIncompressibleTwoPhaseMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalIncompressibleTwoPhaseMixture::thermalIncompressibleTwoPhaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    incompressibleTwoPhaseMixture(U, phi),

    k1_
    (
        "k1_",
        dimEnergy/dimTime/dimLength/dimTemperature,
        subDict(phase1Name_).lookup("k")
    ),
    k2_
    (
        "k2_",
        k1_.dimensions(),
        subDict(phase2Name_).lookup("k")), cp1_ (
        "cp1",
        dimEnergy/dimTemperature/dimMass,
        subDict(phase1Name_).lookup("cp")
    ),
    cp2_
    (
        "cp2",
        dimEnergy/dimTemperature/dimMass,
        subDict(phase2Name_).lookup("cp")
    )//,

   // hEvap_
   // (
   //     "hEvap",
   //     dimEnergy/dimMass,
   //     subDict(phase2Name_).lookup("hEvap")
   // )
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField> Foam::thermalIncompressibleTwoPhaseMixture::k() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
		new volScalarField
        (
            "k",
            limitedAlpha1*k1_
          + (scalar(1) - limitedAlpha1)*k2_
        )
    );
}

Foam::tmp<Foam::surfaceScalarField> 
Foam::thermalIncompressibleTwoPhaseMixture::kfHarmonic() const
{
    const surfaceScalarField limitedAlpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>
    (
		new surfaceScalarField
        (
            "kfHarmonic",
			1.0/
			(
				(1.0/k1_ - 1.0/k2_)*limitedAlpha1f
			  + 1.0/k2_
			)
        )
    );
}

Foam::tmp<Foam::surfaceScalarField> 
Foam::thermalIncompressibleTwoPhaseMixture::kfDensityHarmonic() const
{
    const surfaceScalarField limitedAlpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>
    (
		new surfaceScalarField
        (
            "kfHarmonic",
			(
				k1_*limitedAlpha1f/rho1_ 
			  - (1.0 - limitedAlpha1f)*k2_/rho2_ 
			)/(limitedAlpha1f/rho1_ - (1.0 - limitedAlpha1f)/rho2_)
        )
    );
}

Foam::tmp<Foam::volScalarField> Foam::thermalIncompressibleTwoPhaseMixture::cp() const
{
	const volScalarField limitedAlpha1
	(
		min(max(alpha1_, scalar(0)), scalar(1))
	);

	return tmp<volScalarField>
	(
		new volScalarField
		(
			"cp",
			( 
				cp1_*rho1_*limitedAlpha1 + cp2_*rho2_*(scalar(1) - limitedAlpha1) 
			)/( rho1_*limitedAlpha1 + rho2_*(scalar(1) - limitedAlpha1) )
		)
	);
}

Foam::tmp<Foam::volScalarField> Foam::thermalIncompressibleTwoPhaseMixture::alphaEff() const
{
	const volScalarField limitedAlpha1
	(
		min(max(alpha1_, scalar(0)), scalar(1))
	);

	return tmp<volScalarField>
	(
		new volScalarField
		(
			"alphaEff",
			k()/( cp()*( rho1_*limitedAlpha1 + rho2_*(scalar(1) - limitedAlpha1) ) )
		)
	);
}

bool Foam::thermalIncompressibleTwoPhaseMixture::read()
{
    if (incompressibleTwoPhaseMixture::read())
    {
        subDict(phase1Name_).lookup("k1") >> k1_;
        subDict(phase2Name_).lookup("k2") >> k2_;

        subDict(phase1Name_).lookup("cp") >> cp1_;
        subDict(phase2Name_).lookup("cp") >> cp2_;

      //  subDict(phase2Name_).lookup("hEvap") >> hEvap_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
