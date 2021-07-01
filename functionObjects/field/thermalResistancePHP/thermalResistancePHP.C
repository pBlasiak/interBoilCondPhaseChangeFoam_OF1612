/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "thermalResistancePHP.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(thermalResistancePHP, 0);
    addToRunTimeSelectionTable(functionObject, thermalResistancePHP, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::thermalResistancePHP::writeFileHeader(Ostream& os) const
{
    // Add headers to output data
    writeHeader(os, "Thermal resistance");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "integral");
    os  << endl;
}


void Foam::functionObjects::thermalResistancePHP::calcThermalResistancePHP
(
    const volScalarField& alpha,
    const volScalarField& T,
    const volScalarField& cp,
    const volScalarField& rho,
    volScalarField& thermalResistancePHP
)
{
    surfaceScalarField heatFlux(fvc::interpolate(alpha*cp*rho)*fvc::snGrad(T));

    volScalarField::Boundary& thermalResistancePHPBf = thermalResistancePHP.boundaryFieldRef();

    const surfaceScalarField::Boundary& heatFluxBf = heatFlux.boundaryField();

    const volScalarField::Boundary& TBf = T.boundaryField();

	dimensionedScalar Tbe;
	dimensionedScalar Ae;
	forAll(evapPatchSet_, patchi)
	{
		forAll(TBf, patchID)
		{
			Tbe += mesh_.magSf().boundaryField()[patchi][patchID]*TBf[patchi][patchID];
			Ae  += mesh_.magSf().boundaryField()[patchi][patchID];
		}
	}
	dimensionedScalar TevapAve = Tbe/Ae;
	Info<< "TevapAve = " << TevapAve << endl;
		
	dimensionedScalar Tbc;
	dimensionedScalar Ac;
	forAll(condPatchSet_, patchi)
	{
		forAll(TBf, patchID)
		{
			Tbc += mesh_.magSf().boundaryField()[patchi][patchID]*TBf[patchi][patchID];
			Ac  += mesh_.magSf().boundaryField()[patchi][patchID];
		}
	}
	dimensionedScalar TcondAve = Tbc/Ac;
	Info<< "TcondAve = " << TcondAve << endl;

	dimensionedScalar Q;
	forAll(evapPatchSet_, patchi)
	{
		forAll(heatFluxBf, patchID)
		{
			Q += mesh_.magSf().boundaryField()[patchi][patchID]*heatFluxBf[patchi][patchID];
		}
	}

	thermalResistancePHP = (TevapAve - TcondAve)/Q;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::thermalResistancePHP::thermalResistancePHP
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    evapPatchSet_(),
    condPatchSet_()
{
    volScalarField* thermalResistancePHPPtr
    (
        new volScalarField
        (
            IOobject
            (
                type(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimTemperature/dimPower, 0)
        )
    );

    mesh_.objectRegistry::store(thermalResistancePHPPtr);

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::thermalResistancePHP::~thermalResistancePHP()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::thermalResistancePHP::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    evapPatchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("evapPatch", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (evapPatchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                evapPatchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing evaporator patches: " << nl;
        labelHashSet filteredPatchSet;
        forAllConstIter(labelHashSet, evapPatchSet_, iter)
        {
            label patchi = iter.key();
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall heat-flux on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        evapPatchSet_ = filteredPatchSet;
    }

    condPatchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("condPatch", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (condPatchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                condPatchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing condenser patches: " << nl;
        labelHashSet filteredPatchSet;
        forAllConstIter(labelHashSet, condPatchSet_, iter)
        {
            label patchi = iter.key();
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall heat-flux on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        condPatchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::thermalResistancePHP::execute()
{
    volScalarField& thermalResistancePHP = const_cast<volScalarField&>
    (
        lookupObject<volScalarField>(type())
    );

    const volScalarField& alphaEffp =
        lookupObject<volScalarField>("alphaEff");
    const volScalarField& T =
        lookupObject<volScalarField>("T");
    const volScalarField& cp =
        lookupObject<volScalarField>("cp");
    const volScalarField& rho =
        lookupObject<volScalarField>("rho");

    calcThermalResistancePHP
    (
	    alphaEffp,
		T,
		cp, 
		rho,
        thermalResistancePHP
    );

    return true;
}


bool Foam::functionObjects::thermalResistancePHP::write()
{
    const volScalarField& thermalResistancePHP = lookupObject<volScalarField>(type());

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << thermalResistancePHP.name() << endl;

    thermalResistancePHP.write();

    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf =
        mesh_.magSf().boundaryField();

    forAllConstIter(labelHashSet, evapPatchSet_, iter)
    {
        label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = thermalResistancePHP.boundaryField()[patchi];

        const scalar minHfp = gMin(hfp);
        const scalar maxHfp = gMax(hfp);
        const scalar integralHfp = gSum(magSf[patchi]*hfp);

        if (Pstream::master())
        {
            file()
                << mesh_.time().value()
                << token::TAB << pp.name()
                << token::TAB << minHfp
                << token::TAB << maxHfp
                << token::TAB << integralHfp
                << endl;
        }

        Log << "    min/max(" << pp.name() << ") = "
            << minHfp << ", " << maxHfp << ", " << integralHfp << endl;
    }

    return true;
}


// ************************************************************************* //
