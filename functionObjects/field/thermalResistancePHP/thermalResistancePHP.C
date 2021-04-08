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
    writeHeader(os, "Wall heat-flux");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "integral");
    os  << endl;
}


void Foam::functionObjects::thermalResistancePHP::calcHeatFlux
(
    const volScalarField& alpha,
    const volScalarField& he,
    volScalarField& thermalResistancePHP
)
{
    surfaceScalarField heatFlux(fvc::interpolate(alpha)*fvc::snGrad(he));

    volScalarField::Boundary& thermalResistancePHPBf = thermalResistancePHP.boundaryFieldRef();

    const surfaceScalarField::Boundary& heatFluxBf = heatFlux.boundaryField();

	thermalResistancePHPBf[evapPatch_] = heatFluxBf[evapPatch_];
    //forAll(thermalResistancePHPBf, patchi)
    //{
    //    thermalResistancePHPBf[patchi] = heatFluxBf[patchi];
    //}

    //if (foundObject<volScalarField>("Qr"))
    //{
    //    const volScalarField& Qr = lookupObject<volScalarField>("Qr");

    //    const volScalarField::Boundary& radHeatFluxBf = Qr.boundaryField();

    //    forAll(thermalResistancePHPBf, patchi)
    //    {
    //        thermalResistancePHPBf[patchi] += radHeatFluxBf[patchi];
    //    }
    //}
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
    patchSet_()
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
            dimensionedScalar("0", dimMass/pow3(dimTime), 0)
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

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        forAllConstIter(labelHashSet, patchSet_, iter)
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

        patchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::thermalResistancePHP::execute()
{
    volScalarField& thermalResistancePHP = const_cast<volScalarField&>
    (
        lookupObject<volScalarField>(type())
    );

    const scalarField& alphaEffp =
        patch().lookupPatchField<volScalarField, scalar>(alphaEffName_);
		const scalarField& cp0 =
			patch().lookupPatchField<volScalarField, scalar>("cp");
		const scalarField& rho =
			patch().lookupPatchField<volScalarField, scalar>("rho");
    //if
    //(
    //    foundObject<compressible::turbulenceModel>
    //    (
    //        turbulenceModel::propertiesName
    //    )
    //)
    //{
    //    const compressible::turbulenceModel& turbModel =
    //        lookupObject<compressible::turbulenceModel>
    //        (
    //            turbulenceModel::propertiesName
    //        );

    //    calcHeatFlux
    //    (
    //        turbModel.alphaEff()(),
    //        turbModel.transport().he(),
    //        thermalResistancePHP
    //    );
    //}
    //else if (foundObject<fluidThermo>(fluidThermo::dictName))
    //{
    //    const fluidThermo& thermo =
    //        lookupObject<fluidThermo>(fluidThermo::dictName);

    //    calcHeatFlux
    //    (
    //        thermo.alpha(),
    //        thermo.he(),
    //        thermalResistancePHP
    //    );
    //}
    //else
    //{
    //    FatalErrorInFunction
    //        << "Unable to find compressible turbulence model in the "
    //        << "database" << exit(FatalError);
    //}

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

    forAllConstIter(labelHashSet, patchSet_, iter)
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