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

Class
    Foam::functionObjects::thermalResistancePHP

Group
    grpForcesFunctionObjects

Description
    Calculates and write the thermal resistance of a PHP as the
    dimensionedScalar 'thermalResistancePHP'.

    All wall patches are included by default; to restrict the calculation to
    certain patches, use the optional 'patches' entry.

Usage
    Example of function object specification:
    \verbatim
    thermalResistancePHP1
    {
        type        thermalResistancePHP;
        libs        ("libfieldFunctionObjects.so");
        ...
        patches     (".*Wall");
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property | Description                | Required   | Default value
        type     | type name: thermalResistancePHP    | yes        |
        patches  | list of patches to process | no         | all wall patches
    \endtable

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile
    Foam::functionObjects::timeControl

SourceFiles
    thermalResistancePHP.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_thermalResistancePHP_H
#define functionObjects_thermalResistancePHP_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class thermalResistancePHP Declaration
\*---------------------------------------------------------------------------*/

class thermalResistancePHP
:
    public fvMeshFunctionObject,
    public writeFile
{

protected:

    // Protected data

        //- Optional list of patches to process
        labelHashSet patchSet_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate the heat-flux
        void calcHeatFlux
        (
            const volScalarField& alpha,
            const volScalarField& he,
            volScalarField& thermalResistancePHP
        );


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        thermalResistancePHP(const thermalResistancePHP&);

        //- Disallow default bitwise assignment
        void operator=(const thermalResistancePHP&);


public:

    //- Runtime type information
    TypeName("thermalResistancePHP");


    // Constructors

        //- Construct from Time and dictionary
        thermalResistancePHP
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~thermalResistancePHP();


    // Member Functions

        //- Read the thermalResistancePHP data
        virtual bool read(const dictionary&);

        //- Calculate the thermal resistance 
        virtual bool execute();

        //- Write the thermal resistance
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
