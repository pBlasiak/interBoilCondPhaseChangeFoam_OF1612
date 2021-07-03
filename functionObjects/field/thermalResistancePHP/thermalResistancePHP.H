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
    volScalarField 'thermalResistancePHP'.

    All wall patches are included by default; to restrict the calculation to
    certain evaporator and condenser patches, use the optional 
	'evapPatch' and 'condPatch' entries.

Usage
    Example of function object specification:
    \verbatim
    thermalResistancePHP1
	{
	    // Mandatory entries (unmodifiable)
	    type            thermalResistancePHP;
	    libs            ("libfieldFunctionObjectsPB.so");
	
	    // Optional entries (runtime modifiable)
	    //patches     (<patch1> ... <patchN>); // (wall1 "(wall2|wall3)");
	    evapPatch     (H_WALL); // (wall1 "(wall2|wall3)");
	    condPatch     (C_WALL); // (wall1 "(wall2|wall3)");
	
	    // Optional (inherited) entries
	    writePrecision  8;
	    writeToFile     true;
	    useUserTime     true;
	    region          region0;
	    enabled         true;
	    log             true;
	    timeStart       0;
	    timeEnd         1000;
	    executeControl  timeStep;
	    executeInterval 1;
	    writeControl    timeStep;
	    writeInterval   50;
	}
    \endverbatim

    Where the entries comprise:
    \table
        Property   | Description                        | Required   | Default value
        type       | type name: thermalResistancePHP    | yes        |
        evapPatch  | list of patches to process         | no         | all wall patches
        condPatch  | list of patches to process         | no         | all wall patches
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
        labelHashSet evapPatchSet_;
        labelHashSet condPatchSet_;

		

    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate the thermal resistance of a PHP
        void calcThermalResistancePHP
        (
            const volScalarField& alpha,
            const volScalarField& T,
            const volScalarField& cp,
            const volScalarField& rho,
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
