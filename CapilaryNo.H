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

Global
    FourierNo

Description
    Calculates and outputs the mean and maximum Capilary Numbers.
	According to floquation:
	https://www.cfd-online.com/Forums/openfoam-programming-development/189211-attempt-decrease-spurious-currents-vof.html

\*---------------------------------------------------------------------------*/

scalar maxCapilary
(
    readScalar(runTime.controlDict().lookup("maxCapilary"))
);

scalar CapilaryNum = 0.0;
scalar meanCapilaryNum = 0.0;

    CapilaryNum = gMax(  sqrt(mixture->rho()*pow(CellDims.field(), 3)/
						(4.0*Foam::constant::mathematical::pi*interface.sigma().value()))  );

    meanCapilaryNum = gAverage(  sqrt(mixture->rho()*pow(CellDims.field(), 3)/
						(4.0*Foam::constant::mathematical::pi*interface.sigma().value()))  );

Info<< "Capilary Number mean: " << meanCapilaryNum
    << " max: " << CapilaryNum << endl;

// ************************************************************************* //
