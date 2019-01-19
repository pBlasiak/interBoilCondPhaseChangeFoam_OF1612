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
    Calculates and outputs the mean and maximum Fourier Numbers.

\*---------------------------------------------------------------------------*/

scalar maxFourier
(
    readScalar(runTime.controlDict().lookup("maxFourier"))
);

scalar FourierNum = 0.0;
scalar meanFourierNum = 0.0;
scalar fourierDt = 0.0;

    FourierNum = gMax(  (mixture->alphaEff() + turbulence->nut())
				            *runTime.deltaTValue()/pow(CellDims.field(), 2));

    fourierDt = gMin(  pow(CellDims.field(), 2)/(mixture->alphaEff() + turbulence->nut()));

    meanFourierNum = gAverage(  (mixture->alphaEff() + turbulence->nut())
				            *runTime.deltaTValue()/pow(CellDims.field(), 2));

Info<< "Fourier Number mean: " << meanFourierNum
    << " max: " << FourierNum << endl;
Info<< "Fourier dt: " << fourierDt << endl;

// ************************************************************************* //
