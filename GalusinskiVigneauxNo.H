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
    GalusinskiVigneauxNum

Description
    Calculates and outputs the mean and maximum GalusinskiVigneaux Numbers.

References:
	[1] Galusinski C. and Vigneaux P. 2008 
	    On stability condition for bifluid flows with surface tension: application to
	    microfluidics, J. Comput. Phys., 227, 6140–64

	[2] Suraj S. Deshpande, Lakshman Anumolu and Mario F. Trujillo 
		Evaluating the performance of the two-phase flow solver interFoam

\*---------------------------------------------------------------------------*/

const scalar C1 = 0.01; // from [2]
const scalar C2 = 10.0; // from [2]

//scalar maxGalusinskiVigneaux
//(
//    readScalar(runTime.controlDict().lookup("maxGalusinskiVigneaux"))
//);

scalar GalusinskiVigneauxNum = 0.0;
//scalar meanGalusinskiVigneauxNum = 0.0;
scalarField tauRho(0.0);
scalarField tauMu(0.0);
//scalar GVdt = 0.0;

	tauRho = sqrt(mixture->rho()*pow(CellDims.field(), 3)/interface.sigma().value());
	tauMu  = mixture->mu()*CellDims.field()/interface.sigma().value();

    GalusinskiVigneauxNum = gMax( 1.0/(0.5*(C2*tauMu + sqrt(pow(C2*tauMu, 2) + 4*C1*tauRho*tauRho)))*
								runTime.deltaTValue());
    //GVdt = gMin( 0.5*(C2*tauMu + sqrt(pow(C2*tauMu, 2) + 4*C1*tauRho*tauRho))  );

    //meanGalusinskiVigneauxNum = gAverage( 1.0/(0.5*(C2*tauMu + sqrt(pow(C2*tauMu, 2) + 4*C1*tauRho*tauRho)))*
//								runTime.deltaTValue());

//Info<< "GalusinskiVigneaux Number mean: " << meanGalusinskiVigneauxNum
   Info << "GalusinskiVigneauxNum max: " << GalusinskiVigneauxNum << endl;
//Info<< "GV dt: " << GVdt << endl;

// ************************************************************************* //
