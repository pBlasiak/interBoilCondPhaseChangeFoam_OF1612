surfaceScalarField gradT = fvc::snGrad(T);
    
surfaceScalarField wallHeatFlux = fvc::interpolate(alphaEff*cp*rho)*gradT;
const surfaceScalarField::Boundary& patchHeatFlux =
         wallHeatFlux.boundaryField();

Info<< "\nWall heat fluxes from alphaEff" << endl;
forAll(patchHeatFlux, patchi)
{
   if (typeid(mesh.boundary()[patchi]) == typeid(wallFvPatch))
    {
        Info<< mesh.boundary()[patchi].name()
            << ": Total "
            << gSum
               (
                   mesh.magSf().boundaryField()[patchi]
                  *patchHeatFlux[patchi]
               )
            << " [W] over "
            << gSum
               (
                   mesh.magSf().boundaryField()[patchi]
               )
            << " [m2] ("
            << gSum
               (
                   mesh.magSf().boundaryField()[patchi]
                  *patchHeatFlux[patchi]
               )/
               gSum 
               (
                   mesh.magSf().boundaryField()[patchi]
               )
            << " [W/m2])"
            << endl;
    }
}
