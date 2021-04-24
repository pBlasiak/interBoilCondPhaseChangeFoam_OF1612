gradT = fvc::snGrad(T);

wallHeatFlux = fvc::interpolate(alphaEff*cp*rho)*gradT;

const label evapPatchID = mesh.boundaryMesh().findPatchID(evapPatch);
const label condPatchID = mesh.boundaryMesh().findPatchID(condPatch);
const scalarField Tbe = T.boundaryField()[evapPatchID];
const scalarField Tbc = T.boundaryField()[condPatchID];
const dimensionedScalar Te = gSum(mesh.magSf().boundaryField()[evapPatchID]*Tbe)
	/gSum(mesh.magSf().boundaryField()[evapPatchID]);
const dimensionedScalar Tc = gSum(mesh.magSf().boundaryField()[condPatchID]*Tbc)
	/gSum(mesh.magSf().boundaryField()[condPatchID]);
thermalResistancePHP.value() = 
	(Te.value() - Tc.value())/
	gSum
    (
        mesh.magSf().boundaryField()[evapPatchID]
	   *wallHeatFlux.boundaryField()[evapPatchID]
	);

Info<< "\nThermal resistance of the PHP = " << thermalResistancePHP.value() << " [K/W]" << endl;

if (printWallHeatFluxes)
{
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
}
