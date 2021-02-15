# interBoilCondPhaseChangeFoam
This is a solver for boiling and condensation which is written based on OpenFOAM-v1612+ and solver interPhaseChangeFoam.
## Key features
* Smoothing of alphal field to decrease spurious (parasitic) currents
* Four mass transfer models 
* Galusinsky-Vignoaux criterion for timestep
* Mesh adaption solver

## Installation
It is working on OpenFOAM-v1612+
```bash
git clone https://github.com/NimaSam/phaseChangeHeatFoam/
cd phaseChangeHeatFoam/Application/
1) mkdir -p $FOAM_RUN

2) cd $WM_PROJECT_USER_DIR

3) mkdir -p applications/solvers/multiphase/interBoilCondPhaseChangeFoam
4) mkdir -p src

6) copy to src catalogue
https://github.com/pBlasiak/src_OF1612

7) cd src/transportModels
8) ./Allwmake
9) cd ../TurbulenceModels
10) ./Allwmake
11) cd ../../applications/solvers/multiphase/interBoilCondPhaseChangeFoam
12) copy to applications/solvers/multiphase/interBoilCondPhaseChangeFoam 
https://github.com/pBlasiak/interBoilCondPhaseChangeFoam_OF1612

13) ./Allwmake
```

## Some useful references
* [phaseChangeHeatFoam - Nima Sam's solver](https://github.com/NimaSam/phaseChangeHeatFoam)
* [Discussion on phaseChangeHeatFoam on CFD online](https://www.cfd-online.com/Forums/openfoam-solving/87665-evapphasechangefoam-5.html)
* [interThermalPhaseChangeFoam - article](https://www.sciencedirect.com/science/article/pii/S2352711016300309)
* [Jibran Heider thesis](https://www.researchgate.net/profile/Jibran_Haider/publication/259898900_Numerical_Modelling_of_Evaporation_and_Condensation_Phenomena/links/5738d95308ae9f741b2bda90/Numerical-Modelling-of-Evaporation-and-Condensation-Phenomena.pdf)
* [influence of surface tension implementation on VOF](https://www.sciencedirect.com/science/article/pii/S0301932213000190)
* [sclsVOFFoam](https://bitbucket.org/nunuma/public/src/d03747a27470214d29e43d5ffb4f1f4ef946c45d/OpenFOAM/solvers/2.0/?at=master)
* [S-CLSVOF](http://doras.dcu.ie/20019/1/PhDAAlbadawi.pdf)
* [S-CLSVOF](https://www.sciencedirect.com/science/article/pii/S0045793015003266)
* [S-CLSVOF](https://www.cfd-online.com/Forums/openfoam-solving/129732-clsvof-interfoam.html)
* [spurious currents discussion](https://www.cfd-online.com/Forums/openfoam-programming-development/189211-attempt-decrease-spurious-currents-vof.html)
* [spurious currents discussion](https://github.com/floquation/OF-kva_interfaceProperties)
* [Hakan Nillson's course - coupled LS and VOF](http://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2015/SankarMenon/Report_SankarMenon.pdf)
* [spurious currents suppression method](https://www.tandfonline.com/doi/pdf/10.1080/10407782.2014.916109?needAccess=true)
* [extra term in alphaEqn.H of interPhaseChangeFoam](https://www.cfd-online.com/Forums/openfoam-solving/138606-extra-term-alphaeqn-h-interphasechangefoam-version-2-3-0-a.html)
* [mesh adaption in OpenFOAM](https://www.cfd-online.com/Forums/openfoam-solving/131509-how-use-mesh-adaptation-openfoam.html)
* [origin of spurious currents](https://www.sciencedirect.com/science/article/pii/S0307904X05001666)
* [benchmarks for VOF and solver in OF with smoothing](https://www.sciencedirect.com/science/article/pii/S0045793013002612?via%3Dihub)
* [benchmarks for VOF and solver in OF with smoothing](https://www.cfd-online.com/Forums/openfoam-verification-validation/124363-interfoam-validation-bubble-droplet-flows-microfluidics.html)
* [sclsInterFoam](https://www.researchgate.net/publication/318233965_sCLSinterFoam_OpenFOAM220)
* [interFoam - Berberovic's article](https://journals.aps.org/pre/pdf/10.1103/PhysRevE.79.036306)


## TODO
* Level-Set method coupled with VOF
* Adding new mass transfer models (e.g. sharp model)
* Enabling conjugate heat transfer
* Method to further reduce spurious currents
* Geometric reconstruct of the interface

## Flowchart
* [Flowchart](https://github.com/pBlasiak/interBoilCondPhaseChangeFoam_OF1612/blob/master/interBoilCondPhaseChangeFoam_flowchart.html)

