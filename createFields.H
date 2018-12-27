Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mesh,
	dimensionedScalar("rho0", dimDensity, 0.0)
);
rho.oldTime();

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Creating phaseChangeTwoPhaseMixture\n" << endl;
autoPtr<phaseChangeTwoPhaseMixture> mixture =
    phaseChangeTwoPhaseMixture::New(U, phi);

volScalarField& alpha1(mixture->alpha1());
volScalarField& alpha2(mixture->alpha2());

const dimensionedScalar& rho1 = mixture->rho1();
const dimensionedScalar& rho2 = mixture->rho2();
const dimensionedScalar& k1 = mixture->k1();
const dimensionedScalar& k2 = mixture->k2();
const dimensionedScalar& cp1 = mixture->cp1();
const dimensionedScalar& cp2 = mixture->cp2();
const dimensionedScalar& hEvap = mixture->hEvap();
const dimensionedScalar& pSat = mixture->pSat();
const dimensionedScalar& cutoff = mixture->cutoff();

rho = alpha1*rho1 + alpha2*rho2;
p = p_rgh + rho*gh;
mixture->correct();
volScalarField TSat = mixture->TSat();


//    surfaceScalarField rhoPhi
//    (
//        IOobject
//        (
//            "rhoPhi",
//            runTime.timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        rho1*phi
//    );

    volScalarField alphaEff
    (
        IOobject
        (
            "alphaEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1*k1 + alpha2*k2
    );

    volScalarField cp
    (
        IOobject
        (
            "cp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1*cp1 + alpha2*cp2
    );

// Need to store rho for ddt(rhoCp, U)
    volScalarField rhoCp
    (
        IOobject
        (
            "rhoCp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*alpha1*cp1 + rho2*alpha2*cp2
    );
	rhoCp.oldTime();

    //volScalarField mCond
    //(
    //    IOobject
    //    (
    //        "mCond",
    //        runTime.timeName(),
    //        mesh,
    //        IOobject::NO_READ,
    //        IOobject::NO_WRITE
    //        //IOobject::AUTO_WRITE
    //    ),
    //    //mixture->mDotAlphal()[0]()*(1.0 - alpha1)
    //    mixture->mCondDotAlphal()*(1.0 - alpha1)
    //);

    //volScalarField mVap
    //(
    //    IOobject
    //    (
    //        "mVap",
    //        runTime.timeName(),
    //        mesh,
    //        IOobject::NO_READ,
    //        IOobject::NO_WRITE
    //        //IOobject::AUTO_WRITE
    //    ),
    //    //mixture->mDotAlphal()[1]*alpha1
    //    mixture->mEvapDotAlphal()*alpha1
    //);

// Construct interface from alpha1 distribution
smoothInterfaceProperties interface(alpha1, U, mixture());

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture())
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());
