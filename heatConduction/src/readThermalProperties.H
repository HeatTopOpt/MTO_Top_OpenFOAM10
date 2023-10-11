Info<< "Reading thermal properties\n" << endl;
IOdictionary thermalProperties
(
    IOobject
    (
        "thermalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar k(thermalProperties.lookup("k"));//Thermal Conductivity 
dimensionedScalar kmin(1e-3*k);
dimensionedScalar b2("b2", dimensionSet(0,0,-1,2-Pnorm,0,0,0),1);

volScalarField DT
(
    IOobject
    (
        "DT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    kmin+(k-kmin)*xh*xh*xh
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
    volScalarField Ta
    (
        IOobject
        (
            "Ta",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        T
    );
Info<< "Reading heatsource Q\n" << endl;
volScalarField Q
(
    IOobject
    (
        "Q",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Q", dimensionSet(1,-1,-3,0,0,0,0),1),
    calculatedFvPatchScalarField::typeName
);
volScalarField fsens
(
    IOobject
    (
        "fsens",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
   -3*(k-kmin)*xp*xp*(fvc::grad(T)&fvc::grad(T)),
   zeroGradientFvPatchScalarField::typeName
);
volScalarField fsensh(fsens);
volScalarField gsens
(
    IOobject
    (
        "gsens",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    xh,
    zeroGradientFvPatchScalarField::typeName
);
volScalarField gsensh(gsens);