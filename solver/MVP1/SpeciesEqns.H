{
    volScalarField& diffEffCation = diffCoCation;
    volScalarField& diffEffAnion = diffCoAnion;
    phaseModel& liquidPhase = fluid.phase(0);
    volScalarField alpha_liquid = liquidPhase;

    forAll(mesh.C(), i)
    {
        diffEffCation[i] = diffCationInLiquid.value()*pow(alpha_liquid[i], 2);
        diffEffAnion[i] = diffAnionInLiquid.value()*pow(alpha_liquid[i], 2);

        if (diffEffCation[i] < 1e-16)
        {
            diffEffCation[i] = 0.0;
        }

        if (diffEffAnion[i] < 1e-16)
        {
            diffEffAnion[i] = 0.0;
        }
    }

    Info<< "min/max(diffEffCation) = "
        << min(diffEffCation).value() << ", "
        << max(diffEffCation).value() << endl;
    Info<< "min/max(diffEffAnion) = "
        << min(diffEffAnion).value() << ", "
        << max(diffEffAnion).value() << endl;

    const tmp<volScalarField> tdmdtReaction = fluid.updateDmdt(Ccation, SuOut);
    const volScalarField& dmdtReaction = tdmdtReaction();

    volScalarField precipMolarRate
    (
        IOobject
        (
            "precipMolarRate",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dmdtReaction/molW
    );

    SuCationOut = nuCation*precipMolarRate;
    SuAnionOut = nuAnion*precipMolarRate;

    fvScalarMatrix CcationEqn
    (
        fvm::ddt(Ccation)
      + fvm::div(phi, Ccation, "div(phi,Ccation)")
      - fvm::laplacian(diffEffCation, Ccation, "laplacian(diffEffCation,Ccation)")
      ==
       -nuCation*precipMolarRate
      + fvOptions(Ccation)
    );

    CcationEqn.relax();
    fvOptions.constrain(CcationEqn);
    CcationEqn.solve();
    fvOptions.correct(Ccation);

    fvScalarMatrix CanionEqn
    (
        fvm::ddt(Canion)
      + fvm::div(phi, Canion, "div(phi,Canion)")
      - fvm::laplacian(diffEffAnion, Canion, "laplacian(diffEffAnion,Canion)")
      ==
       -nuAnion*precipMolarRate
      + fvOptions(Canion)
    );

    CanionEqn.relax();
    fvOptions.constrain(CanionEqn);
    CanionEqn.solve();
    fvOptions.correct(Canion);

    Info<< "min/max(Ccation) = "
        << min(Ccation).value() << ", " << max(Ccation).value() << endl;
    Info<< "min/max(Canion) = "
        << min(Canion).value() << ", " << max(Canion).value() << endl;
}
