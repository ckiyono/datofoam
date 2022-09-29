{
	volTensorField tau         = turbulence->muEff()*(fvc::grad(U) + fvc::grad(U)().T() - (2/3)*fvc::div(U)*I);
    volScalarField viscousDiss = fvc::div(tau & U);
	
	volVectorField Urot = (omegaVector ^ mesh.C());
	
    volScalarField& he = thermo.he();
	volScalarField i_comp = 0.5*(magSqr(Urot) - magSqr(U)); // rothalpy complement

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      - fvc::div(phi, i_comp)
      - fvm::laplacian(turbulence->alphaEff(), he)
	  - viscousDiss
     ==
        fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
}
