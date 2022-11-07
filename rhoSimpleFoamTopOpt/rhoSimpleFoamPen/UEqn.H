    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(U)
	  + fvm::Sp(alphaPen, U)                      // Penalization
	  + rho*(omegaVector^(omegaVector^mesh.C()))  // Centrifugal Acceleration
      + 2.0*rho*(omegaVector^U)                   // Coriolis Acceleration
     ==
        fvOptions(rho, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
