// --- Momentum Equation (UEqn.H modified) --- //
fvVectorMatrix UEqn
(
	fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
  + fvc::grad(p)
  + fvm::Sp(alphaPen, U)                      // Penalization
  + rho*(omegaVector^(omegaVector^mesh.C()))  // Centrifugal Acceleration
  + 2.0*rho*(omegaVector^U)                   // Coriolis Acceleration
);

UEqn.relax();

// --- Energy Equation (EEqn.H modified) --- //
volTensorField tau         = turbulence->muEff()*(fvc::grad(U) + fvc::grad(U)().T() - (2/3)*fvc::div(U)*I);
volScalarField viscousDiss = fvc::div(tau & U);

volVectorField Urot = (omegaVector ^ mesh.C());

volScalarField hEnt = thermo.he();
volScalarField i_comp = 0.5*(magSqr(Urot) - magSqr(U)); // rothalpy complement

fvScalarMatrix EEqn
(
	fvm::div(phi, hEnt)
  - fvc::div(phi, i_comp)
  - fvm::laplacian(turbulence->alphaEff(), hEnt)
  - viscousDiss
);

EEqn.relax();

// --- Pressure Equation (pEqn.H modified) --- //
volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(rho)*fvc::flux(HbyA));

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

adjustPhi(phiHbyA, U, p);

fvScalarMatrix pEqn
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rhorAUf, p)
        );
		
pEqn.setReference(pressureControl.refCell(), pressureControl.refValue());
