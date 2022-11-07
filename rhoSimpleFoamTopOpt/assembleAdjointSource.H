densVec dJdvar(Eigen::VectorXd::Zero((3 + 2)*nCells + nFaces)); // Source term for adjoint problem (initializing with zeros)

thermo.T() = thermo.he()/thermo.Cp();
rho        = (p/thermo.he())*(thermo.Cp()/Rgas);

#include "objectiveFunction.H"
scalar Jcurr = J;

Info << "\nAssembling source vector for adjoint system..." << flush;

faceSet inOutFaces
(
	IOobject
	(
		"inOutFaces", 
		searchDir,
		mesh, 
		IOobject::MUST_READ, 
		IOobject::NO_WRITE
	)
);

const labelList& inOutFacesList = inOutFaces.toc();

forAll(inOutFacesList, ID)
{
	int faceID = inOutFacesList[ID];
	int row    = mesh.owner()[faceID];
		
	for(int uCmpt = 0; uCmpt < 3; uCmpt++)
	{
		// Perturb state
		if(uCmpt == 0) U[row] = U[row] + uScale*vector(stepJac, 0      , 0      );
		if(uCmpt == 1) U[row] = U[row] + uScale*vector(0      , stepJac, 0      );
		if(uCmpt == 2) U[row] = U[row] + uScale*vector(0      , 0      , stepJac);
		
		// Correct BC
		U.correctBoundaryConditions();
		p.correctBoundaryConditions();
		thermo.he().correctBoundaryConditions();
	
		thermo.T() = thermo.he()/thermo.Cp();
		rho        = (p/thermo.he())*(thermo.Cp()/Rgas);
		
		#include "objectiveFunction.H"
		
		// Recover state
		U = UState;
		p = pState;
		thermo.he() = hEntState;
		dJdvar[3*row + uCmpt] = (J - Jcurr)/stepJac;
	}
	
	{
		// Perturb state
		p[row] = p[row] + pScale*stepJac;
			
		// Correct BC
		U.correctBoundaryConditions();
		p.correctBoundaryConditions();
		thermo.he().correctBoundaryConditions();
		
		thermo.T() = thermo.he()/thermo.Cp();
		rho        = (p/thermo.he())*(thermo.Cp()/Rgas);
		
		#include "objectiveFunction.H"
	
		// Recover state
		U = UState;
		p = pState;
		thermo.he() = hEntState;
		dJdvar[3*nCells + row] = (J - Jcurr)/stepJac;
	}
	
	{	
		// Perturb state
		thermo.he()[row] = thermo.he()[row] + hScale*stepJac;
			
		// Correct BC
		U.correctBoundaryConditions();
		p.correctBoundaryConditions();
		thermo.he().correctBoundaryConditions();
		
		thermo.T() = thermo.he()/thermo.Cp();
		rho        = (p/thermo.he())*(thermo.Cp()/Rgas);
		
		#include "objectiveFunction.H"
		
		// Recover state
		U = UState;
		p = pState;
		thermo.he() = hEntState;
		dJdvar[4*nCells + row] = (J - Jcurr)/stepJac;
	}
}

Info << "completed!" << endl;
