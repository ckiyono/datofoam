thermo.T() = thermo.he()/thermo.Cp();
rho        = (p/thermo.he())*(thermo.Cp()/Rgas);

// Current state fields (used to reset state)
const volVectorField UState       = U;
const volScalarField pState       = p;
const volScalarField hEntState    = thermo.he();
const surfaceScalarField phiState = phi;

Info << "\nStarting sensitivity calculation!" << endl;

Info << "\nAssembling state Jacobian matrix..." << flush;

spMatrix dRdvar((3 + 2)*nCells + nFaces, (3 + 2)*nCells + nFaces); // Jacobian matrix
std::vector<triplet> jacList;                                      // List to assemble sparse Jacobian (state) matrix

#include "assembleResiduals.H"

volVectorField Ures     = UEqn & U;                          // Residual for current primal solution
volVectorField UresPert = Ures;                              // Residual for perturbed field

volScalarField pRes     = pEqn & p;                          // Residual for current primal solution
volScalarField pResPert = pRes;                              // Residual for perturbed field

volScalarField hRes     = EEqn & hEnt;                       // Residual for current primal solution
volScalarField hResPert = hRes;                              // Residual for perturbed field

surfaceScalarField phiRes     = phiHbyA + pEqn.flux() - phi; // Residual for current primal solution
surfaceScalarField phiResPert = phi;                         // Residual for perturbed field

// Perturb state variables
for(unsigned int group = 0; group < orthoCols.size(); group++)
{
	for(unsigned int colID = 0; colID < orthoCols[group].size(); colID++)
	{
		int col = orthoCols[group][colID]; // Column to perturb
		
		if(col < 3*nCells) // Perturb velocities
		{
			int cellID = col/3;
			
			if( (col%3) == 0 ) U[cellID] = U[cellID] + uScale*vector(stepJac, 0      , 0      );
			if( (col%3) == 1 ) U[cellID] = U[cellID] + uScale*vector(      0, stepJac, 0      );
			if( (col%3) == 2 ) U[cellID] = U[cellID] + uScale*vector(      0, 0      , stepJac);
		}
		
		if( (col >= 3*nCells) and (col < 4*nCells) ) // Perturb pressure
		{
			int cellID = col - 3*nCells;
			p[cellID]  = p[cellID] + pScale*stepJac;
		}
		
		if( (col >= 4*nCells) and (col < 5*nCells) ) // Perturb enthalpy
		{
			int cellID          = col - 4*nCells;
			thermo.he()[cellID] = thermo.he()[cellID] + hScale*stepJac;
		}
		
		if( (col >= 5*nCells) and (col < (5*nCells + nPhi)) ) // Perturb internal phi
		{
			int phiID  = col - 5*nCells;
			phi[phiID] = phi[phiID] + phiScale*stepJac;
		}
		
		if(col >= (5*nCells) + nPhi) // Perturb boundary phi
		{
			int phiID        = col - (5*nCells + nPhi);
			int phiBoundCont = 0;
			
			forAll(mesh.boundary(), patchI)
			{
				if( (phiID >= phiBoundCont) and (phiID < (phiBoundCont + mesh.boundary()[patchI].size())) )
					phi.boundaryFieldRef()[patchI][phiID - phiBoundCont] = phi.boundaryField()[patchI][phiID - phiBoundCont] + phiScale*stepJac;
				
				phiBoundCont = phiBoundCont + mesh.boundary()[patchI].size();
			}
		}
	}
	
	// Correct BC
	U.correctBoundaryConditions();
	p.correctBoundaryConditions();
	thermo.he().correctBoundaryConditions();
	
	thermo.T() = thermo.he()/thermo.Cp();
	rho        = (p/thermo.he())*(thermo.Cp()/Rgas);
	
	#include "assembleResiduals.H"
		
	UresPert   = ((UEqn & U) - (Ures))/stepJac;
	pResPert   = ((pEqn & p) - (pRes))/stepJac;
	hResPert   = ((EEqn & hEnt) - (hRes))/stepJac;		
	phiResPert = (phiHbyA + pEqn.flux() - phi - phiRes)/stepJac;
		
	// Recover state
	U           = UState;
	p           = pState;
	thermo.he() = hEntState;
	phi         = phiState;
	
	thermo.T() = thermo.he()/thermo.Cp();
	rho        = (p/thermo.he())*(thermo.Cp()/Rgas);
	
	// Place Jacobian matrix elements (ALREADY TRANSPOSED)
	for(unsigned int colID = 0; colID < orthoCols[group].size(); colID++)
	{
		int col = orthoCols[group][colID]; // Perturbed column
		int col2 = col;
		
		if((col < 3*nCells) and (col%3==0))        col2 = (col-0)/3;
		if((col < 3*nCells) and (col%3==1))        col2 = (col-1)/3;
		if((col < 3*nCells) and (col%3==2))        col2 = (col-2)/3;
		if((col >= 3*nCells) and (col < 4*nCells)) col2 = col - 3*nCells;
		if((col >= 4*nCells) and (col < 5*nCells)) col2 = col - 4*nCells;
		if(col >= 5*nCells)                        col2 = col - 4*nCells;
		
		for(unsigned int rowID = 0; rowID < colRowNZ[col2].size(); rowID++)
		{
			int row = colRowNZ[col2][rowID];
			
			if(row < nCells) // Collect velocity, pressure and enthalpy residuals
			{
				int cellID = row;
				
				// Velocity, pressure and enthalpy residuals integrated in cell volume
				double URX = UresPert.internalField()[cellID][0]*mesh.V()[cellID];
				double URY = UresPert.internalField()[cellID][1]*mesh.V()[cellID];
				double URZ = UresPert.internalField()[cellID][2]*mesh.V()[cellID];
				double pR = pResPert.internalField()[cellID]*mesh.V()[cellID];
				double hR = hResPert.internalField()[cellID]*mesh.V()[cellID];
				
				if(abs(URX) > spsTol) jacList.push_back(triplet(col, 3*row + 0, URX));
				if(abs(URY) > spsTol) jacList.push_back(triplet(col, 3*row + 1, URY));
				if(abs(URZ) > spsTol) jacList.push_back(triplet(col, 3*row + 2, URZ));
				if(abs(pR) > spsTol)  jacList.push_back(triplet(col, 3*nCells + row, pR));
				if(abs(hR) > spsTol)  jacList.push_back(triplet(col, 4*nCells + row, hR));
			}
			
			if( (row >= nCells) and (row < (nCells + nPhi)) ) // Collect internal phi residuals
			{
				row        = row - nCells;
				int phiID  = row;
				
				// Phi residual in internal faces
				double phiR = phiResPert[phiID];
				if(abs(phiR) > spsTol) jacList.push_back(triplet(col, 5*nCells + row, phiR));
			}
						
			if(row >= (nCells + nPhi)) // Collect boundary phi residuals
			{
				row              = row - nCells - nPhi;
				int phiID        = row;
				int phiBoundCont = 0;
				
				forAll(mesh.boundary(), patchI)
				{
					if( (phiID >= phiBoundCont) and (phiID < (phiBoundCont + mesh.boundary()[patchI].size())) )
					{
						// Phi residual in boundary faces
						double phiR = phiResPert.boundaryField()[patchI][phiID - phiBoundCont];
						if(abs(phiR) > spsTol) jacList.push_back(triplet(col, 5*nCells + nPhi + row, phiR));
					}
					
					phiBoundCont = phiBoundCont + mesh.boundary()[patchI].size();
				}
			}
		}
	}
}

Info << "completed!" << endl;

dRdvar.setFromTriplets(jacList.begin(), jacList.end()); // Construct sparse matrix using triplets
dRdvar.makeCompressed();                                // Eliminate zeros

jacList.clear(); // Clear list

// Print matrix in ordered manner for visualization
word infoMatrix = word(optimizationProperties.lookup("infoMatrix"));
if(infoMatrix == "yes")
{
	int sizeM = (3 + 2)*nCells + nFaces;
	Info << "Matrix: " << endl;
	for(int row = 0; row < sizeM; row++)
	{
		for(int col = 0; col < sizeM; col++)
		{
			Info << " " << dRdvar.coeffRef(row,col) << " ";
		}
		Info << endl;
	}
}
