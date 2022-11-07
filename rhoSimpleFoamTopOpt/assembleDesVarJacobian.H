// Current state fields (used to reset state)
const volScalarField rhoDesState = rhoDes;

Info << "Assembling design variable Jacobian matrix..." << flush;
dRdvar.data().squeeze();                        // Clear matrix
dRdvar.resize(nCells, (3 + 2)*nCells + nFaces); // Redefine size

// Perturb design variables
for(unsigned int group = 0; group < orthoCols.size(); group++)
{
	for(unsigned int colID = 0; colID < orthoCols[group].size(); colID++)
	{
		int cellID     = orthoCols[group][colID]; // Column to perturb
		if((cellID%3 == 0) and (cellID < 3*nCells))
		{
			cellID = cellID/3;
			rhoDes[cellID] = rhoDes[cellID] + stepJac;
		}
	}
	
	#include "interpolationFunction.H"
	
	#include "assembleResiduals.H"
		
	UresPert   = ((UEqn & U) - (Ures))/stepJac;
	pResPert   = ((pEqn & p) - (pRes))/stepJac;
	hResPert   = ((EEqn & hEnt) - (hRes))/stepJac;
	phiResPert = (phiHbyA + pEqn.flux() - phi - phiRes)/stepJac;
		
	// Recover state
	rhoDes = rhoDesState;
	
	// Place Jacobian matrix elements (ALREADY TRANSPOSED)
	for(unsigned int colID = 0; colID < orthoCols[group].size(); colID++)
	{
		int col = orthoCols[group][colID]; // Perturbed column
		
		if((col%3 == 0) and (col < 3*nCells))
		{
			col = col/3;
			
			for(unsigned int rowID = 0; rowID < colRowNZ[col].size(); rowID++)
			{
				int row = colRowNZ[col][rowID];
				
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
}

Info << "completed!\n" << endl;

dRdvar.setFromTriplets(jacList.begin(), jacList.end()); // Construct sparse matrix using triplets
dRdvar.makeCompressed();                                // Eliminate zeros

jacList.clear(); // Clear list

// Print matrix in ordered manner for visualization
if(infoMatrix == "yes")
{
	int sizeM = (3 + 2)*nCells + nFaces;
	Info << "Matrix: " << endl;
	for(int row = 0; row < sizeM; row++)
	{
		for(int col = 0; col < nCells; col++)
		{
			Info << " " << dRdvar.coeffRef(row,col) << " ";
		}
		Info << endl;
	}
}
