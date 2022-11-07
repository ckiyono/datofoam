std::vector<std::vector<int>> cellConnect;       // Cell to cell connectivity list
std::vector<std::vector<int>> cellNeighbours;    // Cell to cell connectivity list level 1 only
std::vector<std::vector<int>> faceConnect;       // Face to face connectivity list
std::vector<std::vector<int>> cellToFaceConnect; // Cell to face connectivity list
std::vector<std::vector<int>> faceToCellConnect; // Face to cell connectivity list
std::vector<int> faceToCellNeighbourAMI;         // Face to cell connectivity list (only used when there is cyclicAMI)
std::vector<int> emptyVec;                       // Empty vector for general use
std::vector<int> m1List;                         // List with -1 to check presence of elements
string cyclicBC = "no";                          // To know if there is cyclicAM

for(int i = 0; i < nCells; i++) cellConnect.push_back(emptyVec);
for(int i = 0; i < nFaces; i++) faceConnect.push_back(emptyVec);
for(int i = 0; i < nCells; i++) cellNeighbours.push_back(emptyVec);
for(int i = 0; i < nFaces; i++) faceToCellNeighbourAMI.push_back(-1);
for(int faceID = 0; faceID < nFaces; faceID++) m1List.push_back(-1);

// Cell list level 1 (this is used to assemble cellConnect)
forAll(mesh.C(), cellL0)
{
	cellNeighbours[cellL0].push_back(cellL0);
	
	forAll(mesh.cellCells()[cellL0], cellL1)
	{
		int neigh1 = mesh.cellCells()[cellL0][cellL1];
		if(std::find(cellNeighbours[cellL0].begin(), cellNeighbours[cellL0].end(), neigh1) == cellNeighbours[cellL0].end())
			cellNeighbours[cellL0].push_back(neigh1);
	}
}

// Adding cells when there are cyclicAMI patches (to cellNeighbours)
scalar rotAngleDeg(readScalar(optimizationProperties.lookup("rotAngle")));
scalar rotAngle = rotAngleDeg*M_PI/180;

for(int i = 0; i < 2; i++)
{
	rotAngle = -rotAngle;
	
	forAll(mesh.boundary(), patchID)
	{
		if(mesh.boundary()[patchID].type()=="cyclicAMI")
		{
			cyclicBC = "yes";
			
			for(int faceCont = 0; faceCont < mesh.boundaryMesh()[patchID].size(); faceCont++)
			{	
				int faceID    = faceCont + mesh.boundaryMesh()[patchID].start();
				int cellOwner = mesh.owner()[faceID];
	
				point refPoint = mesh.Cf()[faceID] + (mesh.Cf()[faceID] - mesh.C()[cellOwner]);
				point rotPoint = point((refPoint.x())*Foam::cos(rotAngle) - (refPoint.y())*Foam::sin(rotAngle), (refPoint.x())*Foam::sin(rotAngle) + (refPoint.y())*Foam::cos(rotAngle), mesh.C()[cellOwner].z());
				
				int cellFound = mesh.findCell(rotPoint);
				
				if((cellFound >= 0) and (cellFound < nCells))
				{
					faceToCellNeighbourAMI[faceID] = cellFound;
					cellNeighbours[cellOwner].push_back(cellFound);
				}
			}
		}
	}
}

// Cell connectivity list (up to level 6)
Info << "\nCreating cell connectivity list (level 6)..." << endl;

forAll(mesh.C(), cellL0)
{
	for(unsigned int cellL1 = 0; cellL1 < cellNeighbours[cellL0].size(); cellL1++)
	{
		int neigh1 = cellNeighbours[cellL0][cellL1];
		if(m1List[neigh1] < 0)
		{
			cellConnect[cellL0].push_back(neigh1);
			m1List[neigh1] = 1;
		}
		
		for(unsigned int cellL2 = 0; cellL2 < cellNeighbours[neigh1].size(); cellL2++)
		{
			int neigh2 = cellNeighbours[neigh1][cellL2];
			if(m1List[neigh2] < 0)
			{
				cellConnect[cellL0].push_back(neigh2);
				m1List[neigh2] = 1;
			}
			
			for(unsigned int cellL3 = 0; cellL3 < cellNeighbours[neigh2].size(); cellL3++)
			{
				int neigh3 = cellNeighbours[neigh2][cellL3];
				if(m1List[neigh3] < 0)
				{
					cellConnect[cellL0].push_back(neigh3);
					m1List[neigh3] = 1;
				}

				for(unsigned int cellL4 = 0; cellL4 < cellNeighbours[neigh3].size(); cellL4++)
				{
					int neigh4 = cellNeighbours[neigh3][cellL4];
					if(m1List[neigh4] < 0)
					{
						cellConnect[cellL0].push_back(neigh4);
						m1List[neigh4] = 1;
					}
					
					for(unsigned int cellL5 = 0; cellL5 < cellNeighbours[neigh4].size(); cellL5++)
					{
						int neigh5 = cellNeighbours[neigh4][cellL5];
						if(m1List[neigh5] < 0)
						{
							cellConnect[cellL0].push_back(neigh5);
							m1List[neigh5] = 1;
						}
						
						for(unsigned int cellL6 = 0; cellL6 < cellNeighbours[neigh5].size(); cellL6++)
						{
							int neigh6 = cellNeighbours[neigh5][cellL6];
							if(m1List[neigh6] < 0)
							{
								cellConnect[cellL0].push_back(neigh6);
								m1List[neigh6] = 1;
							}
						}
					}
				}
			}
		}
	}
	
	if(cellL0%1000==0) Info << "\n...assessing cell: " << cellL0 << "..." << endl;
	std::fill(m1List.begin(), m1List.end(), -1);
}

Info << "\n...completed!" << endl;

// Face connectivity list
Info << "\nCreating face connectivity list..." << endl;

for(int faceI = 0; faceI < nFaces; faceI++)
{
	int cellOwner = mesh.owner()[faceI];
	
	for(unsigned int i = 0; i < cellConnect[cellOwner].size(); i++)
	{
		int cellN = cellConnect[cellOwner][i];
		forAll(mesh.cells()[cellN], faceL)
		{
			int faceID = mesh.cells()[cellN][faceL];
			if((m1List[faceID] < 0) and (faceID < nFaces))
			{
				faceConnect[faceI].push_back(faceID); // Face is pushed back if it is not in the list AND if it not empty (= smaller than nFaces)
				m1List[faceID] = 1;
			}
		}
	}
	
	if(faceI <= (mesh.neighbour().size()-1)) // Condition because not all faces have neighbour cells in polymesh
	{
		int cellNeighbour = mesh.neighbour()[faceI];
	
		if((cellOwner < cellNeighbour) and (cellNeighbour > 0) and (cellNeighbour < nCells))
		{
			for(unsigned int i = 0; i < cellConnect[cellNeighbour].size(); i++)
			{
				int cellN = cellConnect[cellNeighbour][i];		
				forAll(mesh.cells()[cellN], faceL)
				{
					int faceID = mesh.cells()[cellN][faceL];
					if((m1List[faceID] < 0) and (faceID < nFaces))
					{
						faceConnect[faceI].push_back(faceID); // Face is pushed back if it is not in the list AND if it not empty (= smaller than nFaces)
						m1List[faceID] = 1;
					}
				}
			}
		}
	}
	
	if((cyclicBC == "yes") and (faceToCellNeighbourAMI[faceI] >= 0))
	{
		int cellAMI = faceToCellNeighbourAMI[faceI];
	
		for(unsigned int i = 0; i < cellConnect[cellAMI].size(); i++)
		{
			int cellN = cellConnect[cellAMI][i];
			forAll(mesh.cells()[cellN], faceL)
			{
				int faceID = mesh.cells()[cellN][faceL];
				if((m1List[faceID] < 0) and (faceID < nFaces))
				{
					faceConnect[faceI].push_back(faceID); // Face is pushed back if it is not in the list AND if it not empty (= smaller than nFaces)
					m1List[faceID] = 1;
				}
			}
		}
	}
	
	if(faceI%1000==0) Info << "\n...assessing face: " << faceI << "..." << endl;
	std::fill(m1List.begin(), m1List.end(), -1);
}

Info << "\n...completed!" << endl;

// Preliminar cell group assembly
Info << "\nPreliminar cell group assembly..." << endl;

std::vector<std::vector<int>> preOrthoCellGroups; // Preliminar orthogonal group
std::vector<int> forbiddenCells;                  // Neighbour cells that do not allow to push back current cell
std::vector<int> notUsedCells;                    // Not used cells
forAll(mesh.C(), cellID) notUsedCells.push_back(cellID);

forAll(mesh.C(), groupID)
{
	preOrthoCellGroups.push_back(emptyVec); // Always include one empty vector
	
	for(unsigned int i = 0; i < notUsedCells.size(); i++)
	{
		int cellID = notUsedCells[i];
		
		if(preOrthoCellGroups[groupID].size() < 1)
		{
			preOrthoCellGroups[groupID].push_back(cellID);
			notUsedCells.erase(notUsedCells.begin() + i);
			i--;
			m1List[cellID] = 1;
			for(unsigned int j = 0; j < cellConnect[cellID].size(); j++) m1List[cellConnect[cellID][j]] = 1;
		}
	
		if(preOrthoCellGroups[groupID].size() > 0)
		{
			if(m1List[cellID] < 0)
			{
				preOrthoCellGroups[groupID].push_back(cellID);
				notUsedCells.erase(notUsedCells.begin() + i);
				i--;
				m1List[cellID] = 1;
				for(unsigned int j = 0; j < cellConnect[cellID].size(); j++) m1List[cellConnect[cellID][j]] = 1;
			}
		}
		if(notUsedCells.size() == 0) break;
	}
	std::fill(m1List.begin(), m1List.end(), -1);
	Info << "\n...assessing group: " << groupID << " / number of cells: " << preOrthoCellGroups[groupID].size() << "..." << endl;
	if(notUsedCells.size() == 0) break;
}

Info << "\n...completed!" << endl;

// Preliminar face group assembly
Info << "\nPreliminar face group assembly..." << endl;

std::vector<std::vector<int>> preOrthoFaceGroups; // Preliminar orthogonal group
std::vector<int> forbiddenFaces;                  // Neighbour facess that do not allow to push back current face
std::vector<int> notUsedFaces;                    // Not used faces
for(int faceID = 0; faceID < nFaces; faceID++) notUsedFaces.push_back(faceID);

for(int groupID = 0; groupID < nFaces; groupID++)
{
	preOrthoFaceGroups.push_back(emptyVec); // Always include one empty vector
	
	for(unsigned int i = 0; i < notUsedFaces.size(); i++)
	{
		int faceID = notUsedFaces[i];
		
		if(preOrthoFaceGroups[groupID].size() < 1)
		{
			preOrthoFaceGroups[groupID].push_back(faceID);
			notUsedFaces.erase(notUsedFaces.begin() + i);
			i--;
			m1List[faceID] = 1;
			for(unsigned int j = 0; j < faceConnect[faceID].size(); j++) m1List[faceConnect[faceID][j]] = 1;
		}
	
		if(preOrthoFaceGroups[groupID].size() > 0)
		{
			if(m1List[faceID] < 0)
			{			
				preOrthoFaceGroups[groupID].push_back(faceID);
				notUsedFaces.erase(notUsedFaces.begin() + i);
				i--;
				m1List[faceID] = 1;
				for(unsigned int j = 0; j < faceConnect[faceID].size(); j++) m1List[faceConnect[faceID][j]] = 1;
			}
		}
		
		if(notUsedFaces.size() == 0) break;
	}
	std::fill(m1List.begin(), m1List.end(), -1);
	Info << "\n...assessing group: " << groupID << " / number of faces: " << preOrthoFaceGroups[groupID].size() << "..." << endl;
	if(notUsedFaces.size() == 0) break;
}

Info << "\n...completed!" << endl;

// Assembly of final groups for state and design variables
Info << "\nAssembling groups for state and design variables..." << flush;

int preOrthoCellSize = (int) preOrthoCellGroups.size(); // Number of cell groups

for(int i = 0; i < 5*preOrthoCellSize; i++) orthoCols.push_back(emptyVec);

for(unsigned int groupID = 0; groupID < preOrthoCellGroups.size(); groupID++)
{
	for(unsigned int elemID = 0; elemID < preOrthoCellGroups[groupID].size(); elemID++)
	{
		orthoCols[3*groupID + 0].push_back(3*preOrthoCellGroups[groupID][elemID] + 0);                     // Ux (= 3*cell ID + 0)
		orthoCols[3*groupID + 1].push_back(3*preOrthoCellGroups[groupID][elemID] + 1);                     // Uy (= 3*cell ID + 1)
		orthoCols[3*groupID + 2].push_back(3*preOrthoCellGroups[groupID][elemID] + 2);                     // Uz (= 3*cell ID + 2)
		orthoCols[3*preOrthoCellSize + groupID].push_back(3*nCells + preOrthoCellGroups[groupID][elemID]); //  p (= 3*nCells + ID)
		orthoCols[4*preOrthoCellSize + groupID].push_back(4*nCells + preOrthoCellGroups[groupID][elemID]); //  h (= 4*nCells + ID)
	}
}

// Insert col IDs for face state variables
for(unsigned int groupID = 0; groupID < preOrthoFaceGroups.size(); groupID++)
{
	orthoCols.push_back(emptyVec);
	for(unsigned int elemID = 0; elemID < preOrthoFaceGroups[groupID].size(); elemID++)
	{
		orthoCols[5*preOrthoCellSize + groupID].push_back(5*nCells + preOrthoFaceGroups[groupID][elemID]); // faceID (= 5*nCells + local faceID)
	}
}

Info << "completed!" << endl;
Info << "\nTotal number of groups: " << orthoCols.size() << endl;

// Clear lists to re-use to assemble colRowNZ
cellConnect.clear();
faceConnect.clear();

for(int i = 0; i < nCells; i++) cellConnect.push_back(emptyVec);
for(int i = 0; i < nCells; i++) cellToFaceConnect.push_back(emptyVec);
for(int i = 0; i < nFaces; i++) faceConnect.push_back(emptyVec);
for(int i = 0; i < nFaces; i++) faceToCellConnect.push_back(emptyVec);

// Cell connectivity list (up to level 3)
Info << "\nRe-creating cell connectivity list (level 3)..." << endl;

forAll(mesh.C(), cellL0)
{
	for(unsigned int cellL1 = 0; cellL1 < cellNeighbours[cellL0].size(); cellL1++)
	{
		int neigh1 = cellNeighbours[cellL0][cellL1];
		if(m1List[neigh1] < 0)
		{
			cellConnect[cellL0].push_back(neigh1);
			m1List[neigh1] = 1;
		}
		
		for(unsigned int cellL2 = 0; cellL2 < cellNeighbours[neigh1].size(); cellL2++)
		{
			int neigh2 = cellNeighbours[neigh1][cellL2];
			if(m1List[neigh2] < 0)
			{
				cellConnect[cellL0].push_back(neigh2);
				m1List[neigh2] = 1;
			}
			
			for(unsigned int cellL3 = 0; cellL3 < cellNeighbours[neigh2].size(); cellL3++)
			{
				int neigh3 = cellNeighbours[neigh2][cellL3];
				if(m1List[neigh3] < 0)
				{
					cellConnect[cellL0].push_back(neigh3);
					m1List[neigh3] = 1;
				}
			}
		}
	}
	
	if(cellL0%1000==0) Info << "\n...assessing cell: " << cellL0 << "..." << endl;
	std::fill(m1List.begin(), m1List.end(), -1);
}

Info << "\n...completed!" << endl;

// Face connectivity list
Info << "\nRe-creating face connectivity list..." << endl;

for(int faceI = 0; faceI < nFaces; faceI++)
{
	int cellOwner = mesh.owner()[faceI];
	
	for(unsigned int i = 0; i < cellConnect[cellOwner].size(); i++)
	{
		int cellN = cellConnect[cellOwner][i];
		forAll(mesh.cells()[cellN], faceL)
		{
			int faceID = mesh.cells()[cellN][faceL];
			if((m1List[faceID] < 0) and (faceID < nFaces))
			{
				faceConnect[faceI].push_back(faceID); // Face is pushed back if it is not in the list AND if it not empty (= smaller than nFaces)
				m1List[faceID] = 1;
			}
		}
	}
	
	if(faceI <= (mesh.neighbour().size()-1)) // Condition because not all faces have neighbour cells in polymesh
	{
		int cellNeighbour = mesh.neighbour()[faceI];
	
		if((cellOwner < cellNeighbour) and (cellNeighbour > 0) and (cellNeighbour < nCells))
		{
			for(unsigned int i = 0; i < cellConnect[cellNeighbour].size(); i++)
			{
				int cellN = cellConnect[cellNeighbour][i];		
				forAll(mesh.cells()[cellN], faceL)
				{
					int faceID = mesh.cells()[cellN][faceL];
					if((m1List[faceID] < 0) and (faceID < nFaces))
					{
						faceConnect[faceI].push_back(faceID); // Face is pushed back if it is not in the list AND if it not empty (= smaller than nFaces)
						m1List[faceID] = 1;
					}
				}
			}
		}
	}
	
	if((cyclicBC == "yes") and (faceToCellNeighbourAMI[faceI] >= 0))
	{
		int cellAMI = faceToCellNeighbourAMI[faceI];
	
		for(unsigned int i = 0; i < cellConnect[cellAMI].size(); i++)
		{
			int cellN = cellConnect[cellAMI][i];
			forAll(mesh.cells()[cellN], faceL)
			{
				int faceID = mesh.cells()[cellN][faceL];
				if((m1List[faceID] < 0) and (faceID < nFaces))
				{
					faceConnect[faceI].push_back(faceID); // Face is pushed back if it is not in the list AND if it not empty (= smaller than nFaces)
					m1List[faceID] = 1;
				}
			}
		}
	}
	
	if(faceI%1000==0) Info << "\n...assessing face: " << faceI << "..." << endl;
	std::fill(m1List.begin(), m1List.end(), -1);
}

Info << "\n...completed!" << endl;

// Assembly of cellToFaceConnect
Info << "\nCell to face group assembly..." << flush;

forAll(mesh.C(), cellID)
{
	for(unsigned int i = 0; i < cellConnect[cellID].size(); i++)
	{
		int cellN = cellConnect[cellID][i];		
		forAll(mesh.cells()[cellN], faceL)
		{
			int faceID = mesh.cells()[cellN][faceL];
			if((std::find(cellToFaceConnect[cellID].begin(), cellToFaceConnect[cellID].end(), faceID) == cellToFaceConnect[cellID].end()) and (faceID < nFaces))
				cellToFaceConnect[cellID].push_back(faceID);
		}
	}
}

Info << "completed!" << endl;

// Assembly of faceToCellConnect
Info << "\nFace to cell group assembly..." << flush;

for(int faceI = 0; faceI < nFaces; faceI++)
{
	int cellOwner = mesh.owner()[faceI];
	
	for(unsigned int i = 0; i < cellConnect[cellOwner].size(); i++)
	{
		int cellN = cellConnect[cellOwner][i];		
		if(std::find(faceToCellConnect[faceI].begin(), faceToCellConnect[faceI].end(), cellN) == faceToCellConnect[faceI].end())
			faceToCellConnect[faceI].push_back(cellN);
	}
	
	if(faceI < nPhi)
	{
		int cellNeighbour = mesh.neighbour()[faceI];
		for(unsigned int i = 0; i < cellConnect[cellNeighbour].size(); i++)
		{
			int cellN = cellConnect[cellNeighbour][i];		
			if(std::find(faceToCellConnect[faceI].begin(), faceToCellConnect[faceI].end(), cellN) == faceToCellConnect[faceI].end())
				faceToCellConnect[faceI].push_back(cellN);
		}
	}

	if((cyclicBC == "yes") and (faceToCellNeighbourAMI[faceI] >= 0))
	{
		int cellAMI = faceToCellNeighbourAMI[faceI];
	
		for(unsigned int i = 0; i < cellConnect[cellAMI].size(); i++)
		{
			int cellN = cellConnect[cellAMI][i];
			if(std::find(faceToCellConnect[faceI].begin(), faceToCellConnect[faceI].end(), cellN) == faceToCellConnect[faceI].end())
				faceToCellConnect[faceI].push_back(cellN);
		}
	}
}

Info << "completed!" << endl;

// Assembly of possible non-zeros for state and design variables
Info << "\nAssembling groups of non-zeros for state and design variables..." << flush;

for(int i = 0; i < (nCells + nFaces); i++) colRowNZ.push_back(emptyVec);

for(unsigned int cellID = 0; cellID < cellConnect.size(); cellID++)
{
	for(unsigned int cellN = 0; cellN < cellConnect[cellID].size(); cellN++)
	{
		colRowNZ[cellID].push_back(cellConnect[cellID][cellN]);
	}
	
	for(unsigned int faceN = 0; faceN < cellToFaceConnect[cellID].size(); faceN++)
	{
		colRowNZ[cellID].push_back(nCells + cellToFaceConnect[cellID][faceN]);
	}
}

for(int faceID = 0; faceID < nFaces; faceID++)
{
	for(unsigned int faceN = 0; faceN < faceConnect[faceID].size(); faceN++)
	{
		colRowNZ[nCells + faceID].push_back(nCells + faceConnect[faceID][faceN]);
	}
	for(unsigned int cellN = 0; cellN < faceToCellConnect[faceID].size(); cellN++)
	{
		colRowNZ[nCells + faceID].push_back(faceToCellConnect[faceID][cellN]);
	}
}

Info << "completed!" << endl;

// Clear lists
m1List.clear();
cellConnect.clear();
faceConnect.clear();
cellToFaceConnect.clear();
faceToCellConnect.clear();
