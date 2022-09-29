// Initial guess for adjoint (can be enabled)
oldAdj.resize((3 + 2)*nCells + nFaces); // Redefine size
oldAdj.setZero();                       // Zero for all components

/* The 1st 'if' reads a txt file structured as:
/  nGroups
/  nElems elemVal elemVal elemVal elemVal
/  nElems elemVal elemVal
/  nElems elemVal elemVal elemVal elemVal elemVal
/
/  nGroups
/  nElems elemVal elemVal 
/  nElems elemVal elemVal elemVal elemVal
/  nElems elemVal
/ ...
*/

// When start to use OpenFOAM thermophysical models, use:
// thermo.correct();
// rho = thermo.rho();
// in the place of rho update below (also change in other files used for adjoint assembly)
dimensionedScalar pRef  ("pRef", dimMass*dimLength/(dimTime*dimTime*dimLength*dimLength), readScalar(thermophysicalProperties.lookup("pRef")));
dimensionedScalar TRef  ("TRef", dimTemperature, readScalar(thermophysicalProperties.lookup("TRef")));
dimensionedScalar molWeight("molWeight", dimensionSet(1,0,0,0,-1,0,0), readScalar(thermophysicalProperties.subDict("mixture").subDict("specie").lookup("molWeight")));
dimensionedScalar RUniversal("RUniversal", dimensionSet(1,2,-2,-1,-1,0,0), (Foam::constant::thermodynamic::RR));

dimensionedScalar Rgas = RUniversal/molWeight;

thermo.T() = thermo.he()/thermo.Cp();
rho        = (p/thermo.he())*(thermo.Cp()/Rgas);

std::fstream orthoColsFile;
orthoColsFile.open("orthoCols.txt");

if (orthoColsFile.is_open())
{
	Info << "\nOrthogonal groups found! Reading..." << flush;
	
	std::vector<int> listVal; // List with values to be read
	string nGroups;           // Number of groups in one group vector
	string nElems;            // Number of elements in one line
	string elemVals;          // Variable to store elements read
	string skipLine;          // String to go to next line
	
	for(int group = 0; group < 2; group++)
	{
		getline(orthoColsFile, nGroups, '\n');
		int nGroupsCont = std::stoi(nGroups);
		
		for(int k = 0; k < nGroupsCont; k++)
		{
			getline(orthoColsFile, nElems, ' ');
			int nElemsCont  = std::stoi(nElems);
			
			for(int i = 0; i < (nElemsCont-1); i++)
			{
				getline(orthoColsFile, elemVals, ' ');
				listVal.push_back(std::stoi(elemVals));
			}
			
			getline(orthoColsFile, elemVals, '\n');
			listVal.push_back(std::stoi(elemVals));
			
			if(group == 0) orthoCols.push_back(listVal);
			if(group == 1) colRowNZ.push_back(listVal);
			
			listVal.clear();
		}
		getline(orthoColsFile, skipLine); // Go to next line
	}
	Info << "completed!\n"<< endl;
}

else
{
	Info << "\nNo orthogonal groups found...assembling groups of orthogonal columns for state variables..." << endl;

	{
		#include "createGroups.H"
	}

	Info << "\nOrthogonal groups for state and design variables completed!\n" << endl;
	
	std::ofstream orthoColsFileWrite;
	orthoColsFileWrite.open("orthoCols.txt");
	
	// Write orthoCols
    orthoColsFileWrite << orthoCols.size();
	for(unsigned int group = 0; group < orthoCols.size(); group++)
	{
		orthoColsFileWrite << "\n" << orthoCols[group].size() << " ";
		
		for(unsigned int i = 0; i < orthoCols[group].size(); i++)
			orthoColsFileWrite << orthoCols[group][i] << " ";
	}
	
	orthoColsFileWrite << "\n";
	
	// Write colRowNZ
	orthoColsFileWrite << "\n" << colRowNZ.size();
	for(unsigned int group = 0; group < colRowNZ.size(); group++)
	{
		orthoColsFileWrite << "\n" << colRowNZ[group].size() << " ";
		
		for(unsigned int i = 0; i < colRowNZ[group].size(); i++)
			orthoColsFileWrite << colRowNZ[group][i] << " ";
	}
	
	orthoColsFileWrite << "\n";
	
	orthoColsFileWrite.close();	
}
