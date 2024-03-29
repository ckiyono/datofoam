/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    rhoSimpleFoam

Group
    grpCompressibleSolvers

Description
    Steady-state solver for compressible turbulent flow.

\*---------------------------------------------------------------------------*/

#include <iostream> // file manipulation
#include <fstream>  // file manipulation
#include <sstream>  // file manipulation

#include "fvCFD.H"
#include "fluidThermo.H"
#include "turbulentFluidThermoModel.H"
#include "simpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"

#include "cellSet.H"       // to include sets of cells
#include "faceSet.H"       // to include sets of cells
#include "dynamicFvMesh.H" // to include adaptive mesh refinement

#include "IOmanip.H"       // header for setprecision 
#include "Eigen/Eigen"     // header for eigen libs

#include "reverseBytes.H"  // Function to reverse the byte order

typedef Eigen::Triplet<double>                      triplet;  // triplets to sparse matrix assembly
typedef Eigen::SparseMatrix<double,Eigen::RowMajor> spMatrix; // sparse matrix in Eigen (row major storage format)
typedef Eigen::VectorXd                             densVec;  // dense vector in Eigen

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for compressible turbulent flow."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    // #include "createMesh.H"
	#include "createDynamicFvMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	// Initial definitions
	#include "initialDefinitions.H"
	
	// Create groups
	#include "orthoGroups.H"

    for(int optCont = 0; optCont <= nOpt; optCont++)
	{
		// nOptTemp has current optimization iterations, so the material model uses qGray correctly
		std::ofstream nOptTempFile;
		nOptTempFile.open("nOptTemp.txt");
		nOptTempFile << optCont;
		nOptTempFile.close();
		
		// --- Adaptive mesh refinement
		#include "AMR.H"
		
		// --- Primal calculation
		#include "solveAndReadFlow.H"
		
		// --- Adjoint calculation
		#include "assembleStateJacobian.H"
		#include "assembleAdjointSource.H"
		#include "solveAdjoint.H"
		#include "assembleDesVarJacobian.H"
		#include "assembleAdjointFields.H"
		#include "stoppingCriteria.H"
	}
	
	system("rm -f nOptTemp.txt");
	
	Info << "\nOptimization completed!\n" << endl;
	
	runTime.printExecutionTime(Info);
	
    return 0;
}


// ************************************************************************* //
