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

#include "fvCFD.H"
#include "fluidThermo.H"
#include "turbulentFluidThermoModel.H"
#include "simpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"

#include "cellSet.H" // to include sets of cells
#include "faceSet.H"       // to include sets of cells

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
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
	
	#include "finiteDiferenceCheck.H"     // FD verification (can be enabled)
	#include "updateDesVar.H"             // Update rhoDes
	#include "../interpolationFunction.H" // Apply interpolation function

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Pressure-velocity SIMPLE corrector
        #include "UEqn.H"
        #include "EEqn.H"

        if (simple.consistent())
        {
            #include "pcEqn.H"
        }
        else
        {
            #include "pEqn.H"
        }

        turbulence->correct();

        runTime.write();
		if(runTime.write())
		{
			rhoDes.write();
			alphaPen.write();
			designDomain.write();
			if(optConst == "velocProfile") unifDomain.write();
			if((optConst == "pRatio") or (objFunc == "impellerEff"))
			{
				pRatioOut.write();
				pRatioIn.write();
			}
		}

        runTime.printExecutionTime(Info);
		
		#include "../objectiveFunction.H"
		#include "objFuncWrite.H"
		
		Info << "Volume fraction: " << volFrac*100 << "%" << endl;		
		Info << "Objective function value: " << J << "\n" << endl;
    }
		
	Info<< "End\n" << endl;
	
    return 0;
}


// ************************************************************************* //
