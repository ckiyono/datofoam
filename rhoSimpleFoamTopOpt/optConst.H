densVec dCdvar(Eigen::VectorXd::Zero((3 + 2)*nCells + nFaces)); // Source term for adjoint problem (initializing with zeros)

if(optConst == "velocProfile")
{
	Info << "Optimization constraint: velocity profile!" << endl;

	scalar uProf = 0;
	vector uRef (vector(optimizationProperties.lookup("uRef")));  // Velocity reference profile

	faceSet unifDomain
	(
		IOobject
		(
			"unifDomain", 
			searchDir,
			mesh, 
			IOobject::MUST_READ, 
			IOobject::NO_WRITE
		)
	);

	const labelList& unifDomainList = unifDomain.toc();

	forAll(unifDomainList, ID)
	{
		int faceID = unifDomainList[ID];
		vector uCalc = U[mesh.owner()[faceID]] - uRef;
		uProf = uProf + (uCalc & uCalc)*mag(mesh.Sf()[faceID]);
	}

	Info << "\nUniformity profile function value: " << uProf << endl;

	scalar uProfCurr = uProf;
	uProf = 0;

	forAll(unifDomainList, listID)
	{
		int row = mesh.owner()[unifDomainList[listID]];
		
		for(int uCmpt = 0; uCmpt < 3; uCmpt++)
		{
			// Perturb state
			if(uCmpt == 0) U[row] = U[row] + uScale*vector(stepJac, 0      , 0      );
			if(uCmpt == 1) U[row] = U[row] + uScale*vector(0      , stepJac, 0      );
			if(uCmpt == 2) U[row] = U[row] + uScale*vector(0      , 0      , stepJac);
			
			U.correctBoundaryConditions(); // Correct BC
					
			forAll(unifDomainList, ID)
			{
				int faceID = unifDomainList[ID];
				vector uCalc = U[mesh.owner()[faceID]] - uRef;
				uProf = uProf + (uCalc & uCalc)*mag(mesh.Sf()[faceID]);
			}
			
			U = UState; // Recover state
			dCdvar[3*row + uCmpt] = (uProf - uProfCurr)/stepJac;
			uProf = 0;
		}
	}
}

if(optConst == "pRatio")
{
	Info << "Optimization constraint: pressure ratio!" << endl;
	
	scalar pRatioConst = 0;
	scalar pOutCalc    = 0;
	scalar pInCalc     = 0;
	scalar areaOut     = 0;
	scalar areaIn      = 0;
	scalar pRatio (readScalar(optimizationProperties.lookup("pRatio"))); // Desired pressure ratio

	faceSet pRatioOut
	(
		IOobject
		(
			"pRatioOut", 
			searchDir,
			mesh, 
			IOobject::MUST_READ, 
			IOobject::NO_WRITE
		)
	);
	
	faceSet pRatioIn
	(
		IOobject
		(
			"pRatioIn", 
			searchDir,
			mesh, 
			IOobject::MUST_READ, 
			IOobject::NO_WRITE
		)
	);

	const labelList& pRatioOutList = pRatioOut.toc();
	const labelList& pRatioInList  = pRatioIn.toc();

	forAll(mesh.boundary(), patchI)
	{
		forAll(mesh.boundary()[patchI], faceI)
		{
			if(mesh.boundary()[patchI].name() == "outlet")
			{	
				pOutCalc = pOutCalc + p.boundaryField()[patchI][faceI]*mag(mesh.Sf().boundaryField()[patchI][faceI]);
				areaOut = areaOut + mag(mesh.Sf().boundaryField()[patchI][faceI]);
			}
			
			if(mesh.boundary()[patchI].name() == "inlet")
			{
				pInCalc = pInCalc + p.boundaryField()[patchI][faceI]*mag(mesh.Sf().boundaryField()[patchI][faceI]);
				areaIn = areaIn + mag(mesh.Sf().boundaryField()[patchI][faceI]);
			}
		}
	}
	
	pRatioConst = (pOutCalc/areaOut) / (pInCalc/areaIn) - pRatio;

	Info << "\nPressure ratio value: " << pRatioConst + pRatio << endl;

	scalar pRatioConstCurr = pRatioConst;
	pOutCalc = 0;
	pInCalc  = 0;

	forAll(pRatioOutList, listID)
	{
		int row = mesh.owner()[pRatioOutList[listID]];
		
		p[row] = p[row] + pScale*stepJac;
		p.correctBoundaryConditions(); // Correct BC
		
		forAll(mesh.boundary(), patchI)
		{
			forAll(mesh.boundary()[patchI], faceI)
			{
				if(mesh.boundary()[patchI].name() == "outlet")
				{	
					pOutCalc = pOutCalc + p.boundaryField()[patchI][faceI]*mag(mesh.Sf().boundaryField()[patchI][faceI]);
				}
				
				if(mesh.boundary()[patchI].name() == "inlet")
				{
					pInCalc = pInCalc + p.boundaryField()[patchI][faceI]*mag(mesh.Sf().boundaryField()[patchI][faceI]);
				}
			}
		}
		
		pRatioConst = (pOutCalc/areaOut) / (pInCalc/areaIn) - pRatio;
		
		p = pState; // Recover state
		dCdvar[3*nCells + row] = (pRatioConst - pRatioConstCurr)/stepJac;
		pOutCalc = 0;
		pInCalc  = 0;
	
	}
	
	forAll(pRatioInList, listID)
	{
		int row = mesh.owner()[pRatioInList[listID]];
		
		p[row] = p[row] + pScale*stepJac;
		p.correctBoundaryConditions(); // Correct BC
		
		forAll(mesh.boundary(), patchI)
		{
			forAll(mesh.boundary()[patchI], faceI)
			{
				if(mesh.boundary()[patchI].name() == "outlet")
				{	
					pOutCalc = pOutCalc + p.boundaryField()[patchI][faceI]*mag(mesh.Sf().boundaryField()[patchI][faceI]);
				}
				
				if(mesh.boundary()[patchI].name() == "inlet")
				{
					pInCalc = pInCalc + p.boundaryField()[patchI][faceI]*mag(mesh.Sf().boundaryField()[patchI][faceI]);
				}
			}
		}
		
		pRatioConst = (pOutCalc/areaOut) / (pInCalc/areaIn) - pRatio;
		
		p = pState; // Recover state
		dCdvar[3*nCells + row] = (pRatioConst - pRatioConstCurr)/stepJac;
		pOutCalc = 0;
		pInCalc  = 0;
	
	}
}


if((optConst == "velocProfile") or (optConst == "pRatio"))
{
	adjConstSolve = oldAdjConst;
	
	// Writing initial guess for constraint
	std::fstream file_oldAdjConst("oldAdj.dat", std::ios::binary | std::ios::in | std::ios::out | std::ios::trunc);
	
	if(!file_oldAdjConst.is_open()) cout << "\n\nError opening the file!\n\n";
	
	else
	{
		int VEC_CLASS_ID_BE = __builtin_bswap32(VEC_CLASS_ID);
		int nRows_BE        = __builtin_bswap32(nRows);
		
		file_oldAdjConst.write(reinterpret_cast<char*>(&VEC_CLASS_ID_BE), sizeof(int));
		file_oldAdjConst.write(reinterpret_cast<char*>(&nRows_BE), sizeof(int));
		
		for(int i = 0; i < nRows; i++)
		{
			double vecVal = dCdvar(i);
			char* vecValChar = reinterpret_cast<char*>(&vecVal);
			double vecVal_BE = reverseBytes(vecValChar);
			file_oldAdjConst.write(reinterpret_cast<char*>(&vecVal_BE), sizeof(double));
		}
		
		file_oldAdjConst.close();
	}
	
	// Writing source vector for constraint
	std::fstream file_dCdvar("dJdvar.dat", std::ios::binary | std::ios::in | std::ios::out | std::ios::trunc);
	
	if(!file_dCdvar.is_open()) cout << "\n\nError opening the file!\n\n";
	
	else
	{
		int VEC_CLASS_ID_BE = __builtin_bswap32(VEC_CLASS_ID);
		int nRows_BE        = __builtin_bswap32(nRows);
		
		file_dCdvar.write(reinterpret_cast<char*>(&VEC_CLASS_ID_BE), sizeof(int));
		file_dCdvar.write(reinterpret_cast<char*>(&nRows_BE), sizeof(int));
		
		for(int i = 0; i < nRows; i++)
		{
			double vecVal = dCdvar(i);
			char* vecValChar = reinterpret_cast<char*>(&vecVal);
			double vecVal_BE = reverseBytes(vecValChar);
			file_dCdvar.write(reinterpret_cast<char*>(&vecVal_BE), sizeof(double));
		}
		
		file_dCdvar.close();
	}
	
	system("echo Using PETSc...");
	system(petscArg);
	system("echo Completed!");
	Info << endl;
	
	// Reading adjoint solution for constraint
	std::fstream file_adjConstSolve("adjSolve.dat", std::ios::binary | std::ios::in);
	
	if(!file_adjConstSolve.is_open()) cout << "\n\nError opening the file!\n\n";
	
	else
	{
		int VEC_CLASS_ID_BE;
		int nRows_BE;
		
		file_adjConstSolve.read(reinterpret_cast<char*>(&VEC_CLASS_ID_BE), sizeof(int));
		file_adjConstSolve.read(reinterpret_cast<char*>(&nRows_BE), sizeof(int));
		
		VEC_CLASS_ID = __builtin_bswap32 (VEC_CLASS_ID_BE);
		nRows        = __builtin_bswap32 (nRows_BE);
		
		for(int i = 0; i < nRows; i++)
		{
			char val_BE[64];
			file_adjConstSolve.read(reinterpret_cast<char*>(&val_BE), sizeof(double));
			double val = reverseBytes(val_BE);
			adjConstSolve[i] = val;
		}
		
		file_adjConstSolve.close();
	}
	
	oldAdjConst = adjConstSolve;
}
