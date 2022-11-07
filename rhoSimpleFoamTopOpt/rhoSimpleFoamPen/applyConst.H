string optConst (word(optimizationProperties.lookup("optConst"))); // Select the constraint function ("velocProfile" or "pRatio")

faceSet unifDomain
(
	IOobject
	(
		"unifDomain", 
		searchDir,
		mesh, 
		IOobject::READ_IF_PRESENT, 
		IOobject::NO_WRITE
	)
);

if(optConst == "velocProfile")
{
	vector uRef (vector(optimizationProperties.lookup("uRef"))); // Velocity reference profile
	scalar uProf = 0;                                            // Velocity profile constraint

	const labelList& unifDomainList = unifDomain.toc();

	forAll(unifDomainList, ID)
	{
		int faceID = unifDomainList[ID];
		vector uCalc = U[mesh.owner()[faceID]] - uRef;
		uProf = uProf + (uCalc & uCalc)*mag(mesh.Sf()[faceID]);
	}
	
	reduce (uProf, sumOp<scalar>());
	
	scalar wUProf     (readScalar(optimizationProperties.lookup("wUProf")));     // Initial weight
	scalar wUProfMax  (readScalar(optimizationProperties.lookup("wUProfMax")));  // Max weight
	scalar gammaUProf (readScalar(optimizationProperties.lookup("gammaUProf"))); // Multiplier
	scalar lambdaUProf  = 0.0;

	// Read from uProfConst.txt
	if(runTime.value() > 0)
	{
		std::ifstream uProfConstFile;
		string It;
		string wUProfVal;
		string lambdaUProfVal;
		string uProfVal;

		// Read last iteration from uProfConst.txt
		uProfConstFile.open("uProfConst.txt");
		while(!uProfConstFile.eof())
		{
			getline(uProfConstFile, It, ' ');
			getline(uProfConstFile, wUProfVal, ' ');
			getline(uProfConstFile, lambdaUProfVal, ' ');
			getline(uProfConstFile, uProfVal, '\n');
		}
		uProfConstFile.close();
		
		wUProf      = stod(wUProfVal);
		lambdaUProf = stod(lambdaUProfVal);

		wUProf      = min(gammaUProf*wUProf, wUProfMax);
		lambdaUProf = lambdaUProf - 2.0*wUProf*uProf;

		// Velocity profile constraint multiplier
		scalar uProfConstMult = lambdaUProf - 2.0*wUProf*uProf;
		
		volScalarField constSens
		(
			IOobject
			(
				"constSens",
				runTime.timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			),
			mesh
		);
		
		// Update design variables (inside design domain) - divided by cell volume as sensitivity because of RRT
		forAll(desDomCells, ID)
			rhoDes[desDomCells[ID]] = min(max(rhoDes[desDomCells[ID]] + stepSD*(uProfConstMult*constSens[desDomCells[ID]]/mesh.V()[desDomCells[ID]]), 0.0), 1.0);	   
	}
	// File writing //
	std::ofstream UPCFile;
	UPCFile.open("uProfConst.txt", std::ios_base::app);
	UPCFile << "\n" << runTime.value() << " " << wUProf << " " << lambdaUProf << " " << uProf;
	UPCFile.close();
}

faceSet pRatioOut
(
	IOobject
	(
		"pRatioOut", 
		searchDir,
		mesh, 
		IOobject::READ_IF_PRESENT, 
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
		IOobject::READ_IF_PRESENT, 
		IOobject::NO_WRITE
	)
);

if(optConst == "pRatio")
{
	scalar pRatioConst = 0;
	scalar pOutCalc    = 0;
	scalar pInCalc     = 0;
	scalar areaOut     = 0;
	scalar areaIn      = 0;
	scalar pRatio (readScalar(optimizationProperties.lookup("pRatio"))); // Desired pressure ratio

	forAll(mesh.boundary(), patchI)
	{
		forAll(mesh.boundary()[patchI], faceI)
		{
			if(mesh.boundary()[patchI].name() == "outlet")
			{	
				pOutCalc = pOutCalc +
				p.boundaryField()[patchI][faceI]*mag(mesh.Sf().boundaryField()[patchI][faceI]);
				areaOut = areaOut + mag(mesh.Sf().boundaryField()[patchI][faceI]);
			}
			
			if(mesh.boundary()[patchI].name() == "inlet")
			{
				pInCalc = pInCalc +
				p.boundaryField()[patchI][faceI]*mag(mesh.Sf().boundaryField()[patchI][faceI]);
				areaIn = areaIn + mag(mesh.Sf().boundaryField()[patchI][faceI]);
			}
		}
	}
	
	reduce (pOutCalc, sumOp<scalar>());
	reduce (pInCalc, sumOp<scalar>());
	reduce (areaOut, sumOp<scalar>());
	reduce (areaIn, sumOp<scalar>());
	
	pRatioConst = (pOutCalc/areaOut) / (pInCalc/areaIn) - pRatio;
	
	scalar wPRatio     (readScalar(optimizationProperties.lookup("wPRatio")));     // Initial weight
	scalar wPRatioMax  (readScalar(optimizationProperties.lookup("wPRatioMax")));  // Max weight
	scalar gammaPRatio (readScalar(optimizationProperties.lookup("gammaPRatio"))); // Multiplier
	scalar lambdaPRatio  = 0.0;

	// Read from pRatioConst.txt
	if(runTime.value() > 0)
	{
		std::ifstream pRatioConstFile;
		string It;
		string wPRatioVal;
		string lambdaPRatioVal;
		string pRatioVal;

		// Read last iteration from pRatioConst.txt
		pRatioConstFile.open("pRatioConst.txt");
		while(!pRatioConstFile.eof())
		{
			getline(pRatioConstFile, It, ' ');
			getline(pRatioConstFile, wPRatioVal, ' ');
			getline(pRatioConstFile, lambdaPRatioVal, ' ');
			getline(pRatioConstFile, pRatioVal, '\n');
		}
		pRatioConstFile.close();
		
		wPRatio      = stod(wPRatioVal);
		lambdaPRatio = stod(lambdaPRatioVal);

		wPRatio      = min(gammaPRatio*wPRatio, wPRatioMax);
		lambdaPRatio = lambdaPRatio - 2.0*wPRatio*pRatioConst;

		// Pressure ratio constraint multiplier
		scalar pRatioConstMult = lambdaPRatio - 2.0*wPRatio*pRatioConst;
		
		volScalarField constSens
		(
			IOobject
			(
				"constSens",
				runTime.timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			),
			mesh
		);
		
		// Update design variables (inside design domain) - divided by cell volume as sensitivity because of RRT
		forAll(desDomCells, ID)
			rhoDes[desDomCells[ID]] = min(max(rhoDes[desDomCells[ID]] + stepSD*(pRatioConstMult*constSens[desDomCells[ID]]/mesh.V()[desDomCells[ID]]), 0.0), 1.0);	   
	}
	// File writing //
	std::ofstream PRCFile;
	PRCFile.open("pRatioConst.txt", std::ios_base::app);
	PRCFile << "\n" << runTime.value() << " " << wPRatio << " " << lambdaPRatio << " " << pRatioConst;
	PRCFile.close();	
}

if(optConst == "statP")
{
	scalar statPConst  = 0;
	scalar pOutCalc    = 0;
	scalar pInCalc     = 0;
	scalar areaOut     = 0;
	scalar areaIn      = 0;
	scalar statP (readScalar(optimizationProperties.lookup("statP"))); // Desired static pressure difference

	forAll(mesh.boundary(), patchI)
	{
		forAll(mesh.boundary()[patchI], faceI)
		{
			if(mesh.boundary()[patchI].name() == "outlet")
			{	
				pOutCalc = pOutCalc +
				p.boundaryField()[patchI][faceI]*mag(mesh.Sf().boundaryField()[patchI][faceI]);
				areaOut = areaOut + mag(mesh.Sf().boundaryField()[patchI][faceI]);
			}
			
			if(mesh.boundary()[patchI].name() == "inlet")
			{
				pInCalc = pInCalc +
				p.boundaryField()[patchI][faceI]*mag(mesh.Sf().boundaryField()[patchI][faceI]);
				areaIn = areaIn + mag(mesh.Sf().boundaryField()[patchI][faceI]);
			}
		}
	}
	
	reduce (pOutCalc, sumOp<scalar>());
	reduce (pInCalc, sumOp<scalar>());
	reduce (areaOut, sumOp<scalar>());
	reduce (areaIn, sumOp<scalar>());
	
	statPConst = ((pOutCalc/areaOut) - (pInCalc/areaIn)) - statP;
	
	scalar wStatP     (readScalar(optimizationProperties.lookup("wStatP")));     // Initial weight
	scalar wStatPMax  (readScalar(optimizationProperties.lookup("wStatPMax")));  // Max weight
	scalar gammaStatP (readScalar(optimizationProperties.lookup("gammaStatP"))); // Multiplier
	scalar lambdaStatP  = 0.0;

	// Read from statPConst.txt
	if(runTime.value() > 0)
	{
		std::ifstream statPConstFile;
		string It;
		string wStatPVal;
		string lambdaStatPVal;
		string statPVal;

		// Read last iteration from statPConst.txt
		statPConstFile.open("statPConst.txt");
		while(!statPConstFile.eof())
		{
			getline(statPConstFile, It, ' ');
			getline(statPConstFile, wStatPVal, ' ');
			getline(statPConstFile, lambdaStatPVal, ' ');
			getline(statPConstFile, statPVal, '\n');
		}
		statPConstFile.close();
		
		wStatP      = stod(wStatPVal);
		lambdaStatP = stod(lambdaStatPVal);

		wStatP      = min(gammaStatP*wStatP, wStatPMax);
		lambdaStatP = lambdaStatP - 2.0*wStatP*statPConst;

		// Pressure ratio constraint multiplier
		scalar statPConstMult = lambdaStatP - 2.0*wStatP*statPConst;
		
		volScalarField constSens
		(
			IOobject
			(
				"constSens",
				runTime.timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			),
			mesh
		);
		
		// Update design variables (inside design domain) - divided by cell volume as sensitivity because of RRT
		forAll(desDomCells, ID)
			rhoDes[desDomCells[ID]] = min(max(rhoDes[desDomCells[ID]] + stepSD*(statPConstMult*constSens[desDomCells[ID]]/mesh.V()[desDomCells[ID]]), 0.0), 1.0);	   
	}
	// File writing //
	std::ofstream SPCFile;
	SPCFile.open("statPConst.txt", std::ios_base::app);
	SPCFile << "\n" << runTime.value() << " " << wStatP << " " << lambdaStatP << " " << statPConst;
	SPCFile.close();	
}
