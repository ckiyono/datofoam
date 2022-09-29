int itOptAMR (readScalar(optimizationProperties.lookup("itOptAMR"))); // Apply AMR at each itOptAMR iterations (disabled if <= 0)
string oneOrAll (word(optimizationProperties.lookup("oneOrAll"))); // If it is "one", only 1 AMR. If it is "all" use in all its
if((oneOrAll == "one") and (optCont > itOptAMR)) itOptAMR = 0;

if((itOptAMR > 0) and (optCont > 0) and (optCont >= itOptAMR) and (optCont % itOptAMR == 0))
{
	string optConst (word(optimizationProperties.lookup("optConst"))); // Select the constraint function ("velocProfile" or "pRatio")
	
	volScalarField constSens
	(
		IOobject
		(
			"constSens",
			runTime.timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("zero",dimensionSet(0,0,0,0,0,0,0),0)
	);
	
	volScalarField magGradRhoDes = mag(fvc::grad(rhoDes));
	system("echo ");
	system("echo Adaptive mesh refinement enabled! Performing mesh operations...");
	system("echo ");
	system("rm -f orthoCols.txt");
	orthoCols.clear();
	colRowNZ.clear();
	
	while(simple.loop())
	{
		mesh.update(); // Apply AMR
		runTime.writeNow();
		psiC.write();
		psiU.write();
		psiE.write();
		Sensitivities.write();
		rhoDes.write();
		alphaPen.write();
		if((optConst == "velocProfile") or (optConst == "pRatio")) constSens.write();
		break;
	}
	
	// General mesh parameters update
	nCells = mesh.nCells(); // Total no of cells
	nPhi   = phi.size();    // No of internal face
	nFaces = nPhi;          // Total no of faces
	forAll(mesh.boundary(), patchI) nFaces = nFaces + mesh.boundary()[patchI].size(); // Counting with boundary faces
	
	system("topoSet");       // Create cellSet again
	#include "orthoGroups.H" // Create groups for new mesh
	system("echo Remesh and regroup completed!");
	system("echo ");
	
	searchDir = runTime.timeName() + "/polyMesh/sets";	
	
	cellSet designDomain
	(
		IOobject
		(
			"designDomain", 
			searchDir,
			mesh, 
			IOobject::MUST_READ, 
			IOobject::NO_WRITE
		)
	);
	designDomain.write();
	
	if(optConst == "velocProfile")
	{
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
		unifDomain.write();
	}
	
	if((optConst == "pRatio") or (objFunc == "impellerEff"))
	{
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
		pRatioOut.write();
		
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
		pRatioIn.write();
	}
}
