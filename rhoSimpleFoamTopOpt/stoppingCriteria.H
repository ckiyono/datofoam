{
	// Setting the correct dir to search the cellSet
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

	const labelList& desDomCells = designDomain.toc();
	
	scalar volTarget (readScalar(optimizationProperties.lookup("volTarget"))); // Volume fraction (0 to 1)
	scalar wVol      (readScalar(optimizationProperties.lookup("wVol")));      // Initial weight for volume constraint

	// Integral of rhoDes over the design domain
	scalar integralDomain = 0.0;
	forAll (desDomCells, ID) integralDomain += rhoDes[desDomCells[ID]]*mesh.V()[desDomCells[ID]];
	reduce (integralDomain, sumOp<scalar>());

	// Volume integral over design domain
	scalar domainVolume = 0.0;
	forAll(desDomCells, ID){domainVolume += mesh.V()[desDomCells[ID]];}
	reduce (domainVolume, sumOp<scalar>());

	scalar volFrac = integralDomain/domainVolume; // Volume fraction
	
	#include "objectiveFunction.H"

	JOld  = JAct;
	JAct  = J;
	JDiff = mag((JAct - JOld)/(JOld + 1.0e-9));

	Info << "Objective function difference = " << JDiff*100.0 << " %" << endl;

	// With volume constraint
	if((wVol > 1.0e-9) and (mag(JDiff) < optTol) and (mag(volFrac - volTarget) < optTol)) break;

	// Without volume constraint
	if((wVol < 1.0e-9) and (mag(JDiff) < optTol)) break;
}
