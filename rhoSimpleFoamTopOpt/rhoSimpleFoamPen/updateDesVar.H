scalar volTarget     (readScalar(optimizationProperties.lookup("volTarget"))); // Volume fraction (0 to 1)
scalar wVol          (readScalar(optimizationProperties.lookup("wVol")));      // Initial weight for volume constraint
scalar wVolMax       (readScalar(optimizationProperties.lookup("wVolMax")));   // Max weight for volume constraint
scalar gammaVol      (readScalar(optimizationProperties.lookup("gammaVol")));  // Multiplier for wVol
scalar stepSD        (readScalar(optimizationProperties.lookup("stepSD")));    // Step for Steepest Descent
int    itOptAMR      (readScalar(optimizationProperties.lookup("itOptAMR")));  // Apply AMR at each itOptAMR iterations (disabled if <= 0)
string parallelSolve (word(optimizationProperties.lookup("parallelSolve")));   // If parallel is enabled, look for cellSet in other dir
string objFunc (word(optimizationProperties.lookup("objFunc"))); // Select the objective function ("entropyVar" or "impellerEff" or "impellerTorque")

// Open nOptTemp.txt and read current optimization iteration (CAREFUL, NOT THE SAME CODE AS IN interpolationFunction.H)
std::ifstream tempoFile;
string optCont;
tempoFile.open("nOptTemp.txt");
getline(tempoFile, optCont, '\n');
tempoFile.close();
int optContVal = std::stoi(optCont);

// Setting the correct dir to search the cellSet
string searchDir = "constant/polyMesh/sets";
if((itOptAMR > 0) and (optContVal >= itOptAMR))
{
	if(parallelSolve=="yes") searchDir = runTime.timeName() + "/polyMesh/sets";
	else searchDir = runTime.timeName();
}

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

const labelList& desDomCells = designDomain.toc();

// Integral of rhoDes over the design domain
scalar integralDomain = 0.0;
forAll (desDomCells, ID) integralDomain += rhoDes[desDomCells[ID]]*mesh.V()[desDomCells[ID]];
reduce (integralDomain, sumOp<scalar>());

// Volume integral over design domain
scalar domainVolume = 0.0;
forAll(desDomCells, ID){domainVolume += mesh.V()[desDomCells[ID]];}
reduce (domainVolume, sumOp<scalar>());

scalar lambdaVol  = 0.0;
scalar volFrac    = integralDomain/domainVolume; // Volume fraction
scalar cVol       = volFrac - volTarget;         // Volume constraint

// Read from volumeConstraint.txt
if(runTime.value() > 0)
{
	std::ifstream volConstFile;
	string It;
	string wVolVal;
	string lambdaVolVal;
	string volFracVal;

	// Read last iteration from volumeConstraint.txt
	volConstFile.open("volumeConstraint.txt");
	while(!volConstFile.eof())
	{
		getline(volConstFile, It, ' ');
		getline(volConstFile, wVolVal, ' ');
		getline(volConstFile, lambdaVolVal, ' ');
		getline(volConstFile, volFracVal, '\n');
	}
	volConstFile.close();
	
	wVol      = stod(wVolVal);
	lambdaVol = stod(lambdaVolVal);

	wVol      = min(gammaVol*wVol, wVolMax);
	lambdaVol = lambdaVol - 2.0*wVol*cVol;

	// Volume constraint sensitivity
	scalar volConstMult = lambdaVol - 2.0*wVol*cVol;
	
	// Update design variables (inside design domain) - divided by cell volume as sensitivity because of RRT
	forAll(desDomCells, ID)
		rhoDes[desDomCells[ID]] = min(max(rhoDes[desDomCells[ID]] - stepSD*(Sensitivities[desDomCells[ID]] - (1.0/domainVolume)*volConstMult), 0.0), 1.0);
}
// File writing //
std::ofstream vcFile;
vcFile.open("volumeConstraint.txt", std::ios_base::app);
vcFile << "\n" << runTime.value() << " " << wVol << " " << lambdaVol << " " << volFrac;
vcFile.close();

#include "applyConst.H"
