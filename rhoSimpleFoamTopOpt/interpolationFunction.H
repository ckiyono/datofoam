// Interpolation function from Borrvall and Petersson (2003)
scalar qGray    (readScalar(optimizationProperties.lookup("qGray")));    // Grayness factor
scalar qGrayMax (readScalar(optimizationProperties.lookup("qGrayMax"))); // Max grayness factor
int    qIter    (readScalar(optimizationProperties.lookup("qIter")));    // Change at qIter iterations
int    qMult    (readScalar(optimizationProperties.lookup("qMult")));    // Multiply qGray by qMult

scalar alphaMin (readScalar(optimizationProperties.lookup("alphaMin")));
scalar alphaMax (readScalar(optimizationProperties.lookup("alphaMax")));

// Open nOptTemp.txt and read current optimization iteration
std::ifstream tempFile;
string nCurr;
tempFile.open("nOptTemp.txt");
getline(tempFile, nCurr, '\n');
tempFile.close();
int nCurrVal = std::stoi(nCurr);

// Change qGray
qGray = max(min((qGray*Foam::pow(qMult, (nCurrVal/qIter))), qGrayMax), qGray);

forAll(mesh.C(), cellID)
	alphaPen[cellID] = alphaMax + (alphaMin - alphaMax)*rhoDes[cellID]*(1.0 + qGray)/(rhoDes[cellID] + qGray);
