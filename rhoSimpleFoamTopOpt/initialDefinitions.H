// Optimization parameters
int nOpt (readScalar(optimizationProperties.lookup("nOpt")));          // Number of optimization iterations
double spsTol  (readScalar(optimizationProperties.lookup("spsTol")));  // Tolerance to insert in sparse matrix and vector
scalar stepJac (readScalar(optimizationProperties.lookup("stepJac"))); // Step for Jacobian calculation

// Definitions for stopping criteria
scalar optTol (readScalar(optimizationProperties.lookup("optTol")));
scalar JOld  = 0;
scalar JAct  = 0;
scalar JDiff = 0;

// General mesh parameters
int nCells = mesh.nCells(); // Total no of cells
int nPhi   = phi.size();    // No of internal face
int nFaces = nPhi;          // Total no of faces
forAll(mesh.boundary(), patchI) nFaces = nFaces + mesh.boundary()[patchI].size(); // Counting with boundary faces

// Initial guess vector
densVec oldAdj (Eigen::VectorXd::Zero((3 + 2)*nCells + nFaces));
densVec oldAdjConst (Eigen::VectorXd::Zero((3 + 2)*nCells + nFaces));
std::vector<std::vector<int>> orthoCols;  // Groups of orthogonal variables
std::vector<std::vector<int>> colRowNZ;   // Rows with non-zeros for each column
string searchDir = "constant/polyMesh/sets";

// Objective functions
string objFunc (word(optimizationProperties.lookup("objFunc"))); // Select the objective function ("entropyVar" or "impellerEff" or "impellerTorque")

// Constraint
string optConst (word(optimizationProperties.lookup("optConst"))); // Select the constraint function ("velocProfile" or "pRatio")

// Parallel run
string parallelSolve (word(optimizationProperties.lookup("parallelSolve")));

// Scaling for pertubation
scalar uScale   (readScalar(optimizationProperties.lookup("uScale")));   // Scaling for velocity perturbation
scalar pScale   (readScalar(optimizationProperties.lookup("pScale")));   // Scaling for pressure perturbation
scalar phiScale (readScalar(optimizationProperties.lookup("phiScale"))); // Scaling for flux perturbation
scalar hScale   (readScalar(optimizationProperties.lookup("hScale")));   // Scaling for enthalpy perturbation
