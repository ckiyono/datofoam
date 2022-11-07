densVec adjSolve (Eigen::VectorXd::Zero((3 + 2)*nCells + nFaces));
densVec adjConstSolve (Eigen::VectorXd::Zero((3 + 2)*nCells + nFaces));

Info << "\nSolving adjoint problem...\n" << endl;

{
    string dafoamArg (word(optimizationProperties.lookup("dafoamArg"))); // Argument for dafoam
	std::replace(petscArg.begin(), petscArg.end(), '|', ' ');          // Replace all '|' to ' '
	std::replace(petscArg.begin(), petscArg.end(), '*', '|');          // Replace all '*' to '|'

    system("echo Using DAFoam...\n");
	system(dafoamArg);
	system("echo Completed!");

    #include "optConst.H"  // Calculate and write RHS for adjoint problem of constraint
}