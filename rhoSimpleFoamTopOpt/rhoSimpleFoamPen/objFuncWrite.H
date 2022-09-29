// File writing //
std::ofstream JFile;
JFile.open("objectiveFunction.txt", std::ios_base::app);
JFile << "\n" << runTime.value() << " " << J;
JFile.close();
