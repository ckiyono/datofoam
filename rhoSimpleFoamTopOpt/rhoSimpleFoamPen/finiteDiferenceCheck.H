word FDcheck = word(optimizationProperties.lookup("FDcheck"));              // Enable FD check
scalar FDcellID  (readScalar(optimizationProperties.lookup("FDcellID")));   // ID of cell  (only if FDcheck = yes)
scalar FDcellVal (readScalar(optimizationProperties.lookup("FDcellVal" ))); // Value to FD (only if FDcheck = yes)

if(FDcheck=="yes")
{
	Info << "\nFinite Differences verification enabled!" << endl;
	Info << "Cell: " << FDcellID << " Value: " << FDcellVal << "\n" << endl;

	rhoDes[FDcellID] = FDcellVal;
}

