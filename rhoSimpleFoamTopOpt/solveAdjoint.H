densVec adjSolve (Eigen::VectorXd::Zero((3 + 2)*nCells + nFaces));
densVec adjConstSolve (Eigen::VectorXd::Zero((3 + 2)*nCells + nFaces));

Info << "\nSolving adjoint problem...\n" << endl;

{
	adjSolve = oldAdj;
	
	string petscArg (word(optimizationProperties.lookup("petscArg"))); // Argument for petsc solver
	std::replace(petscArg.begin(), petscArg.end(), '|', ' ');          // Replace all '|' to ' '
	std::replace(petscArg.begin(), petscArg.end(), '*', '/');          // Replace all '*' to '/'

	//petscArg = "$PETSC_DIR/$PETSC_ARCH/bin/" + petscArg;
	
	/*   PETSc file structures for matrix and vector
	PetscInt    MAT_FILE_CLASSID
	PetscInt    number of rows
	PetscInt    number of columns
	PetscInt    total number of nonzeros
	PetscInt    *number nonzeros in each row
	PetscInt    *column indices of all nonzeros (starting index is zero)
	PetscScalar *values of all nonzeros
	
	PetscInt    VEC_FILE_CLASSID
	PetscInt    number of rows
	PetscScalar *values of all entries
	*/
	
	Info << "Writing data for adjoint system..." << flush;
	
	int MAT_CLASS_ID = 1211216; // ID for petsc to recognize that it is a matrix
	int VEC_CLASS_ID = 1211214; // ID for petsc to recognize that it is a vector
	int nRows = dRdvar.rows();
	int nCols = dRdvar.cols();
	int nnz   = dRdvar.nonZeros();
	
	// Writing matrix
	std::fstream file_dRdvar("dRdvar.dat", std::ios::binary | std::ios::in | std::ios::out | std::ios::trunc);
	
	if(!file_dRdvar.is_open()) cout << "\n\nError opening the file!\n\n";
	
	else
	{
		// PETSc binary works with big endian...thus, the _BE indicates big endian storage
		int MAT_CLASS_ID_BE = __builtin_bswap32(MAT_CLASS_ID);
		int nRows_BE        = __builtin_bswap32(nRows);
		int nCols_BE        = __builtin_bswap32(nCols);
		int nnz_BE          = __builtin_bswap32(nnz);
		
		file_dRdvar.write(reinterpret_cast<char*>(&MAT_CLASS_ID_BE), sizeof(int));
		file_dRdvar.write(reinterpret_cast<char*>(&nRows_BE), sizeof(int));
		file_dRdvar.write(reinterpret_cast<char*>(&nCols_BE), sizeof(int));
		file_dRdvar.write(reinterpret_cast<char*>(&nnz_BE), sizeof(int));
		
		for(int i = 0; i < nRows; i++)
		{
			int nnzPerRow         = dRdvar.row(i).nonZeros();
			int nnzPerRow_BE = __builtin_bswap32(nnzPerRow);
			file_dRdvar.write(reinterpret_cast<char*>(&nnzPerRow_BE), sizeof(int));
		}
		
		for(int k = 0; k < dRdvar.outerSize(); k++)
		{
			for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(dRdvar,k); it; ++it)
			{
				int sysCol    = it.col();
				int sysCol_BE = __builtin_bswap32(sysCol);
				file_dRdvar.write(reinterpret_cast<char*>(&sysCol_BE), sizeof(int));
			}
		}
		
		for(int k = 0; k < dRdvar.outerSize(); k++)
		{
			for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(dRdvar,k); it; ++it)
			{
				double sysVal = it.value();
				char* sysValChar = reinterpret_cast<char*>(&sysVal);
				double sysVal_BE = reverseBytes(sysValChar);
				file_dRdvar.write(reinterpret_cast<char*>(&sysVal_BE), sizeof(double));
			}
		}
		
		file_dRdvar.close();
	}
	
	// Writing source vector
	std::fstream file_dJdvar("dJdvar.dat", std::ios::binary | std::ios::in | std::ios::out | std::ios::trunc);
	
	if(!file_dJdvar.is_open()) cout << "\n\nError opening the file!\n\n";
	
	else
	{
		int VEC_CLASS_ID_BE = __builtin_bswap32(VEC_CLASS_ID);
		int nRows_BE        = __builtin_bswap32(nRows);
		
		file_dJdvar.write(reinterpret_cast<char*>(&VEC_CLASS_ID_BE), sizeof(int));
		file_dJdvar.write(reinterpret_cast<char*>(&nRows_BE), sizeof(int));
		
		for(int i = 0; i < nRows; i++)
		{
			double vecVal = dJdvar(i);
			char* vecValChar = reinterpret_cast<char*>(&vecVal);
			double vecVal_BE = reverseBytes(vecValChar);
			file_dJdvar.write(reinterpret_cast<char*>(&vecVal_BE), sizeof(double));
		}
		
		file_dJdvar.close();
	}
	
	// Writing initial guess
	std::fstream file_oldAdj("oldAdj.dat", std::ios::binary | std::ios::in | std::ios::out | std::ios::trunc);
	
	if(!file_oldAdj.is_open()) cout << "\n\nError opening the file!\n\n";
	
	else
	{
		int VEC_CLASS_ID_BE = __builtin_bswap32(VEC_CLASS_ID);
		int nRows_BE        = __builtin_bswap32(nRows);
		
		file_oldAdj.write(reinterpret_cast<char*>(&VEC_CLASS_ID_BE), sizeof(int));
		file_oldAdj.write(reinterpret_cast<char*>(&nRows_BE), sizeof(int));
		
		for(int i = 0; i < nRows; i++)
		{
			double vecVal = oldAdj(i);
			char* vecValChar = reinterpret_cast<char*>(&vecVal);
			double vecVal_BE = reverseBytes(vecValChar);
			file_oldAdj.write(reinterpret_cast<char*>(&vecVal_BE), sizeof(double));
		}
		
		file_oldAdj.close();
	}
	
	system("echo Using PETSc...\n");
	system(petscArg);
	system("echo Completed!");
	Info << endl;
	
	// Reading adjoint solution
	std::fstream file_adjSolve("adjSolve.dat", std::ios::binary | std::ios::in);
	
	if(!file_adjSolve.is_open()) cout << "\n\nError opening the file!\n\n";
	
	else
	{
		int VEC_CLASS_ID_BE;
		int nRows_BE;
		
		file_adjSolve.read(reinterpret_cast<char*>(&VEC_CLASS_ID_BE), sizeof(int));
		file_adjSolve.read(reinterpret_cast<char*>(&nRows_BE), sizeof(int));
		
		VEC_CLASS_ID = __builtin_bswap32 (VEC_CLASS_ID_BE);
		nRows        = __builtin_bswap32 (nRows_BE);
		
		for(int i = 0; i < nRows; i++)
		{
			char val_BE[64];
			file_adjSolve.read(reinterpret_cast<char*>(&val_BE), sizeof(double));
			double val = reverseBytes(val_BE);
			adjSolve[i] = val;
		}
		
		file_adjSolve.close();
	}
	
	#include "optConst.H"  // Calculate and write RHS for adjoint problem of constraint
}

oldAdj = adjSolve;
