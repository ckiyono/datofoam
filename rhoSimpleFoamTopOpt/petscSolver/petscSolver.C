static char help[] = "Adjoint solver\n";

#include <iostream>
#include <fstream>
#include <petscksp.h>

using namespace std;

int main(int argc,char **args)
{	
	KSP            ksp;  // linear solver context
	Vec            x, b; // old solution, solution, RHS
	Mat            A;    // linear system matrix
	PetscErrorCode ierr = 0;
	PetscViewer    AViewer, bViewer, oldXViewer, xViewer;
	PetscReal      norm;
	PetscMPIInt    size;
	
	ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRMPI(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"MPI size %D\n",size);CHKERRQ(ierr);
	
	// Reading dRdvar.dat
	ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
	ierr = MatSetType(A,MATAIJ);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"dRdvar.dat",FILE_MODE_READ,&AViewer);CHKERRQ(ierr);
	ierr = MatLoad(A,AViewer);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&AViewer);CHKERRQ(ierr);
	
	// Reading dJdvar.dat
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"dJdvar.dat",FILE_MODE_READ,&bViewer);CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
	ierr = VecLoad(b,bViewer);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&bViewer);CHKERRQ(ierr);
	ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
	
	// Reading oldAdj.dat
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"oldAdj.dat",FILE_MODE_READ,&oldXViewer);CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
	ierr = VecLoad(x,oldXViewer);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&oldXViewer);CHKERRQ(ierr);
	
	//Solving Ax = b via PETSc KSP method
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
	ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
	
	//Writing adjSolve.dat
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"adjSolve.dat",FILE_MODE_WRITE,&xViewer);CHKERRQ(ierr);
	ierr = VecView(x,xViewer);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&xViewer);CHKERRQ(ierr);
		
	//- Destroying vectors and matrices -//
	ierr = VecDestroy(&b);CHKERRQ(ierr);
	ierr = VecDestroy(&x);CHKERRQ(ierr);
	ierr = MatDestroy(&A);CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
	
	ierr = PetscFinalize();
  
	return ierr;
  
}
/* ------------------------------------------------------------------- */
