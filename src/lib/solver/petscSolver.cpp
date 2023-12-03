#include "PetscSolver.h"
#include "petscmat.h"
#include "petscksp.h"
#include <iostream>
using namespace std;



PetscSolver::PetscSolver()
{
}


PetscSolver::~PetscSolver()
{
  //PetscPrintf(MPI_COMM_WORLD, "SolverPetsc::~SolverPetsc() \n");
  //cout << "SolverPetsc::~SolverPetsc() " << endl;

  //PetscPrintf(MPI_COMM_WORLD, "SolverPetsc::~SolverPetsc() \n");
  //cout << "SolverPetsc::~SolverPetsc() " << endl;
}

// Initialise Petsc solver
// size_local  - number of local rows/columns in the matrix
// size_global - number of global rows/columns in the matrix
// diag_nnz - number of nonzeros in the diagonal matrix
// offdiag_nnz - number of nonzeros in the off-diagonal matrix
// 
int PetscSolver::initialise(int size_local, int size_global, int* diag_nnz, int* offdiag_nnz)
{
    nRow = nCol = size_global;

    int dummy = 50;

    // Create PETSc vector
    errpetsc = VecCreate(PETSC_COMM_WORLD, &solnVec);
    CHKERRQ(errpetsc);

    errpetsc = VecSetSizes(solnVec, size_local, size_global);
    CHKERRQ(errpetsc);

    errpetsc = VecSetFromOptions(solnVec);
    CHKERRQ(errpetsc);

    errpetsc = VecDuplicate(solnVec, &rhsVec);
    CHKERRQ(errpetsc);

    errpetsc = VecSetOption(rhsVec, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);

    errpetsc = VecDuplicate(solnVec, &reacVec);
    CHKERRQ(errpetsc);


    // Create PETSc matrix
    errpetsc = MatCreate(PETSC_COMM_WORLD, &mtx);
    CHKERRQ(errpetsc);

    errpetsc = MatSetSizes(mtx, size_local, size_local, size_global, size_global);
    CHKERRQ(errpetsc);

    errpetsc = MatSetFromOptions(mtx);
    CHKERRQ(errpetsc);

    errpetsc = MatMPIAIJSetPreallocation(mtx, dummy, diag_nnz, dummy, offdiag_nnz);
    CHKERRQ(errpetsc);

    errpetsc = MatSeqAIJSetPreallocation(mtx, dummy, diag_nnz);
    CHKERRQ(errpetsc);

    errpetsc = MatSetOption(mtx, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    CHKERRQ(errpetsc);

    errpetsc = MatSetOption(mtx, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
    CHKERRQ(errpetsc);

    errpetsc = MatSetOption(mtx, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    CHKERRQ(errpetsc);

    // Create the KSP context
    errpetsc = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRQ(errpetsc);

    // Set the operators for the KSP context
    errpetsc = KSPSetOperators(ksp, mtx, mtx);
    CHKERRQ(errpetsc);

    //  Set whether to use non-zero initial guess or not
    //  KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    //  KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);


    // Set KSP options from the input file
    // This is convenient as it allows to choose different options
    // from the input files instead of recompiling the code
    errpetsc = KSPSetFromOptions(ksp);    CHKERRQ(errpetsc);

    errpetsc = KSPGetPC(ksp, &pc);    CHKERRQ(errpetsc);

    errpetsc = PCSetFromOptions(pc);    CHKERRQ(errpetsc);

    currentStatus = SOLVER_EMPTY;

    return 0;
}



int PetscSolver::factorise()
{
  if (currentStatus != ASSEMBLY_OK)
  {
    cerr << " SolverPetsc::factorise ... assemble matrix first!" << endl;
    return -1;
  }


  currentStatus = FACTORISE_OK;

  return 0;
}


int PetscSolver::solve()
{  
  //PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ...  \n\n");
  if (currentStatus != FACTORISE_OK)
  {
    cerr << " SolverPetsc::solve ... factorise matrix first!" << endl;
    return -1;
  }

  errpetsc = MatAssemblyBegin(mtx,MAT_FINAL_ASSEMBLY); CHKERRQ(errpetsc);
  errpetsc = MatAssemblyEnd(mtx,MAT_FINAL_ASSEMBLY); CHKERRQ(errpetsc);

  //errpetsc = MatView(mtx,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(errpetsc);


  //PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... Matrix Assembly ...  \n\n");

  errpetsc = VecAssemblyBegin(rhsVec); CHKERRQ(errpetsc);
  errpetsc = VecAssemblyEnd(rhsVec); CHKERRQ(errpetsc);
  
  //errpetsc = VecView(rhsVec,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(errpetsc);
  //PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... rhsVec Assembly ...  \n\n");

  errpetsc = VecAssemblyBegin(solnVec); CHKERRQ(errpetsc);
  errpetsc = VecAssemblyEnd(solnVec); CHKERRQ(errpetsc);

  //PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... soln Assembly ...  \n\n");

  errpetsc = VecZeroEntries(solnVec); CHKERRQ(errpetsc);

  PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... vec zero ...  \n\n");

  KSPSolve(ksp, rhsVec, solnVec); CHKERRQ(errpetsc);

  PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... KSP solve ...  \n\n");

  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);

  PetscInt its;

  errpetsc = KSPGetIterationNumber(ksp, &its); CHKERRQ(errpetsc);

  if(reason<0)
  {
    PetscPrintf(MPI_COMM_WORLD, "\n Divergence... %d iterations. \n\n", its);
    cout <<  reason << endl;
    exit(1);
    return -1;
  }
  else
  {
    PetscPrintf(MPI_COMM_WORLD, "\n Convergence in %d iterations.\n\n", its);
  }

  return 0;
}


int PetscSolver::factoriseAndSolve()
{
  char fct[] = "SolverPetsc::factoriseAndSolve";

  if(currentStatus != ASSEMBLY_OK)
  {
    cerr << " SolverPetsc::factoriseAndSolve ... assemble matrix first!" << endl;
    return -1;
  }

  factorise();
  solve();

  return 0;
}


int PetscSolver::assembleMatrixAndVectorSerial(vector<int>& forAssyElem, vector<int>& forAssyElemRHS, MatrixXd& Klocal, VectorXd& Flocal)
{
  int  size1 = forAssyElem.size();
  int  size2 = forAssyElemRHS.size();
  MatrixXdRM Klocal2 = Klocal;

  VecSetValues(rhsVec, size1, &forAssyElemRHS[0], &Flocal[0], INSERT_VALUES);
  MatSetValues(mtx,    size1, &forAssyElem[0], size1, &forAssyElemRHS[0], &Klocal2(0,0), ADD_VALUES);

  return 0;
}


int PetscSolver::zeroMtx()
{

  errpetsc = MatAssemblyBegin(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(errpetsc);
  errpetsc = MatAssemblyEnd(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(errpetsc);

  VecAssemblyBegin(rhsVec);
  VecAssemblyEnd(rhsVec);
  //PetscPrintf(MPI_COMM_WORLD, " 2 \n\n\n");
  VecZeroEntries(rhsVec);

  VecAssemblyBegin(reacVec);
  VecAssemblyEnd(reacVec);

  VecZeroEntries(reacVec);

  return 0;
}

int PetscSolver::free()
{
  //if(FREED) return 0;

  //PetscPrintf(MPI_COMM_WORLD, "SolverPetsc::free() \n");

  errpetsc = VecDestroy(&solnVec);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;
  //errpetsc = VecDestroy(&solnVecPrev);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = VecDestroy(&rhsVec);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = VecDestroy(&reacVec);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = MatDestroy(&mtx);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = PCReset(pc);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = KSPDestroy(&ksp);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = KSPReset(ksp);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;

  //PetscPrintf(MPI_COMM_WORLD, "SolverPetsc::free() \n");

  return 0;
}
