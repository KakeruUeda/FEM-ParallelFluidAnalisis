#include "FEM.h"
using namespace std;

void FEM::Stokes(){

  PetscPrintf(MPI_COMM_WORLD, "\n\n\n Solve Stokes equation \n");

  int  aa, bb, ee, ii, jj, kk, count, row, col, jpn, n1, n2, size1, size2;

  double  norm_rhs, timer;

  VectorXd  reacVec(numOfNodeGlobalFluid*numOfDofsNode);
  jpn = numOfNodeInElm*numOfDofsNode;
  VectorXd  Flocal(jpn);
  MatrixXd  Klocal(jpn, jpn);
  
  PetscScalar *arrayTempSolnFluid;

  Vec            vec_SEQ;
  VecScatter     ctx;
  VecScatterCreateToAll(solverPetscFluid->solnVec, &ctx, &vec_SEQ);
  MPI_Barrier(MPI_COMM_WORLD);

  solverPetscFluid->zeroMtx();
  reacVec.setZero();
  
  assignBoundaryConditions();
  
  MPI_Barrier(MPI_COMM_WORLD);
  PetscPrintf(MPI_COMM_WORLD, "\n ****** ELEMENT LOOP START ******* \n", timer);
  
  timer = MPI_Wtime();
  
  for(int ic=0;ic<numOfElmGlobalFluid;ic++){
    if(elmFluid[ic]->getSubdomainId() == myId){
      
      Klocal.setZero();
      Flocal.setZero();

      if(phiEXFluid[ic]>0.999){
        calcStokesMatrix(ic,Klocal,Flocal);
      }else{
        calcStokesMatrixXFEM(ic,Klocal,Flocal);
      }
      
      size1 = elmFluid[ic]->nodeForAssyBCsFluid.size();
      applyBoundaryConditions(Flocal, size1, ic);
      solverPetscFluid->assembleMatrixAndVectorSerial(elmFluid[ic]->nodeForAssyBCsFluid, elmFluid[ic]->nodeForAssyFluid, Klocal, Flocal);
    }
  }
    
  MPI_Barrier(MPI_COMM_WORLD);

  timer = MPI_Wtime() - timer;
  //computerTimeAssembly += timerVal;
  PetscPrintf(MPI_COMM_WORLD, "\n Time for matrix assembly = %f seconds \n", timer);

  PetscPrintf(MPI_COMM_WORLD, "\n ****** ELEMENT LOOP END ****** \n\n", timer);

  VecAssemblyBegin(solverPetscFluid->rhsVec);
  VecAssemblyEnd(solverPetscFluid->rhsVec);

  //VecNorm(solverPetsc->rhsVec, NORM_2, &norm_rhs);
  solverPetscFluid->currentStatus = ASSEMBLY_OK;
  
  MPI_Barrier(MPI_COMM_WORLD);
  PetscPrintf(MPI_COMM_WORLD, "\n ****** SOLVING LINEAR EQUATION ****** \n", timer);
  timer = MPI_Wtime();

  solverPetscFluid->factoriseAndSolve();
  
  MPI_Barrier(MPI_COMM_WORLD);
  timer = MPI_Wtime() - timer;
  //computerTimeSolver += timerVal;
  PetscPrintf(MPI_COMM_WORLD, "\n Time for PETSc solver = %f seconds \n", timer);
  PetscPrintf(MPI_COMM_WORLD, "\n ****** SOLVING LINEAR EQUATION END ****** \n\n", timer);
  

  /////////////////////////////////////////////////////////////////////////////
  // get the solution vector onto all the processors
  /////////////////////////////////////////////////////////////////////////////
  
  VecScatterBegin(ctx, solverPetscFluid->solnVec, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(ctx, solverPetscFluid->solnVec, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);

  VecGetArray(vec_SEQ, &arrayTempSolnFluid);

  // update solution vector
  for(ii=0; ii<numOfDofsGlobalFluid; ii++){
    SolnDataFluid.soln[assyForSolnFluid[ii]]   +=  arrayTempSolnFluid[ii];
  }

  MPI_Barrier(MPI_COMM_WORLD);
  
  VecRestoreArray(vec_SEQ, &arrayTempSolnFluid);
  
  VecScatterDestroy(&ctx);
  VecDestroy(&vec_SEQ);
}


void FEM::assignBoundaryConditions(){
  
  int ii, n1, n2, n3;
  DirichletBCsFluid.resize(numOfDofsNode*numOfNodeGlobalFluid);
  for(ii=0; ii<numOfBdNodeFluid; ii++){
    n1 = int(DirichletBCsFluid_tmp[ii][0]);
    n2 = int(DirichletBCsFluid_tmp[ii][1]);
    n3 = n1*numOfDofsNode+n2;
    DirichletBCsFluid[n3] = DirichletBCsFluid_tmp[ii][2];
  }
  return;

}

void FEM::applyBoundaryConditions(VectorXd& Flocal, const int size, const int ic){

  int ii;
  for(ii=0; ii<size; ii++){
    if(elmFluid[ic]->nodeForAssyBCsFluid[ii] == -1) {
      Flocal(ii) = DirichletBCsFluid[elmFluid[ic]->globalDOFnumsFluid[ii]];    
    }
  }

}