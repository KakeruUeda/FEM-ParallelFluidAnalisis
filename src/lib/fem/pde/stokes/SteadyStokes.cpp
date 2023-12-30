#include "FEM.h"

void FEM::Stokes(){

  PetscPrintf(MPI_COMM_WORLD, "\n\n\n STEADY STOKES EQUATION \n");

  int  aa, bb, ee, ii, jj, kk, count, row, col, jpn, n1, n2, size1, size2;

  double  norm_rhs, timer;

  uFluid.resize(numOfNodeGlobalFluid,0); u.resize(numOfNodeGlobal,0);
  vFluid.resize(numOfNodeGlobalFluid,0); v.resize(numOfNodeGlobal,0);
  wFluid.resize(numOfNodeGlobalFluid,0); w.resize(numOfNodeGlobal,0);
  pFluid.resize(numOfNodeGlobalFluid,0); p.resize(numOfNodeGlobal,0);


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

  assignBCs();
  applyBCs();

  MPI_Barrier(MPI_COMM_WORLD);
  PetscPrintf(MPI_COMM_WORLD, "\n ****** ELEMENT LOOP START ******* \n", timer);
  
  timer = MPI_Wtime();
  //#pragma omp parallel for
  for(int ic=0;ic<numOfElmGlobalFluid;ic++){
    if(elmFluid[ic]->getSubdomainId() == myId){
      
      Klocal.setZero();
      Flocal.setZero();

      switch(bd){
        case BOUNDARY::XFEM:
          if(phiEXFluid[ic]>0.999){
            MatAssySTT(ic,Klocal,Flocal);
          }else{
            XFEM_MatAssySTT(ic,Klocal,Flocal);
          }
          break;

        case BOUNDARY::DARCY:
          if(phiFluid[ic]>0.999){
            MatAssySTT(ic,Klocal,Flocal);
          }else{
            Darcy_MatAssySTT(ic,Klocal,Flocal);
          }
          break;

        default:
          cout << "undefined boundary method" << endl;
          exit(1);
      }

      solverPetscFluid->setValue(elmFluid[ic]->nodeForAssyBCsFluid, elmFluid[ic]->nodeForAssyFluid, Klocal, Flocal);
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
  

  //get the solution vector 
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



