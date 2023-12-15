#include "FEM.h"

void FEM::UnsteadyNavierStokes(){

  PetscPrintf(MPI_COMM_WORLD, "\n\n\n UNSTEADY NAVIER STOKES EQUATION \n");

  int  aa, bb, ee, ii, jj, kk, count, row, col, jpn, n1, n2, size1, size2;
  double  norm_rhs, timer;

  double norm,norm0;
  double tmp = 1e12;

  uFluid.resize(numOfNodeGlobalFluid,0e0);
  vFluid.resize(numOfNodeGlobalFluid,0e0);
  wFluid.resize(numOfNodeGlobalFluid,0e0);
  pFluid.resize(numOfNodeGlobalFluid,0e0);

  uf.resize(timeMax,vector<double>(numOfNodeGlobalFluid,0e0));
  vf.resize(timeMax,vector<double>(numOfNodeGlobalFluid,0e0));
  wf.resize(timeMax,vector<double>(numOfNodeGlobalFluid,0e0));
  pf.resize(timeMax,vector<double>(numOfNodeGlobalFluid,0e0));

  u.resize(numOfNodeGlobal,0e0);
  v.resize(numOfNodeGlobal,0e0);
  w.resize(numOfNodeGlobal,0e0);
  p.resize(numOfNodeGlobal,0e0);

  VectorXd  reacVec(numOfNodeGlobalFluid*numOfDofsNode);
  jpn = numOfNodeInElm*numOfDofsNode;
  VectorXd  Flocal(jpn);
  MatrixXd  Klocal(jpn, jpn);
  
  PetscScalar *arrayTempSolnFluid;

  Vec            vec_SEQ;
  VecScatter     ctx;
  VecScatterCreateToAll(solverPetscFluid->solnVec, &ctx, &vec_SEQ);
  MPI_Barrier(MPI_COMM_WORLD);
  
  assignBCs();


  for(int t_itr=0; t_itr<timeMax; t_itr++){

    solverPetscFluid->zeroMtx();
    reacVec.setZero();
    
    applyBCs();

    MPI_Barrier(MPI_COMM_WORLD);
    timer = MPI_Wtime();
    
    for(int ic=0;ic<numOfElmGlobalFluid;ic++){
      if(elmFluid[ic]->getSubdomainId() == myId){

        Klocal.setZero();
        Flocal.setZero();

        if(phiEXFluid[ic]>0.999){
          MatAssyUSNS(Klocal,Flocal,ic,t_itr);
        }else{
          //XFEM_SteadyNavierStokesMatrixXFEM(ic,Klocal,Flocal);
        }
        solverPetscFluid->assembleMatrixAndVectorSerial(elmFluid[ic]->nodeForAssyBCsFluid, elmFluid[ic]->nodeForAssyFluid, Klocal, Flocal);
      }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    timer = MPI_Wtime() - timer;
    //computerTimeAssembly += timerVal;
    PetscPrintf(MPI_COMM_WORLD, "\n ****** Time for matrix assembly = %f seconds ****** \n", timer);

    VecAssemblyBegin(solverPetscFluid->rhsVec);
    VecAssemblyEnd(solverPetscFluid->rhsVec);

    //VecNorm(solverPetsc->rhsVec, NORM_2, &norm_rhs);
    solverPetscFluid->currentStatus = ASSEMBLY_OK;
  
    MPI_Barrier(MPI_COMM_WORLD);

    timer = MPI_Wtime();

    solverPetscFluid->factoriseAndSolve();
  
    MPI_Barrier(MPI_COMM_WORLD);
    timer = MPI_Wtime() - timer;
    //computerTimeSolver += timerVal;

    PetscPrintf(MPI_COMM_WORLD, "\n ****** Time for PETSc solver = %f seconds ****** \n", timer);
    
    for(ii=0; ii<numOfDofsGlobalFluid; ii++){
      SolnDataFluid.soln[assyForSolnFluid[ii]] = 0e0;
    }

    VecScatterBegin(ctx, solverPetscFluid->solnVec, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, solverPetscFluid->solnVec, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);

    VecGetArray(vec_SEQ, &arrayTempSolnFluid);

    // update solution vector
    for(ii=0; ii<numOfDofsGlobalFluid; ii++){
      SolnDataFluid.soln[assyForSolnFluid[ii]]   +=  arrayTempSolnFluid[ii];
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    VecRestoreArray(vec_SEQ, &arrayTempSolnFluid);

    postCaluculation_timeItr(t_itr);

    if(myId == 0){
      double time_now = t_itr*dt;
      printf(" Time = %e \n",time_now);
    }
  }

  VecScatterDestroy(&ctx);
  VecDestroy(&vec_SEQ);
}