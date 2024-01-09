#include "FEM.h"
using namespace std;

void FEM::SteadyNavierStokes(){

  PetscPrintf(MPI_COMM_WORLD, "\n\n\n STEADY NAVIER STOKES EQUATION \n");

  int  aa, bb, ee, ii, jj, kk, count, row, col, jpn, n1, n2, size1, size2;
  double  norm_rhs, timer;

  double norm,norm0;
  double tmp = 1e12;

  jpn = numOfNodeInElm*numOfDofsNode;
  VectorXd  Flocal(jpn);
  MatrixXd  Klocal(jpn, jpn);
  
  PetscScalar *arrayTempSolnFluid;

  Vec            vec_SEQ;
  VecScatter     ctx;
  VecScatterCreateToAll(solverPetscFluid->solnVec, &ctx, &vec_SEQ);
  MPI_Barrier(MPI_COMM_WORLD);
  
  assignBCs();
  solverPetscFluid->initialAssembly();

  for(ii=0; ii<numOfNodeGlobalFluid; ii++){
    if(bd_iu_fluid[ii][0] == 0) uFluid[ii] = bd_u_fluid[ii][0];
    if(bd_iu_fluid[ii][1] == 0) vFluid[ii] = bd_u_fluid[ii][1];
    if(bd_iu_fluid[ii][2] == 0) wFluid[ii] = bd_u_fluid[ii][2];
    if(bd_ip_fluid[ii] == 0) pFluid[ii] = bd_p_fluid[ii];
  }
  for(ii=0; ii<numOfDofsNode*numOfNodeGlobalFluid; ii++){
    DirichletBCsFluid[ii] = 0e0;
  }

  for(int loop=0;loop<NRitr;loop++)
  {  
    solverPetscFluid->setZeroValue();
    
    applyBCs();

    MPI_Barrier(MPI_COMM_WORLD);
    timer = MPI_Wtime();
    
    for(int ic=0;ic<numOfElmGlobalFluid;ic++)
    {
      if(elmFluid[ic]->getSubdomainId() == myId)
      {
        Klocal.setZero();
        Flocal.setZero();
        
        switch(bd)
        {
          case BOUNDARY::XFEM:
            if(phiEXFluid[ic]>0.999){
              MatAssySNS(ic,Klocal,Flocal);
            }else{
              PetscPrintf(MPI_COMM_WORLD, "XFEM not yet implemented");
              exit(1);
            }
            break;

          case BOUNDARY::DARCY:
            if(phiVOFFluid[ic]>0.999){
              MatAssySNS(ic,Klocal,Flocal);
            }else{
              PetscPrintf(MPI_COMM_WORLD, "Darcy not yet implemented");
              exit(1);
            }
            break;

          default:
            PetscPrintf(MPI_COMM_WORLD, "undefined boundary method");
            exit(1);
        }
        solverPetscFluid->setValue(elmFluid[ic]->nodeForAssyBCsFluid, elmFluid[ic]->nodeForAssyFluid, Klocal, Flocal);
      }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    timer = MPI_Wtime() - timer;
    //computerTimeAssembly += timerVal;
    PetscPrintf(MPI_COMM_WORLD, "\n ****** Time for matrix assembly = %f seconds ****** \n", timer);

    //VecAssemblyBegin(solverPetscFluid->rhsVec);
    //VecAssemblyEnd(solverPetscFluid->rhsVec);

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

    norm = solverPetscFluid->vectorNorm(numOfDofsGlobalFluid,SolnDataFluid.soln);
    if(loop==0) norm0 = norm;
    
    postCaluculation_itr(loop);

    if(myId == 0){
      printf(" NR iter.=%d norm/norm0=%e\n",loop,norm/norm0);
    }

    if(norm/norm0<NRtolerance) break;
    tmp = norm/norm0;
  }

  VecScatterDestroy(&ctx);
  VecDestroy(&vec_SEQ);
}

