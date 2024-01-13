#include "FEM.h"

void FEM::UnsteadyNavierStokes()
{
  PetscPrintf(MPI_COMM_WORLD, "\n\n\n UNSTEADY NAVIER STOKES EQUATION \n");

  int  aa, bb, ee, ii, jj, kk, count, row, col, jpn, n1, n2, size1, size2;
  double  norm_rhs, timer;

  double norm,norm0;
  double tmp = 1e12;

  jpn = numOfNodeInElm * numOfDofsNode;
  VectorXd  Flocal(jpn);
  MatrixXd  Klocal(jpn, jpn);
  
  PetscScalar *arrayTempSolnFluid;

  Vec            vec_SEQ;
  VecScatter     ctx;

  VecScatterCreateToAll(solverPetscFluid->solnVec, &ctx, &vec_SEQ);
  
  MPI_Barrier(MPI_COMM_WORLD);
  assignBCs();
  
  MPI_Barrier(MPI_COMM_WORLD);
  setMatAndVecZero();

  MPI_Barrier(MPI_COMM_WORLD);
  solverPetscFluid->initialAssembly();  

  for(int t_itr=0; t_itr<timeMax; t_itr++)
  { 
    MPI_Barrier(MPI_COMM_WORLD);
    solverPetscFluid->setZeroValue();

    MPI_Barrier(MPI_COMM_WORLD);
    if(pulsatile_flow){
      if(t_itr > pulse_begin_itr){
        assignPulsatileBCs(t_itr);
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    applyBCs();
    
    MPI_Barrier(MPI_COMM_WORLD);
    timer = MPI_Wtime();

    for(int ic=0; ic<numOfElmGlobalFluid; ic++)
    {
      if(elmFluid[ic]->getSubdomainId() == myId)
      {
        Klocal.setZero();
        Flocal.setZero();

        switch(bd)
        {
          case BOUNDARY::XFEM:
            if(phiEXFluid[ic]>0.999){
              MatAssyUSNS(Klocal, Flocal, ic, t_itr);
            }else{
              XFEM_MatAssyUSNS(Klocal, Flocal, ic, t_itr);
            }
            break;

          case BOUNDARY::DARCY:
            if(phiVOFFluid[ic]>0.999){
              MatAssyUSNS(Klocal, Flocal, ic, t_itr);
            }else{
              Darcy_MatAssyUSNS(Klocal, Flocal, ic, t_itr);
            }
            break;

          default:
            PetscPrintf(MPI_COMM_WORLD, "undefined boundary method");
            exit(1);
        }
        solverPetscFluid->setValue(elmFluid[ic]->nodeForAssyBCsFluid, elmFluid[ic]->nodeForAssyFluid, Klocal, Flocal);
      
      }
    }

    timer = MPI_Wtime() - timer;
    PetscPrintf(MPI_COMM_WORLD, "\n ****** Time for matrix assembly = %f seconds ****** \n", timer);

    solverPetscFluid->currentStatus = ASSEMBLY_OK;
    
    MPI_Barrier(MPI_COMM_WORLD); 
    timer = MPI_Wtime();
    
    solverPetscFluid->solve();
  
    timer = MPI_Wtime() - timer;

    PetscPrintf(MPI_COMM_WORLD, "\n ****** Time for PETSc solver = %f seconds ****** \n", timer);
    
    for(ii=0; ii<numOfDofsGlobalFluid; ii++){
      SolnDataFluid.soln[assyForSolnFluid[ii]] = 0e0;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    VecScatterBegin(ctx, solverPetscFluid->solnVec, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, solverPetscFluid->solnVec, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);

    VecGetArray(vec_SEQ, &arrayTempSolnFluid);

    // update solution vector
    for(ii=0; ii<numOfDofsGlobalFluid; ii++){
      SolnDataFluid.soln[assyForSolnFluid[ii]] += arrayTempSolnFluid[ii];
    }
    
    VecRestoreArray(vec_SEQ, &arrayTempSolnFluid);

    postCaluculation_timeItr(t_itr);

    if(myId == 0){
      double time_now = t_itr * dt;
      printf("\n TIME = %f \n", time_now);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }

  VecScatterDestroy(&ctx);
  VecDestroy(&vec_SEQ);
}