#include "FEM.h"

void FEM::Stokes()
{
  PetscPrintf(MPI_COMM_WORLD, "\n\n\n STEADY STOKES EQUATION \n");

  int  aa, bb, ee, ii, jj, kk, count, row, col, jpn, n1, n2, size1, size2;

  double  norm_rhs, timer;

  jpn = numOfNodeInElm*numOfDofsNode;
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
  
  MPI_Barrier(MPI_COMM_WORLD);
  applyBCs();

  MPI_Barrier(MPI_COMM_WORLD);
  PetscPrintf(MPI_COMM_WORLD, "\n ****** ELEMENT LOOP START ******* \n", timer);
  
  timer = MPI_Wtime();
  int i = 0;

  for(int ic=0;ic<numOfElmGlobalFluid;ic++)
  {
    if(elmFluid[ic]->getSubdomainId() == myId)
    {  
      Klocal.setZero();
      Flocal.setZero();

      VDOUBLE1D sdf_current(numOfNodeInElm,0e0);
      double sdf_extention = 2e0 * sqrt(dx * dx + dy * dy);
      double sdf_ave = 0e0;
      bool xfem_bd;
  
      for(int i=0;i<numOfNodeInElm;i++){
        sdf_current[i] = sdfFluid[elmFluid[ic]->nodeNumsPrevFluid[i]];
        sdf_ave += sdf_current[i];
      }

      sdf_ave = sdf_ave / numOfNodeInElm;
      
      if(fabs(sdf_ave) < sdf_extention){
        xfem_bd = true;
      }else{
        xfem_bd = false;
      }

      switch(bd){
        case BOUNDARY::XFEM:
          if(xfem_bd){
            XFEM_MatAssySTT(ic,Klocal,Flocal);
          }else{
            MatAssySTT(ic,Klocal,Flocal);
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

  PetscPrintf(MPI_COMM_WORLD, "\n Time for matrix assembly = %f seconds \n", timer);

  PetscPrintf(MPI_COMM_WORLD, "\n ****** ELEMENT LOOP END ****** \n\n", timer);

  VecAssemblyBegin(solverPetscFluid->rhsVec);
  VecAssemblyEnd(solverPetscFluid->rhsVec);

  solverPetscFluid->currentStatus = ASSEMBLY_OK;
  
  MPI_Barrier(MPI_COMM_WORLD);

  PetscPrintf(MPI_COMM_WORLD, "\n ****** SOLVING LINEAR EQUATION ****** \n", timer);
  timer = MPI_Wtime();

  solverPetscFluid->factoriseAndSolve();
  
  MPI_Barrier(MPI_COMM_WORLD);
  timer = MPI_Wtime() - timer;

  PetscPrintf(MPI_COMM_WORLD, "\n Time for PETSc solver = %f seconds \n", timer);
  
  //get the solution vector 
  VecScatterBegin(ctx, solverPetscFluid->solnVec, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(ctx, solverPetscFluid->solnVec, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);

  VecGetArray(vec_SEQ, &arrayTempSolnFluid);

  for(ii=0; ii<numOfDofsGlobalFluid; ii++){
    SolnDataFluid.soln[assyForSolnFluid[ii]]   +=  arrayTempSolnFluid[ii];
  }

  MPI_Barrier(MPI_COMM_WORLD);
  
  VecRestoreArray(vec_SEQ, &arrayTempSolnFluid);
  
  VecScatterDestroy(&ctx);
  VecDestroy(&vec_SEQ);
  
}



