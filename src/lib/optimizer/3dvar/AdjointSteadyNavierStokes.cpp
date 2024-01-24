#include "VariationalDataAssimilation.h"

void VarDA::adjoint_SteadyNavierStokes(VDOUBLE2D &externalForce)
{
  PetscPrintf(MPI_COMM_WORLD, "\n\n\n ADJOINT STEADY NAVIER STOKES EQUATION \n");

  int  aa, bb, ee, ii, jj, kk, count, row, col, jpn, jpn2, n1, n2, size1, size2;
  double  norm_rhs, timer;

  double norm, norm0;
  double tmp = 1e12;

  jpn = numOfNodeInElm * numOfDofsNode;
  VectorXd  Flocal(jpn);
  MatrixXd  Klocal(jpn, jpn);

  jpn2 = 4 * numOfDofsNode;
  VectorXd FlocalBd(jpn2);
  MatrixXd KlocalBd(jpn2, jpn2);
  
  PetscScalar *arrayTempSolnFluid;

  Vec            vec_SEQ;
  VecScatter     ctx;
  
  VecScatterCreateToAll(adjointSolverPetscFluid->solnVec, &ctx, &vec_SEQ);

  MPI_Barrier(MPI_COMM_WORLD);
  assignBCs();
  
  MPI_Barrier(MPI_COMM_WORLD);
  setAdjointMatAndVecZero();

  MPI_Barrier(MPI_COMM_WORLD);
  adjointSolverPetscFluid->initialAssembly();  

  MPI_Barrier(MPI_COMM_WORLD);
  setNRInitialValue(); 

  MPI_Barrier(MPI_COMM_WORLD);
  applyBCs();

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
            AdjointMatAssySNS(ic,Klocal,Flocal);
          }else{
            PetscPrintf(MPI_COMM_WORLD, "XFEM not yet implemented");
            exit(1);
          }
          break;

        case BOUNDARY::DARCY:
          if(phiVOFFluid[ic]>0.999){
            AdjointMatAssySNS(ic,Klocal,Flocal);
          }else{
            Darcy_AdjointMatAssySNS(ic,Klocal,Flocal);
          }
          break;

        default:
          PetscPrintf(MPI_COMM_WORLD, "undefined boundary method");
          exit(1);
      }
      adjointSolverPetscFluid->setValue(elmFluid[ic]->nodeForAssyBCsFluid, elmFluid[ic]->nodeForAssyFluid, Klocal, Flocal);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for(int ic=0; ic<control_elm_node.size(); ic++)
  {
    FlocalBd.setZero();
    KlocalBd.setZero();
  }

  return;
}