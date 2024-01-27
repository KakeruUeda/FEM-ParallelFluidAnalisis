#include "VariationalDataAssimilation.h"

void VarDA::adjointSteadyNavierStokes(VDOUBLE2D &externalForce)
{
  PetscPrintf(MPI_COMM_WORLD, "\n\n\n ADJOINT STEADY NAVIER STOKES EQUATION \n");

  int  aa, bb, ee, ii, jj, kk, count, row, col, jpn, jpn2, n1, n2, size1, size2;
  double  norm_rhs, timer;

  double norm, norm0;
  double tmp = 1e12;

  jpn = numOfNodeInElm * numOfDofsNode;
  VectorXd  Flocal(jpn);
  MatrixXd  Klocal(jpn, jpn);

  jpn2 = numOfNodeInElm * numOfDofsNode + 4 * 3;
  VectorXd FlocalBd(jpn2);
  MatrixXd KlocalBd(jpn2, jpn2);
  
  PetscScalar *arrayTempSolnAdjointFluid;

  Vec            vec_SEQ;
  VecScatter     ctx;
  
  VecScatterCreateToAll(adjointSolverPetscFluid->solnVec, &ctx, &vec_SEQ);

  MPI_Barrier(MPI_COMM_WORLD);
  assignBCsAdjoint();
  
  MPI_Barrier(MPI_COMM_WORLD);
  setAdjointMatAndVecZero();

  MPI_Barrier(MPI_COMM_WORLD);
  adjointSolverPetscFluid->initialAssembly();  

  MPI_Barrier(MPI_COMM_WORLD);
  applyBCsAdjoint();
      cout << "2" << endl;

  for(int ic=0; ic<numOfElmGlobalFluid; ic++)
  {
    if(elmFluid[ic]->getSubdomainId() == myId)
    {
      Klocal.setZero();
      Flocal.setZero();
      KlocalBd.setZero();
      FlocalBd.setZero();

      if(!control_elm_prev_type[ic])
      {
        switch(bd)
        {
          case BOUNDARY::XFEM:
            if(phiEXFluid[ic]>0.999){
              AdjointMatAssySNS(ic, Klocal, Flocal, externalForce);
            }else{
              PetscPrintf(MPI_COMM_WORLD, "XFEM not yet implemented");
              exit(1);
            }
            break;

          case BOUNDARY::DARCY:
            if(phiVOFFluid[ic]>0.999){
              AdjointMatAssySNS(ic, Klocal, Flocal, externalForce);
            }else{
              Darcy_AdjointMatAssySNS(ic, Klocal, Flocal, externalForce);
            }
            break;

          default:
            PetscPrintf(MPI_COMM_WORLD, "undefined boundary method");
            exit(1);
        }
        adjointSolverPetscFluid->setValue(elmFluid[ic]->nodeForAssyBCsAdjointFluid, elmFluid[ic]->nodeForAssyAdjointFluid, Klocal, Flocal);
      }
      else
      {
        switch(bd)
        {
          case BOUNDARY::XFEM:
            if(phiEXFluid[ic]>0.999){
              AdjointMatAssySNS(ic, KlocalBd, FlocalBd, externalForce);
            }else{
              PetscPrintf(MPI_COMM_WORLD, "XFEM not yet implemented");
              exit(1);
            }
            break;

          case BOUNDARY::DARCY:
            if(phiVOFFluid[ic]>0.999){
              AdjointBdMatAssySNS(ic, KlocalBd, FlocalBd, externalForce);
            }else{
              Darcy_AdjointBdMatAssySNS(ic, KlocalBd, FlocalBd, externalForce);
            }
            break;

          default:
            PetscPrintf(MPI_COMM_WORLD, "undefined boundary method");
            exit(1);
        }
        calcBoundaryIntegral(ic, KlocalBd, FlocalBd);
        adjointSolverPetscFluid->setValue(elmFluid[ic]->nodeForAssyBCsAdjointFluid, elmFluid[ic]->nodeForAssyAdjointFluid, KlocalBd, FlocalBd);
      }    
    
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  adjointSolverPetscFluid->solve();

  for(ii=0; ii<numOfDofsAdjointGlobalFluid; ii++){
    SolnDataFluid.solnAdjoint[assyForSolnAdjointFluid[ii]] = 0e0;
  }

  VecScatterBegin(ctx, adjointSolverPetscFluid->solnVec, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(ctx, adjointSolverPetscFluid->solnVec, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);
  VecGetArray(vec_SEQ, &arrayTempSolnAdjointFluid);

  for(ii=0; ii<numOfDofsAdjointGlobalFluid; ii++){
    SolnDataFluid.solnAdjoint[assyForSolnAdjointFluid[ii]] += arrayTempSolnAdjointFluid[ii];
  }

  MPI_Barrier(MPI_COMM_WORLD); 
  VecRestoreArray(vec_SEQ, &arrayTempSolnAdjointFluid);

  scatterPysicalvVariables();

  return;
}

void VarDA::scatterPysicalvVariables()
{
  for(int in=0; in<numOfNodeGlobalFluid; in++)
  {
    int nn = nodeDofsMapFluid[in];
    
    if(!control_node_type[in])
    {
      u_adjoint_fluid[in] = SolnDataFluid.solnAdjoint[nn];
      v_adjoint_fluid[in] = SolnDataFluid.solnAdjoint[nn+1];
      w_adjoint_fluid[in] = SolnDataFluid.solnAdjoint[nn+2];
      p_adjoint_fluid[in] = SolnDataFluid.solnAdjoint[nn+3];
      
      u_adjoint[sortNode[in]] = u_adjoint_fluid[in];
      v_adjoint[sortNode[in]] = v_adjoint_fluid[in];
      w_adjoint[sortNode[in]] = w_adjoint_fluid[in];
      p_adjoint[sortNode[in]] = p_adjoint_fluid[in];
    }
    else
    {
      u_adjoint_fluid[in] = SolnDataFluid.solnAdjoint[nn];
      v_adjoint_fluid[in] = SolnDataFluid.solnAdjoint[nn+1];
      w_adjoint_fluid[in] = SolnDataFluid.solnAdjoint[nn+2];
      p_adjoint_fluid[in] = SolnDataFluid.solnAdjoint[nn+3];
      
      lambda_u_fluid[in]  = SolnDataFluid.solnAdjoint[nn+4];
      lambda_v_fluid[in]  = SolnDataFluid.solnAdjoint[nn+5];
      lambda_w_fluid[in]  = SolnDataFluid.solnAdjoint[nn+6];
      
      u_adjoint[sortNode[in]] = u_adjoint_fluid[in];
      v_adjoint[sortNode[in]] = v_adjoint_fluid[in];
      w_adjoint[sortNode[in]] = w_adjoint_fluid[in];
      p_adjoint[sortNode[in]] = p_adjoint_fluid[in];
      
      lambda_u[sortNode[in]]  = lambda_u_fluid[in];
      lambda_v[sortNode[in]]  = lambda_v_fluid[in];
      lambda_w[sortNode[in]]  = lambda_w_fluid[in];
    }

  }

  return;
}
