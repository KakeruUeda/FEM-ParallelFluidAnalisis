#include "FEM.h"
using namespace std;

void FEM::stokes()
{
  int  aa, bb, ee, ii, jj, kk, count, row, col, jpn, n1, n2, size1, size2;

  double  norm_rhs, timer;

  VectorXd  reacVec(numOfNodeGlobal*numOfDofsNode);
  jpn = numOfNodeInElm*numOfDofsNode;
  VectorXd  Flocal(jpn);
  MatrixXd  Klocal(jpn, jpn);
  
  PetscScalar *arrayTempSoln;

  Vec            vec_SEQ;
  VecScatter     ctx;
  VecScatterCreateToAll(solverPetsc->solnVec, &ctx, &vec_SEQ);
  MPI_Barrier(MPI_COMM_WORLD);
  
  PetscPrintf(MPI_COMM_WORLD, " First assembly start \n");
  timer = MPI_Wtime();
  
  solverPetsc->zeroMtx();
  reacVec.setZero();
  
  PetscPrintf(MPI_COMM_WORLD, " First assembly done \n\n");
  timer = MPI_Wtime() - timer;
  PetscPrintf(MPI_COMM_WORLD, " ***** Time for first assembly  = %f seconds ***** \n\n", timer);
  
  assignBoundaryConditions();
  
  MPI_Barrier(MPI_COMM_WORLD);

  PetscPrintf(MPI_COMM_WORLD, " Element loop start \n");
  
  timer = MPI_Wtime();

  for(int ic=0;ic<numOfElmGlobal;ic++){
    if(elm[ic]->getSubdomainId() == myId){
      
      Klocal.setZero();
      Flocal.setZero();
      
      //if(phi[ic]>1e-5){
      //cout << "2" << endl;
      //  calcStokesMatrix(ic,Klocal,Flocal);
      //}else{
      // calcStokesMatrix(ic,Klocal,Flocal);
      //}
      //else{
      //  continue;
      //  for(int i=0; i<elm[ic]->nodeForAssyBCs.size(); i++){
      //    if(((i+1)%4 != 0) && (elm[ic]->nodeForAssyBCs[i] != -1)){
      //      cout << " ic = " << ic << " i = " << i <<  " elm[ic]->nodeForAssyBCs[i] = " << elm[ic]->nodeForAssyBCs[i] << " myId = " << myId << endl;
      //    }
      //  }
      //
      
      if(phiEX[ic]>0.999){
        calcStokesMatrix(ic,Klocal,Flocal);
      }else if(phiEX[ic]>1e-5){
        calcStokesMatrixXFEM(ic,Klocal,Flocal);
      }else{
        calcStokesMatrix(ic,Klocal,Flocal);
      }

      size1 = elm[ic]->nodeForAssyBCs.size();
      
      applyBoundaryConditions(Flocal, size1, ic);
      solverPetsc->assembleMatrixAndVectorSerial(elm[ic]->nodeForAssyBCs, elm[ic]->nodeForAssy, Klocal, Flocal);
    }
  }
    
  MPI_Barrier(MPI_COMM_WORLD);

  timer = MPI_Wtime() - timer;
  //computerTimeAssembly += timerVal;
  PetscPrintf(MPI_COMM_WORLD, " Element loop done \n\n");
  PetscPrintf(MPI_COMM_WORLD, " ***** Time for matrix assembly = %f seconds ***** \n\n", timer);

  VecAssemblyBegin(solverPetsc->rhsVec);
  VecAssemblyEnd(solverPetsc->rhsVec);

  //VecNorm(solverPetsc->rhsVec, NORM_2, &norm_rhs);
  solverPetsc->currentStatus = ASSEMBLY_OK;
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  timer = MPI_Wtime();
  solverPetsc->factoriseAndSolve();
  
  MPI_Barrier(MPI_COMM_WORLD);
  timer = MPI_Wtime() - timer;
  //computerTimeSolver += timerVal;
  PetscPrintf(MPI_COMM_WORLD, " ***** Time for PETSc solver = %f seconds ***** \n\n", timer);
  
  ////// gather the solution vectors //////
  VecScatterBegin(ctx, solverPetsc->solnVec, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(ctx, solverPetsc->solnVec, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);

  VecGetArray(vec_SEQ, &arrayTempSoln);

  for(ii=0; ii<numOfDofsGlobal; ii++){
    SolnData.soln[assyForSoln[ii]]   +=  arrayTempSoln[ii];
  }

  MPI_Barrier(MPI_COMM_WORLD);
  
  VecRestoreArray(vec_SEQ, &arrayTempSoln);
  
  VecScatterDestroy(&ctx);
  VecDestroy(&vec_SEQ);
}


void FEM::assignBoundaryConditions()
{
  int ii, n1, n2, n3;
  DirichletBCs.resize(numOfDofsNode*numOfNodeGlobal);
  for(ii=0; ii<numOfBdNode; ii++){
    n1 = int(DirichletBCs_tmp[ii][0]);
    n2 = int(DirichletBCs_tmp[ii][1]);
    n3 = n1*numOfDofsNode+n2;
    DirichletBCs[n3] = DirichletBCs_tmp[ii][2];
  }
  return;
}

void FEM::applyBoundaryConditions(VectorXd& Flocal, const int size, const int ic)
{
  int ii;
  for(ii=0; ii<size; ii++){
    if(elm[ic]->nodeForAssyBCs[ii] == -1) {
      Flocal(ii) = DirichletBCs[elm[ic]->globalDOFnums[ii]];    
    }
  }
}