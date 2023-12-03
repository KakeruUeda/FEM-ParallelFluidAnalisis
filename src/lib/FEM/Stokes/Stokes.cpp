#include "FEM.h"
using namespace std;

void FEM::Stokes(){

  PetscPrintf(MPI_COMM_WORLD, " Solve Stokes equation \n\n\n");

  int  stepsCompleted=0;
  int  aa, bb, ee, ii, jj, kk, count, row, col, jpn, n1, n2, size1, size2;

  double  norm_rhs, fact, fact1, fact2, timerVal;

  VectorXd  reacVec(numOfNodeGlobal*numOfDofsNode);
  jpn = numOfNodeInElm*numOfDofsNode;
  VectorXd  Flocal(jpn);
  MatrixXd  Klocal(jpn, jpn);
  
  PetscScalar *arrayTempSoln;

  Vec            vec_SEQ;
  VecScatter     ctx;
  VecScatterCreateToAll(solverPetsc->solnVec, &ctx, &vec_SEQ);
  MPI_Barrier(MPI_COMM_WORLD);

  solverPetsc->zeroMtx();
  reacVec.setZero();

  MPI_Barrier(MPI_COMM_WORLD);
  
  assignBoundaryConditions();

  for(int ic=0;ic<numOfElmGlobal;ic++){
    if(elm[ic]->getSubdomainId() == myId){

      Klocal.setZero();
      Flocal.setZero();
      
      calcStokesMatrix(ic,Klocal,Flocal);
      
      size1 = elm[ic]->nodeForAssyBCs.size();
      applyBoundaryConditions(Flocal, size1, ic);
      solverPetsc->assembleMatrixAndVectorSerial(elm[ic]->nodeForAssyBCs, elm[ic]->nodeForAssy, Klocal, Flocal);
    
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  
  VecAssemblyBegin(solverPetsc->rhsVec);
  VecAssemblyEnd(solverPetsc->rhsVec);

  //VecNorm(solverPetsc->rhsVec, NORM_2, &norm_rhs);
  solverPetsc->currentStatus = ASSEMBLY_OK;
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  PetscPrintf(MPI_COMM_WORLD, "Assembly done. Solving the matrix system. \n");
  
  solverPetsc->factoriseAndSolve();

  PetscPrintf(MPI_COMM_WORLD, "Solved. \n");
  
  MPI_Barrier(MPI_COMM_WORLD);

  /////////////////////////////////////////////////////////////////////////////
  // get the solution vector onto all the processors
  /////////////////////////////////////////////////////////////////////////////
  
  VecScatterBegin(ctx, solverPetsc->solnVec, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(ctx, solverPetsc->solnVec, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);

  VecGetArray(vec_SEQ, &arrayTempSoln);

  // update solution vector
  for(ii=0; ii<numOfDofsGlobal; ii++){
    SolnData.soln[assyForSoln[ii]]   +=  arrayTempSoln[ii];
  }

  MPI_Barrier(MPI_COMM_WORLD);
  
  VecRestoreArray(vec_SEQ, &arrayTempSoln);
  
  getSolution();
  
  VecScatterDestroy(&ctx);
  VecDestroy(&vec_SEQ);
}


void FEM::assignBoundaryConditions(){

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

void FEM::applyBoundaryConditions(VectorXd& Flocal, const int size, const int ic){

  int ii, value;
  for(ii=0; ii<size; ii++){
    value = DirichletBCs[elm[ic]->globalDOFnums[ii]];
    if(elm[ic]->nodeForAssyBCs[ii] == -1) {
      Flocal(ii) = value;    
    }
  }

}


void FEM::getSolution(){

  int ii;
  vector<double> u(numOfNodeGlobal,0);
  vector<double> v(numOfNodeGlobal,0);
  vector<double> w(numOfNodeGlobal,0);
  vector<double> p(numOfNodeGlobal,0);

  for(ii=0; ii<numOfNodeGlobal; ii++){
      int nnn = nodeMap[ii];
      int kkk = nnn*numOfDofsNode;
      u[ii] = SolnData.soln[kkk];
      v[ii] = SolnData.soln[kkk+1];
      w[ii] = SolnData.soln[kkk+2];
      p[ii] = SolnData.soln[kkk+3];
  }
        
  string vtiFile;
  vtiFile = "resutls.vti";
  export_vti_result(vtiFile,u,v,w,p);

  vtiFile = "resutls_2D.vti";
  export_vti_result_2D(vtiFile,u,v,p);

}