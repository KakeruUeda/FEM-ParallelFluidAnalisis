#include "FEM.h"
using namespace std;

void FEM::Stokes()
{
  PetscPrintf(MPI_COMM_WORLD, " Solve Stokes equation \n\n\n");

  int  stepsCompleted=0;
  int  aa, bb, ee, ii, jj, kk, count, row, col, ind, n1, n2, size1, size2;

  double  norm_rhs, fact, fact1, fact2, timerVal;

  VectorXd  reacVec(numOfNodeGlobal*ndof);
  ind = numOfNodeInElm*ndof;
  VectorXd  Flocal(ind);
  MatrixXd  Klocal(ind, ind);
  
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
    cout << "ic = " << ic << " elm[ic]->getSubdomainId() = " << elm[ic]->getSubdomainId() << endl;
    if(elm[ic]->getSubdomainId() == this_mpi_proc)
    {
      Klocal.setZero();
      Flocal.setZero();
      calcStokesMatrix(ic,Klocal,Flocal);

      /*
      if(this_mpi_proc == 0){
        ofstream outKlocal0("Klocaltmp.dat");
        ofstream outFlocal0("Flocaltmp.dat");
        
        outKlocal0 << Klocal << endl;
        outKlocal0 << endl;
        
        outFlocal0 << Flocal << endl;
        outFlocal0 << endl;
        
        outKlocal0.close();
        outFlocal0.close();
        //exit(1);
      }
      */

      size1 = elm[ic]->forAssyVec.size();
      for(ii=0; ii<size1; ii++)
      {
        aa = elm[ic]->forAssyVec[ii];
        fact = SolnData.solnApplied[elm[ic]->globalDOFnums[ii]];
        //Flocal(ii) = fact;
        if(aa == -1) // this DOF has a prescribed value
        {
          fact = SolnData.solnApplied[elm[ic]->globalDOFnums[ii]];
          //Flocal(ii) = fact;
          // check if fact is zero. We don't need to
          // execute the for loop if fact is zero.
          if( abs(fact) > 1.0e-10)
          {
            //Flocal(ii) = fact;
            for(jj=0; jj<size1; jj++)
            {
              if(elm[ic]->forAssyVec[jj] != -1)
              {
                //cout << "fact = " << fact << endl;
                Flocal(jj) -= Klocal(jj,ii) * fact;
              }
            }
          }
        }
      }    
      solverPetsc->assembleMatrixAndVectorSerial(elm[ic]->forAssyVec, Klocal, Flocal);
      /*
      if(this_mpi_proc == 0){
        ofstream Assy("forAssyVec0.dat"); 
        for(ic=0; ic<numOfElmGlobal; ic++){
          int size = elm[ic]->forAssyVec.size();
          for(ii=0; ii<size; ii++){
            if(ii != 0 && ii%4 == 0) Assy << " ";
            Assy << elm[ic]->forAssyVec[ii] << " ";
          } 
          Assy << endl;
        }
        Assy.close();
      }
      if(this_mpi_proc == 1){
        ofstream Assy1("forAssyVec1.dat"); 
        for(ic=0; ic<numOfElmGlobal; ic++){
          int size = elm[ic]->forAssyVec.size();
          for(ii=0; ii<size; ii++){
            if(ii != 0 && ii%4 == 0) Assy1 << " ";
            Assy1 << elm[ic]->forAssyVec[ii] << " ";
          } 
          Assy1 << endl;
        }
        Assy1.close();
      }
      */
      
      /*
      if(this_mpi_proc == 0){
        ofstream outKlocal0("Klocal.dat");
        ofstream outFlocal0("Flocal.dat");
        
        outKlocal0 << Klocal << endl;
        outKlocal0 << endl;
        
        outFlocal0 << Flocal << endl;
        outFlocal0 << endl;
        
        outKlocal0.close();
        outFlocal0.close();
        exit(1);
      }
      */
    }
  }

  //MPI_Barrier(MPI_COMM_WORLD);
  for(ii=0; ii<numOfBdNode; ii++){
    n1 = int(DirichletBCs[ii][0]);
    n2 = int(DirichletBCs[ii][1]);

    jj = n1*ndof+n2;

    //SolnData.soln[jj]  += SolnData.solnApplied[jj];
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
  for(ii=0; ii<ntotdofs_global; ii++)
  {
    SolnData.soln[assyForSoln[ii]]   +=  arrayTempSoln[ii];
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  VecRestoreArray(vec_SEQ, &arrayTempSoln);
  
  vector<double> u(numOfNodeGlobal,0);
  vector<double> v(numOfNodeGlobal,0);
  vector<double> w(numOfNodeGlobal,0);
  vector<double> p(numOfNodeGlobal,0);

  for(ii=0; ii<numOfNodeGlobal; ii++){
    kk = ii*ndof;
    u[ii] = SolnData.soln[kk];
    v[ii] = SolnData.soln[kk+1];
    w[ii] = SolnData.soln[kk+2];
    p[ii] = SolnData.soln[kk+3];
    if(u[ii] != 0 || v[ii] != 0 || w[ii] != 0 || p[ii] != 0){
      cout << "i = " << ii <<  " u = " << u[ii] << " v = " << v[ii] << " w = " << w[ii] << " p = "<< p[ii] << endl;
    }
  }
  
  if(this_mpi_proc == 0){
    ofstream res("results.dat");
    for(int i=0; i<numOfNodeGlobal; i++){
      res << "i = " << i <<  " u = " << u[i] << " v = " << v[i] << " w = " << w[i] << " p = "<< p[i] << endl;
    }
    res.close();
  }
  if(this_mpi_proc == 0){
    ofstream res2("resultsOrig.dat");
    vector<double> vec(4,0);
    for(int i=0; i<numOfNodeGlobal; i++){
      int nnn = node_map_get_new[i];
      int kkk = nnn*ndof;

      vec[0] = SolnData.soln[kkk];
      vec[1] = SolnData.soln[kkk+1];
      vec[2] = SolnData.soln[kkk+2];
      vec[3] = SolnData.soln[kkk+3];

      res2 << "i = " << i <<  " u = " << vec[0] << " v = " << vec[1] << " w = " << vec[2] << " p = "<< vec[3] << endl;
    }
    res2.close();
    }
        
    string vtiFile;
    vtiFile = "resutls.vti";
    export_vti_result(vtiFile,u,v,w,p);

    vtiFile = "resutls_2D.vti";
    export_vti_result_2D(vtiFile,u,v,p);
}


void FEM::assignBoundaryConditions()
{
    int ii, n1, n2, ind;
    for(ii=0; ii<numOfBdNode; ii++)
    {
        n1 = int(DirichletBCs[ii][0]);
        n2 = int(DirichletBCs[ii][1]);
        ind = n1*ndof+n2;
        //SolnData.soln[ind] = DirichletBCs[ii][2] * timeFact;
        SolnData.solnApplied[ind] = DirichletBCs[ii][2];
        //cout << ii << '\t' <<  SolnData.solnApplied[ind] << endl;

        //solnApplied[ind] = analy.computeValue(n2, node_coords[n1][0], node_coords[n1][1], 0.0, timeNow) - soln[ind];
    }
    //exit(1);
    //printVector(solnApplied);

    return;
}