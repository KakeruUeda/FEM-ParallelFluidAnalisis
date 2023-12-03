#include "FEM.h"

void FEM::postCaluculation(){

  if(myId>0) return;

  int ii, kk, nn, mm;

  vector<double> u(numOfNodeGlobal,0);
  vector<double> v(numOfNodeGlobal,0);
  vector<double> w(numOfNodeGlobal,0);
  vector<double> p(numOfNodeGlobal,0);

  for(ii=0; ii<numOfNodeGlobal; ii++){
      nn = nodeMap[ii];
      mm = nn*numOfDofsNode;
      u[ii] = SolnData.soln[mm];
      v[ii] = SolnData.soln[mm+1];
      w[ii] = SolnData.soln[mm+2];
      p[ii] = SolnData.soln[mm+3];
  }
        
  string vtiFile;
  vtiFile = "resutls.vti";
  export_vti_result(vtiFile,u,v,w,p);

  vtiFile = "resutls_2D.vti";
  export_vti_result_2D(vtiFile,u,v,p);


  double Q1 = 0e0;
  double Q2 = 0e0;

  for(kk=0; kk<nz+1; kk++){
    for(ii=0; ii<nx+1; ii++){
      Q1 += v[kk*(nx+1)*(ny+1)+ii]*dx*dy;
      Q2 += v[kk*(nx+1)*(ny+1)+(nx+1)*(ny)+ii]*dx*dy;
    }
  }

  printf("Q1 = %lf Q2 = %lf \n", Q1,Q2);

}