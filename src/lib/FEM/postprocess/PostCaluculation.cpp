#include "FEM.h"

void FEM::postCaluculation(){

  if(myId>0) return;

  int ii, kk, nn, mm;

  vector<double> uFluid(numOfNodeGlobalFluid,0);
  vector<double> vFluid(numOfNodeGlobalFluid,0);
  vector<double> wFluid(numOfNodeGlobalFluid,0);
  vector<double> pFluid(numOfNodeGlobalFluid,0);

  vector<double> u(numOfNodeGlobal,0);
  vector<double> v(numOfNodeGlobal,0);
  vector<double> w(numOfNodeGlobal,0);
  vector<double> p(numOfNodeGlobal,0);


  for(ii=0; ii<numOfNodeGlobalFluid; ii++){
    nn = nodeMapFluid[ii];
    mm = nn*numOfDofsNode;
    uFluid[ii] = SolnDataFluid.soln[mm];
    vFluid[ii] = SolnDataFluid.soln[mm+1];
    wFluid[ii] = SolnDataFluid.soln[mm+2];
    pFluid[ii] = SolnDataFluid.soln[mm+3];

    u[sortNode[ii]] = uFluid[ii];
    v[sortNode[ii]] = vFluid[ii];
    w[sortNode[ii]] = wFluid[ii];   
    p[sortNode[ii]] = pFluid[ii];
  }

        
  string vtiFile;
  vtiFile = "resutls.vti";
  export_vti_result(vtiFile,u,v,w,p);

  vtiFile = "resutls_2D.vti";
  export_vti_result_2D(vtiFile,u,v,p);


  double QINLET = 0e0;
  double QHALF = 0e0;
  double QOUTLET = 0e0;

  for(kk=0; kk<nz+1; kk++){
    for(ii=0; ii<nx+1; ii++){
      QINLET += v[kk*(nx+1)*(ny+1)+ii]*dx*dz;
      QHALF += v[kk*(nx+1)*(ny+1)+(nx+1)*(ny/2)+ii]*dx*dz;
      QOUTLET += v[kk*(nx+1)*(ny+1)+(nx+1)*(ny)+ii]*dx*dz;
    }
  }

  printf("QINLET = %lf QHALF = %lf QOUTLET = %lf \n", QINLET,QHALF,QOUTLET);

}