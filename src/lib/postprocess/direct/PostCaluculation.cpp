#include "FEM.h"


void FEM::postCaluculation_itr(const int loop){

  if(myId>0) return;

  int ii, kk, nn, mm;
  double relaxationParam = 5e-1;

  for(ii=0; ii<numOfNodeGlobalFluid; ii++){
    nn = nodeMapFluid[ii];
    mm = nn*numOfDofsNode;
    uFluid[ii] = SolnDataFluid.soln[mm] * relaxationParam;
    vFluid[ii] = SolnDataFluid.soln[mm+1] * relaxationParam;
    wFluid[ii] = SolnDataFluid.soln[mm+2] * relaxationParam;
    pFluid[ii] = SolnDataFluid.soln[mm+3] * relaxationParam;
    if(uFluid[ii] != 0) cout << uFluid[ii] << endl;

    u[sortNode[ii]] = uFluid[ii];
    v[sortNode[ii]] = vFluid[ii];
    w[sortNode[ii]] = wFluid[ii];   
    p[sortNode[ii]] = pFluid[ii];
  }
 
  string vtiFile;
  
  vtiFile = outputDir + "/NS_itr_"+to_string(loop)+".vti";
  export_vti_result(vtiFile,u,v,w,p);

  vtiFile = outputDir + "/NS_itr_"+to_string(loop)+"_2D.vti";
  export_vti_result_2D(vtiFile,u,v,p);
}

void FEM::postCaluculation(){

  if(myId>0) return;

  int ii, kk, nn, mm;


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
  
  vtiFile = outputDir + "/result.vti";
  export_vti_result(vtiFile,u,v,w,p);

  vtiFile = outputDir + "/result_2D.vti";
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