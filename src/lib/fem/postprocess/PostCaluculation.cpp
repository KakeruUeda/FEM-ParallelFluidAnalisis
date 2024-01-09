#include "FEM.h"

void FEM::postCaluculation_itr(const int loop){

  int ii, kk, nn, mm;

  for(ii=0; ii<numOfNodeGlobalFluid; ii++){
    nn = nodeMapFluid[ii];
    mm = nn*numOfDofsNode;
    if(loop < NRitr_initial){
      uFluid[ii] += SolnDataFluid.soln[mm] * relaxationParam_initial;
      vFluid[ii] += SolnDataFluid.soln[mm+1] * relaxationParam_initial;
      wFluid[ii] += SolnDataFluid.soln[mm+2] * relaxationParam_initial;
      pFluid[ii] += SolnDataFluid.soln[mm+3] * relaxationParam_initial;
    }else{
      uFluid[ii] += SolnDataFluid.soln[mm] * relaxationParam;
      vFluid[ii] += SolnDataFluid.soln[mm+1] * relaxationParam;
      wFluid[ii] += SolnDataFluid.soln[mm+2] * relaxationParam;
      pFluid[ii] += SolnDataFluid.soln[mm+3] * relaxationParam;
    }

    u[sortNode[ii]] = uFluid[ii];
    v[sortNode[ii]] = vFluid[ii];
    w[sortNode[ii]] = wFluid[ii];   
    p[sortNode[ii]] = pFluid[ii];
  }
  if(myId==0){
    string vtiFile;
    
    vtiFile = outputDir + "/SNS_"+to_string(loop)+".vti";
    export_vti_result(vtiFile,u,v,w,p);
  
    vtiFile = outputDir + "/SNS_2D_"+to_string(loop)+".vti";
    export_vti_result_2D(vtiFile,u,v,p);
  }

}


void FEM::postCaluculation_timeItr(const int t_itr){

  int ii, kk, nn, mm;

  for(ii=0; ii<numOfNodeGlobalFluid; ii++)
  {  
    nn = nodeMapFluid[ii];
    mm = nn * numOfDofsNode;

    uf[1][ii] = 0e0; 
    vf[1][ii] = 0e0;
    wf[1][ii] = 0e0;
    pf[1][ii] = 0e0;

    uf[1][ii] = uf[0][ii];
    vf[1][ii] = vf[0][ii];
    wf[1][ii] = wf[0][ii];
    pf[1][ii] = pf[0][ii];

    uf[0][ii] = 0e0; 
    vf[0][ii] = 0e0;
    wf[0][ii] = 0e0;
    pf[0][ii] = 0e0;

    uf[0][ii] = SolnDataFluid.soln[mm];
    vf[0][ii] = SolnDataFluid.soln[mm+1];
    wf[0][ii] = SolnDataFluid.soln[mm+2];
    pf[0][ii] = SolnDataFluid.soln[mm+3];

    u[sortNode[ii]] = 0e0;
    v[sortNode[ii]] = 0e0;
    w[sortNode[ii]] = 0e0;
    p[sortNode[ii]] = 0e0;

    u[sortNode[ii]] = uf[0][ii];
    v[sortNode[ii]] = vf[0][ii];
    w[sortNode[ii]] = wf[0][ii];
    p[sortNode[ii]] = pf[0][ii];
  }


  if(myId == 0){
    string vtiFile;
    
    vtiFile = outputDir + "/USNS_"+to_string(t_itr)+".vti";
    export_vti_result(vtiFile,u,v,w,p);
  
    vtiFile = outputDir + "/USNS_2D_"+to_string(t_itr)+".vti";
    export_vti_result_2D(vtiFile,u,v,p);
  }
}

void FEM::postCaluculation(){

  if(myId>0) return;

  int ii, kk, nn, mm;

  for(ii=0; ii<numOfNodeGlobalFluid; ii++)
  {
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