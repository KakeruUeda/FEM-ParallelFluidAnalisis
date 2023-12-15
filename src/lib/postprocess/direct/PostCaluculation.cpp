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
    mm = nn*numOfDofsNode;
    
    uFluid[ii] = SolnDataFluid.soln[mm];
    vFluid[ii] = SolnDataFluid.soln[mm+1];
    wFluid[ii] = SolnDataFluid.soln[mm+2];
    pFluid[ii] = SolnDataFluid.soln[mm+3];

    uf[t_itr][ii] = SolnDataFluid.soln[mm];
    vf[t_itr][ii] = SolnDataFluid.soln[mm+1];
    wf[t_itr][ii] = SolnDataFluid.soln[mm+2];
    pf[t_itr][ii] = SolnDataFluid.soln[mm+3];

    u[sortNode[ii]] = uf[t_itr][ii];
    v[sortNode[ii]] = vf[t_itr][ii];
    w[sortNode[ii]] = wf[t_itr][ii];   
    p[sortNode[ii]] = pf[t_itr][ii];
  }

  if(t_itr >= 2){
    uf[t_itr-2].clear();
    vf[t_itr-2].clear();
    wf[t_itr-2].clear();
    pf[t_itr-2].clear();
  }


  if(myId==0){
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