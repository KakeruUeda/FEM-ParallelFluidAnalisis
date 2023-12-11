#include "FEM.h"
using namespace std;

void FEM::readInput(){

  string str,base_label,label;

  int tmp[3];
  double dummy[3];

  /////domain/////
  base_label = "/Domain";
  label = base_label + "/nx";
  if ( !tp.getInspectedVector(label,tmp,3)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  nx=tmp[0];
  ny=tmp[1];
  nz=tmp[2];

  label = base_label + "/Lx";
  if ( !tp.getInspectedVector(label,dummy,3)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  Lx=dummy[0];
  Ly=dummy[1];
  Lz=dummy[2];
  
  dx=Lx/nx;
  dy=Ly/ny;
  dz=Lz/nz; 
  numOfNodeGlobal=(nx+1)*(ny+1)*(nz+1);
  numOfElmGlobal=nx*ny*nz;
  numOfNodeInElm=8;

  /////param/////
  base_label = "/Parameter";
  label = base_label + "/mu";
  if ( !tp.getInspectedValue(label, mu)){
    cout << label << " is not set" << endl;
    exit(0);
  }

  label = base_label + "/OMPnumThreads";
  if ( !tp.getInspectedValue(label, numOfOMP)){
    cout << label << " is not set" << endl;
    exit(0);
  }

  bd_u.resize(numOfNodeGlobal,vector<double>(3,0));
  bd_iu.resize(numOfNodeGlobal,vector<int>(3,0));

  bd_p.resize(numOfNodeGlobal,0);
  bd_ip.resize(numOfNodeGlobal,0);

  for(int i=0;i<numOfNodeGlobal;i++){
    for(int j=0;j<3;j++){
      bd_iu[i][j] = 1;
      bd_u[i][j] = 0;
    }
  }

  for(int i=0;i<numOfNodeGlobal;i++){
    bd_ip[i] = 1;
    bd_p[i] = 0;
  }

  string bdType;
  base_label = "/Boundary";

  label = base_label + "/bottom/type";
  if ( !tp.getInspectedValue(label,bdType)){
    cout << label << " is not set" << endl;
    exit(0);
  }

  if(bdType=="v"){
    double value[3];
    label = base_label + "/bottom/value";
    if ( !tp.getInspectedVector(label,value,3)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = 0;
    for(int i=0;i<nz+1;i++){
      for(int j=0;j<nx+1;j++){
        for(int k=0;k<3;k++){
          bd_iu[i*(nx+1)*(ny+1)+tmp+j][k] = 0;
          bd_u[i*(nx+1)*(ny+1)+tmp+j][k] = value[k];
        }
      }
    }
  }else if(bdType=="p"){
    double value;
    label = base_label + "/bottom/value";
    if ( !tp.getInspectedValue(label,value)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = 0;
    for(int i=0;i<nz+1;i++){
      for(int j=0; j<nx+1; j++){
        bd_ip[i*(nx+1)*(ny+1)+tmp+j] = 0;
        bd_p[i*(nx+1)*(ny+1)+tmp+j] = value;
      }
    }
  }else if(bdType=="free"){
  }else{
    cout << bdType << " is undefined. Exit..." << endl;
  }


    label = base_label + "/top/type";
  if ( !tp.getInspectedValue(label,bdType)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  //ofstream top["top.dat");

  if(bdType=="v"){
    double value[3];
    label = base_label + "/top/value";
    if ( !tp.getInspectedVector(label,value,3)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = (nx+1)*ny;
    for(int i=0;i<nz+1;i++){
      for(int j=0; j<nx+1; j++){
        for(int k=0;k<3;k++){
          bd_iu[i*(nx+1)*(ny+1)+tmp+j][k] = 0;
          bd_u[i*(nx+1)*(ny+1)+tmp+j][k]  = value[k];
        }
      }
    }
  }else if(bdType=="p"){
    double value;
    label = base_label + "/top/value";
    if ( !tp.getInspectedValue(label,value)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = (nx+1)*ny;
    for(int i=0;i<nz+1;i++){
      for(int j=0; j<nx+1; j++){
        bd_ip[i*(nx+1)*(ny+1)+tmp+j] = 0;
        bd_p[i*(nx+1)*(ny+1)+tmp+j] = value;
      }
    }
  }else if(bdType=="free"){
  }else{
    cout << bdType << " is undefined. Exit..." << endl;
  }

  label = base_label + "/left/type";
  if ( !tp.getInspectedValue(label,bdType)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  
  if(bdType=="v"){
    double value[3];
    label = base_label + "/left/value";
    if ( !tp.getInspectedVector(label,value,3)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp=0;
    for(int i=0; i<nz+1; i++){
      for(int j=0; j<ny+1; j++){
        for(int k=0; k<3; k++){
          bd_iu[(nx+1)*(ny+1)*i+(nx+1)*j][k] = 0;
          bd_u[(nx+1)*(ny+1)*i+(nx+1)*j][k] = value[k];
        }
      }
    }
  }else if(bdType=="p"){
    double value;
    label = base_label + "/left/value";
    if ( !tp.getInspectedValue(label,value)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = 0;
    for(int i=0;i<nz+1;i++){
      for(int j=0; j<ny+1; j++){
        bd_ip[(nx+1)*(ny+1)*i+(nx+1)*j] = 0;
        bd_p[(nx+1)*(ny+1)*i+(nx+1)*j] = value;
      }
    }
  }else if(bdType=="free"){
  }else{
    cout << bdType << " is undefined. Exit..." << endl;
  }

  label = base_label + "/right/type";
  if ( !tp.getInspectedValue(label,bdType)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  
  if(bdType=="v"){
    double value[3];
    label = base_label + "/right/value";
    if ( !tp.getInspectedVector(label,value,3)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = nx;
    for(int i=0;i<nz+1;i++){
      for(int j=0; j<ny+1; j++){
        for(int k=0;k<3;k++){
          bd_iu[(nx+1)*(ny+1)*i+(nx+1)*j+tmp][k] = 0;
          bd_u[(nx+1)*(ny+1)*i+(nx+1)*j+tmp][k] = value[k];
        }
      }
    }
  }else if(bdType=="v_poiseuille"){
    int tmp = nx;
    for(int i=0;i<ny+1;i++){
      for(int j=0;j<2;j++){
        bd_iu[(nx+1)*i+tmp][j] = 0;
        //valueはinputImageInfoで設定
      }
    }
  }else if(bdType=="p"){
    double value;
    label = base_label + "/right/value";
    if ( !tp.getInspectedValue(label,value)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = 0;
    for(int i=0; i<nz+1;i++){
      for(int j=0; j<ny+1; j++){
        bd_ip[(nx+1)*(ny+1)*i+(nx+1)*j+tmp] = 0;
        bd_p[(nx+1)*(ny+1)*i+(nx+1)*j+tmp] = value;
      }
    }
  }else if(bdType=="free"){
  }else{
    cout << bdType << " is undefined. Exit..." << endl;
  }



  
  label = base_label + "/front/type";
  if ( !tp.getInspectedValue(label,bdType)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  
  if(bdType=="v"){
    double value[3];
    label = base_label + "/front/value";
    if ( !tp.getInspectedVector(label,value,3)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp=0;
    for(int i=0; i<ny+1; i++){
      for(int j=0; j<nx+1; j++){
        for(int k=0; k<3; k++){
          bd_iu[j+(nx+1)*i][k] = 0;
          bd_u[j+(nx+1)*i][k] = value[k];
        }
      }
    }
  }else if(bdType=="p"){
    double value;
    label = base_label + "/front/value";
    if ( !tp.getInspectedValue(label,value)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = 0;
    for(int i=0;i<ny+1;i++){
      for(int j=0; j<nx+1; j++){
        bd_ip[j+(nx+1)*i] = 0;
        bd_p[j+(nx+1)*i] = value;
      }
    }
  }else if(bdType=="free"){
  }else{
    cout << bdType << " is undefined. Exit..." << endl;
  }

  
  label = base_label + "/back/type";
  if ( !tp.getInspectedValue(label,bdType)){
    cout << label << " is not set" << endl;
    exit(0);
  }
  
  if(bdType=="v"){
    double value[3];
    label = base_label + "/back/value";
    if ( !tp.getInspectedVector(label,value,3)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp=0;
    for(int i=0; i<ny+1; i++){
      for(int j=0; j<nx+1; j++){
        for(int k=0; k<3; k++){
          bd_iu[j+(nx+1)*i+(nx+1)*(ny+1)*nz][k] = 0;
          bd_u[j+(nx+1)*i+(nx+1)*(ny+1)*nz][k] = value[k];
        }
      }
    }
  }else if(bdType=="p"){
    double value;
    label = base_label + "/back/value";
    if ( !tp.getInspectedValue(label,value)){
      cout << label << " is not set" << endl;
      exit(0);
    }
    int tmp = 0;
    for(int i=0;i<ny+1;i++){
      for(int j=0; j<nx+1; j++){
        bd_ip[j+(nx+1)*i+(nx+1)*(ny+1)*nz] = 0;
        bd_p[j+(nx+1)*i+(nx+1)*(ny+1)*nz] = value;
      }
    }
  }else if(bdType=="free"){
  }else{
    cout << bdType << " is undefined. Exit..." << endl;
  }


  ////setImage////
  phi.resize(numOfElmGlobal);
  phiEX.resize(numOfElmGlobal);
  sdf.resize(numOfNodeGlobal);

  for(int i=0;i<numOfElmGlobal;i++){
    phi[i] = 1e0;
    phiEX[i] = 1e0;
  }

  string imageFile;
  label = "/Domain/image";
  if ( !tp.getInspectedValue(label,imageFile)){
    cout << label << " is not set" << endl;
    exit(0);
  }

  double value;
  FILE *fp; 
  if ((fp = fopen(imageFile.c_str(), "r")) == NULL)
  {
    cout << imageFile << " open error" << endl;
    exit(1);
  }
  for(int k=0; k<nz; k++){
    for(int j=0;j<ny;j++){
      for(int i=0;i<nx;i++){
        int tmp1,tmp2,tmp3;
        fscanf(fp,"%d %d %d %lf\n",&tmp1,&tmp2,&tmp3,&value);
        phi[tmp1+tmp2*nx+tmp3*nx*ny] = value;
      }
    }
  }
  fclose(fp);

  label = "/Domain/imageEX";
  if ( !tp.getInspectedValue(label,imageFile)){
    cout << label << " is not set" << endl;
    exit(0);
  }


  if ((fp = fopen(imageFile.c_str(), "r")) == NULL)
  {
    cout << imageFile << " open error" << endl;
    exit(1);
  }
  for(int k=0; k<nz; k++){
    for(int j=0;j<ny;j++){
      for(int i=0;i<nx;i++){
        int tmp1,tmp2,tmp3;
        fscanf(fp,"%d %d %d %lf\n",&tmp1,&tmp2,&tmp3,&value);
        phiEX[tmp1+tmp2*nx+tmp3*nx*ny] = value;
      }
    }
  }
  fclose(fp);

  label = "/Domain/sdf";
  if ( !tp.getInspectedValue(label,imageFile)){
    cout << label << " is not set" << endl;
    exit(0);
  }

  if ((fp = fopen(imageFile.c_str(), "r")) == NULL)
  {
    cout << imageFile << " open error" << endl;
    exit(1);
  }
  for(int k=0; k<nz+1; k++){
    for(int j=0;j<ny+1;j++){
      for(int i=0;i<nx+1;i++){
        int tmp1,tmp2,tmp3;
        fscanf(fp,"%d %d %d %lf\n",&tmp1,&tmp2,&tmp3,&value);
        //value = -value;
        if(fabs(value)<1e-10) {
          value = 0e0;
        }
        sdf[tmp1+tmp2*(nx+1)+tmp3*(nx+1)*(ny+1)] = value;
      }
    }
  }
  fclose(fp);
  


}