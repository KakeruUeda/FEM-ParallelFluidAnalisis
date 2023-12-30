#include "FEM.h"

void FEM::initialize()
{ 
  readInput(); 
  setDomain();
  setBoundary();
  setFluidDomain();
  prepareMatrix();

  //binaryTreeSubDivision();
  //interfacePartition();
}


void FEM::setDomain()
{
  x.resize(numOfNodeGlobal,vector<double>(3, 0));
  element.resize(numOfElmGlobal);

  for(int ic=0;ic<numOfElmGlobal;ic++){
    element[ic].resize(8);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  int tmp2=0;
  for(int k=0; k<nz+1; k++){
    for(int i=0; i<ny+1; i++){
      for(int j=0; j<nx+1; j++){
        x[tmp2][0]=j*dx;
        x[tmp2][1]=i*dy;
        x[tmp2][2]=k*dz; 
        tmp2++;
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  tmp2=0;
  for(int k=0;k<nz;k++){
    for(int i=0;i<ny;i++){
      for(int j=0;j<nx;j++){
        element[tmp2][0]=j   +i*(nx+1) +k*(nx+1)*(ny+1);
        element[tmp2][1]=j+1 +i*(nx+1) +k*(nx+1)*(ny+1);
        element[tmp2][2]=j+1 +(i+1)*(nx+1) +k*(nx+1)*(ny+1);
        element[tmp2][3]=j   +(i+1)*(nx+1) +k*(nx+1)*(ny+1);
        element[tmp2][4]=j   +i*(nx+1) +(k+1)*(nx+1)*(ny+1);
        element[tmp2][5]=j+1 +i*(nx+1) +(k+1)*(nx+1)*(ny+1);
        element[tmp2][6]=j+1 +(i+1)*(nx+1) +(k+1)*(nx+1)*(ny+1);
        element[tmp2][7]=j   +(i+1)*(nx+1) +(k+1)*(nx+1)*(ny+1);
        tmp2++;
      }
    }
  }
  prepare();
}

void FEM::prepare()
{
  /*
  elm = new ElementBaseFEM* [numOfElmGlobal];

  for(int ee=0;ee<numOfElmGlobal;++ee){
    elm[ee] = new ElementBaseFEM();
    elm[ee]->SolnData = &(SolnData);
  }

  if(solverPetsc != NULL) delete solverPetsc;
  solverPetsc = NULL;

  solverPetsc = (PetscSolver*) new PetscSolver;

  SolnData.initialise(numOfNodeGlobal*numOfDofsNode);
  */

}


void FEM::setBoundary(){
  
  for(int ic=0;ic<numOfElmGlobal;ic++){
    
    if(phi[ic]>1e-10) continue;

    for(int i=0;i<numOfNodeInElm;i++){
      for(int j=0;j<3;j++){
        bd_iu[element[ic][i]][j] = 0;
        bd_u[element[ic][i]][j] = 0;
        //bd_ip[element[ic][i]] = 0;
        //bd_p[element[ic][i]] = 0;
      }
    }
  }
  
  /*
  int count = 0;
  vector<double> vecDbTmp(3,0);

  for(int i=0;i<numOfNodeGlobal;i++){
    for(int k=0; k<3; k++){
      if(bd_iu[i][k] == 0){
        vecDbTmp[0] = i;
        vecDbTmp[1] = k;
        vecDbTmp[2] = bd_u[i][k];
        DirichletBCs_tmp.push_back(vecDbTmp);
        count++;
      }
      if(bd_ip[i] == 0){
        vecDbTmp[0] = i;
        vecDbTmp[1] = 3;
        vecDbTmp[2] = bd_p[i];
        DirichletBCs_tmp.push_back(vecDbTmp);
        count++;
      }
    }
  }

  ///pressure///

  numOfBdNode = count;
  */

}


void FEM::setFluidDomain()
{
  numOfElmGlobalFluid=0;
  
  for(int ic=0; ic<numOfElmGlobal; ic++){
    if(phi[ic]<1e-12){
      continue;
    }
    sortElm.push_back(ic);
    elementFluidPrev.push_back(element[ic]);
    numOfElmGlobalFluid++;
  }


  sortNode.resize(numOfElmGlobalFluid*numOfNodeInElm);
  
  int tmp2 = 0;
  for(int ic=0; ic<numOfElmGlobalFluid; ic++){
    for(int ii=0; ii<numOfNodeInElm; ii++){
      sortNode[tmp2] = elementFluidPrev[ic][ii];
      tmp2++;
    }
  }
  sort(sortNode.begin(), sortNode.end());
  sortNode.erase(unique(sortNode.begin(), sortNode.end()), sortNode.end());
  
  vector<int> sortNodeNew(numOfNodeGlobal,0); 
  

  ///// sortNodeNew: start with 0 /////
  numOfNodeGlobalFluid = 0;
  for(auto it = sortNode.begin(); it != sortNode.end(); it++){
    sortNodeNew[*it] = numOfNodeGlobalFluid;
    numOfNodeGlobalFluid++;
  }

  int n1, n2;
  elementFluid.resize(numOfElmGlobalFluid,vector<int>(numOfNodeInElm,0));
  for(int ii=0; ii<numOfElmGlobalFluid; ii++){
    for(int jj=0; jj<numOfNodeInElm; jj++){
      n1 = elementFluidPrev[ii][jj];
      elementFluid[ii][jj] = sortNodeNew[n1];
    }
  }

  /*
  if(myId == 0){
    for(int i=0; i<numOfElmFluidGlobal; i++){
      for(int j=0; j<numOfNodeInElm; j++){
        //cout << "i = " << i << " j= " << j << " elementFluidPrev node = " << elementFluidPrev[i][j] << endl;
      }
    }
    for(int i=0; i<numOfElmFluidGlobal; i++){
      for(int j=0; j<numOfNodeInElm; j++){
        //cout << "i = " << i << " j= " << j << " elementFluid node = " << elementFluid[i][j] << endl;
      }
    }
    //exit(1);
  }
  */

  bd_iu_fluid.resize(numOfNodeGlobalFluid,vector<int>(3,0));
  bd_u_fluid.resize(numOfNodeGlobalFluid,vector<double>(3,0e0));
  bd_ip_fluid.resize(numOfNodeGlobalFluid,0);
  bd_p_fluid.resize(numOfNodeGlobalFluid,0e0);

  phiFluid.resize(numOfElmGlobalFluid,0e0);
  phiEXFluid.resize(numOfElmGlobalFluid,0e0);
  phiVOFFluid.resize(numOfElmGlobalFluid,0e0);
  sdfFluid.resize(numOfElmGlobalFluid,0e0);
  sdfFluid_node.resize(numOfNodeGlobalFluid,0e0);
  
  for(int ii=0; ii<numOfNodeGlobalFluid; ii++){
    for(int kk=0; kk<3; kk++){
      bd_iu_fluid[ii][kk] = bd_iu[sortNode[ii]][kk];
      bd_u_fluid[ii][kk] = bd_u[sortNode[ii]][kk];
    }
    sdfFluid_node[ii] = sdf_node[sortNode[ii]];
    bd_ip_fluid[ii] = bd_ip[sortNode[ii]];
    bd_p_fluid[ii] = bd_p[sortNode[ii]];
  }

  for(int ii=0; ii<numOfElmGlobalFluid; ii++){
    sdfFluid[ii] = sdf[sortElm[ii]];
    phiFluid[ii] = phi[sortElm[ii]];
    phiEXFluid[ii] = phiEX[sortElm[ii]];
    phiVOFFluid[ii] = phiVOF[sortElm[ii]];
  }

  xFluid.resize(numOfNodeGlobalFluid,vector<double>(3, 0));
  for(int ii=0; ii<numOfNodeGlobalFluid; ii++){
    for(int kk=0; kk<3; kk++){
      xFluid[ii][kk] = x[sortNode[ii]][kk];
    }
  }

  /*
  if(myId == 0){
    for(int i=0; i<numOfNodeGlobalFluid; i++){
      cout << "i = " << i << " " ;
      for(int j=0; j<3; j++){
        cout << bd_iu_fluid[i][j] << " " ;
      }
      cout << endl;
    }
    for(int i=0; i<numOfNodeGlobalFluid; i++){
      cout << "i = " << i << " " ;
      for(int j=0; j<3; j++){
        cout << bd_u_fluid[i][j] << " " ;
      }
      cout << endl;
    }
    //exit(1);
  }
  */

  int count = 0;
  vector<double> vecDbTmp(3,0);

  for(int ii=0;ii<numOfNodeGlobalFluid;ii++){
    for(int kk=0; kk<3; kk++){
      if(bd_iu_fluid[ii][kk] == 0){
        vecDbTmp[0] = ii;
        vecDbTmp[1] = kk;
        vecDbTmp[2] = bd_u_fluid[ii][kk];
        DirichletBCsFluid_tmp.push_back(vecDbTmp);
        count++;
      }
    }
  }
  for(int ii=0;ii<numOfNodeGlobalFluid;ii++){
    if(bd_ip_fluid[ii] == 0){
      vecDbTmp[0] = ii;
      vecDbTmp[1] = 3;
      vecDbTmp[2] = bd_p_fluid[ii];
      DirichletBCsFluid_tmp.push_back(vecDbTmp);
      count++;
    }
  }

  numOfBdNodeFluid = count;


  //// FLUID ONLY PREPARE ////
  elmFluid = new ElementBaseFEM* [numOfElmGlobalFluid];

  for(int ee=0;ee<numOfElmGlobalFluid;ee++){
    elmFluid[ee] = new ElementBaseFEM();
    elmFluid[ee]->SolnDataFluid = &(SolnDataFluid);
  }

  if(solverPetscFluid != NULL) delete solverPetscFluid;
  solverPetscFluid = NULL;

  solverPetscFluid = (PetscSolver*) new PetscSolver;

  SolnDataFluid.initialise(numOfNodeGlobalFluid*numOfDofsNode);

}


void FEM::binaryTreeSubDivision()
{
  sub_div_total = sub_div * sub_div * sub_div;

  for(int ic=0; ic<numOfElmGlobalFluid; ic++)
  { 
    if(phiEXFluid[ic] > 0.999) continue;
    //cout << "numOfElmGlobalFluid = " << numOfElmGlobalFluid << " ic = " << ic << endl;
    if(elmFluid[ic]->getSubdomainId() == myId)
    {
      vector<vector<double>> x_parent;
      vector<double> sdf_parent;
      
      x_parent.resize(numOfNodeInElm, vector<double>(3, 0e0));
      sdf_parent.resize(numOfNodeInElm);
      
      for(int i=0; i<numOfNodeInElm; i++){
        sdf_parent[i] = sdfFluid_node[elmFluid[ic]->nodeNumsPrevFluid[i]];
        //if(sdf_parent[i] == 0) cout << sdf_parent[i] << endl;
        for(int j=0;j<3;j++){
          x_parent[i][j] = xFluid[elmFluid[ic]->nodeNumsPrevFluid[i]][j];
        }
      }
      int depth = 0;
      int tmp = 0;
      vector<double> x_center_parent(3,0e0);
      x_center_parent[0] = (x_parent[0][0] + x_parent[1][0]) / 2e0;
      x_center_parent[1] = (x_parent[0][1] + x_parent[3][1]) / 2e0;
      x_center_parent[2] = (x_parent[0][2] + x_parent[4][2]) / 2e0;
      
      gererateSubElms(sdf_parent, x_parent, x_center_parent, ic, tmp, depth);
    }
  }
}

void FEM::gererateSubElms(vector<double> &sdf_parent, vector<vector<double>> &x_sub, vector<double> &x_center_parent, const int &ic, int &tmp, int &depth)
{
  depth = depth + 1;

  for(int kk=0; kk<sub_div; kk++){
    for(int jj=0; jj<sub_div; jj++){
      for(int ii=0; ii<sub_div; ii++){
        
        vector<vector<double>> x_sub_sub;
        vector<double> sdf_sub_sub;
        
        x_sub_sub.resize(numOfNodeInElm, vector<double>(3, 0e0));
        sdf_sub_sub.resize(numOfNodeInElm, 0e0);
        getSubSubCoordinates(x_sub, x_sub_sub, sdf_sub_sub, ii, jj, kk);
        
        vector<double> N(numOfNodeInElm);
        
        int check = 0;
        
        for(int i=0; i<numOfNodeInElm; i++)
        {
          double point_x = (x_sub_sub[i][0] - x_center_parent[0]) / (dx/2e0);
          double point_y = (x_sub_sub[i][1] - x_center_parent[1]) / (dy/2e0);
          double point_z = (x_sub_sub[i][2] - x_center_parent[2]) / (dz/2e0);
          
          double esp = 1e-8;
          if(point_x > 1e0+esp || point_x < -1e0-esp){ 
            cout << "point_x interpolation error..." << point_x ; 
            exit(1); 
          }
          if(point_y > 1e0+esp || point_y < -1e0-esp){ 
            cout << "point_y interpolation error..." << point_y ; 
            exit(1);
          }
          if(point_z > 1e0+esp || point_z < -1e0-esp){ 
            cout << "point_z interpolation error..." << point_z ; 
            exit(1);
          }
 
          ShapeFunction3D::C3D8_N(N, point_x, point_y, point_z);
          
          sdf_sub_sub[i] = 0e0;
          for(int j=0; j<numOfNodeInElm; j++){
            sdf_sub_sub[i] += N[j] * sdf_parent[j];
          }

          if(sdf_sub_sub[i] < 0e0) check = check + 1;
        }

        if(depth < max_depth){
          if(check != 0 && check != numOfNodeInElm){
            gererateSubElms(sdf_parent, x_sub_sub, x_center_parent, ic, tmp, depth);
            if(ii == sub_div-1 && jj == sub_div-1 && kk == sub_div-1){
              depth = depth - 1;
            }        
          }else{
            makeSubElmsData(sdf_sub_sub, x_sub_sub, x_center_parent, ic, tmp, depth);
            if(ii == sub_div-1 && jj == sub_div-1 && kk == sub_div-1){
              depth = depth - 1;
            }        
          }
        }else{
          makeSubElmsData(sdf_sub_sub, x_sub_sub, x_center_parent, ic, tmp, depth);
          if(ii == sub_div-1 && jj == sub_div-1 && kk == sub_div-1){
            depth = depth - 1;
          }        
        }

        //vector<double>().swap(sdf_sub_sub);
        //vector<vector<double>>().swap(x_sub_sub);
        //vector<double>().swap(N);
      }
    }
  }

}


void FEM::getSubSubCoordinates(vector<vector<double>> &x_sub, vector<vector<double>> &x_sub_sub, vector<double> &sdf_sub_sub, const int &ii, const int &jj, const int &kk)
{
  
  double x_orig = x_sub[0][0];
  double y_orig = x_sub[0][1];
  double z_orig = x_sub[0][2];

  double dx_sub_sub = (x_sub[1][0] - x_sub[0][0]) / (double)sub_div;
  double dy_sub_sub = (x_sub[3][1] - x_sub[0][1]) / (double)sub_div;
  double dz_sub_sub = (x_sub[4][2] - x_sub[0][2]) / (double)sub_div;

  x_sub_sub[0][0] = x_orig + ii * dx_sub_sub;
  x_sub_sub[0][1] = y_orig + jj * dy_sub_sub;
  x_sub_sub[0][2] = z_orig + kk * dz_sub_sub;

  x_sub_sub[1][0] = x_orig + ii * dx_sub_sub + dx_sub_sub;
  x_sub_sub[1][1] = y_orig + jj * dy_sub_sub;
  x_sub_sub[1][2] = z_orig + kk * dz_sub_sub;
  
  x_sub_sub[2][0] = x_orig + ii * dx_sub_sub + dx_sub_sub;
  x_sub_sub[2][1] = y_orig + jj * dy_sub_sub + dy_sub_sub;
  x_sub_sub[2][2] = z_orig + kk * dz_sub_sub;

  x_sub_sub[3][0] = x_orig + ii * dx_sub_sub;
  x_sub_sub[3][1] = y_orig + jj * dy_sub_sub + dy_sub_sub;
  x_sub_sub[3][2] = z_orig + kk * dz_sub_sub;

  x_sub_sub[4][0] = x_orig + ii * dx_sub_sub;
  x_sub_sub[4][1] = y_orig + jj * dy_sub_sub;
  x_sub_sub[4][2] = z_orig + kk * dz_sub_sub + dz_sub_sub;

  x_sub_sub[5][0] = x_orig + ii * dx_sub_sub + dx_sub_sub;
  x_sub_sub[5][1] = y_orig + jj * dy_sub_sub;
  x_sub_sub[5][2] = z_orig + kk * dz_sub_sub + dz_sub_sub;
  
  x_sub_sub[6][0] = x_orig + ii * dx_sub_sub + dx_sub_sub;
  x_sub_sub[6][1] = y_orig + jj * dy_sub_sub + dy_sub_sub;
  x_sub_sub[6][2] = z_orig + kk * dz_sub_sub + dz_sub_sub;

  x_sub_sub[7][0] = x_orig + ii * dx_sub_sub;
  x_sub_sub[7][1] = y_orig + jj * dy_sub_sub + dy_sub_sub;
  x_sub_sub[7][2] = z_orig + kk * dz_sub_sub + dz_sub_sub;

}


void FEM::makeSubElmsData(vector<double> &sdf_sub_sub, vector<vector<double>> &x_sub_sub, vector<double> &x_center_parent, const int &ic, int &tmp, int &depth)
{
  elmFluid[ic]->sub_elm_node.emplace_back();
  elmFluid[ic]->sub_elm.emplace_back();

  for(int i=0; i<numOfNodeInElm; i++)
  {
    double sdf_tmp = sdf_sub_sub[i];  
    double x_tmp = x_sub_sub[i][0];  
    double y_tmp = x_sub_sub[i][1];  
    double z_tmp = x_sub_sub[i][2]; 
    elmFluid[ic]->sub_elm_node[tmp].emplace_back(); 
    elmFluid[ic]->sub_elm_node[tmp][i].sub_elm_sdf.push_back(sdf_tmp);
    elmFluid[ic]->sub_elm_node[tmp][i].sub_elm_x.push_back(x_tmp);
    elmFluid[ic]->sub_elm_node[tmp][i].sub_elm_x.push_back(y_tmp);
    elmFluid[ic]->sub_elm_node[tmp][i].sub_elm_x.push_back(z_tmp);
  }
  
  int GP = 2;

  for(int kk=0; kk<GP; kk++){
    for(int jj=0; jj<GP; jj++){
      for(int ii=0; ii<GP; ii++){
        
        double sub_gx_tmp, sub_gy_tmp, sub_gz_tmp;
        double weight;
        
        subGaussPoint(x_sub_sub, sub_gx_tmp, sub_gy_tmp, sub_gz_tmp, x_center_parent, ii, jj, kk);
        subGaussWeight(weight, ii, jj, kk, depth);
        
        elmFluid[ic]->sub_elm[tmp].sub_gx.push_back(sub_gx_tmp);
        elmFluid[ic]->sub_elm[tmp].sub_gy.push_back(sub_gy_tmp);
        elmFluid[ic]->sub_elm[tmp].sub_gz.push_back(sub_gz_tmp);
       
        elmFluid[ic]->sub_elm[tmp].sub_weight.push_back(weight);
      }
    }
  }

  tmp++;
}

void FEM::subGaussPoint(vector<vector<double>> &x_sub_sub, double &sub_gx_tmp, double &sub_gy_tmp, double &sub_gz_tmp, vector<double> &x_center_parent, const int &ii, const int &jj, const int &kk)
{
  Gauss gauss(2);

  double x_parent_center = x_center_parent[0];
  double y_parent_center = x_center_parent[1];
  double z_parent_center = x_center_parent[2];

  double x_sub_center = (x_sub_sub[1][0] + x_sub_sub[0][0]) / 2e0;
  double y_sub_center = (x_sub_sub[3][1] + x_sub_sub[0][1]) / 2e0;
  double z_sub_center = (x_sub_sub[4][2] + x_sub_sub[0][2]) / 2e0;

  double dx_sub = (x_sub_sub[1][0] - x_sub_sub[0][0]);
  double dy_sub = (x_sub_sub[3][1] - x_sub_sub[0][1]);
  double dz_sub = (x_sub_sub[4][2] - x_sub_sub[0][2]);
  
  double dx_sub_half = (x_sub_sub[1][0] - x_sub_sub[0][0]) / 2e0;
  double dy_sub_half = (x_sub_sub[3][1] - x_sub_sub[0][1]) / 2e0;
  double dz_sub_half = (x_sub_sub[4][2] - x_sub_sub[0][2]) / 2e0;

  double x_length_from_center = gauss.point[ii] * dx_sub_half;
  double y_length_from_center = gauss.point[jj] * dy_sub_half;
  double z_length_from_center = gauss.point[kk] * dz_sub_half;
  
  double gauss_gx_real = x_sub_center + x_length_from_center;
  double gauss_gy_real = y_sub_center + y_length_from_center;
  double gauss_gz_real = z_sub_center + z_length_from_center;

  sub_gx_tmp = (gauss_gx_real - x_parent_center) / (dx/2e0);
  sub_gy_tmp = (gauss_gy_real - y_parent_center) / (dy/2e0);
  sub_gz_tmp = (gauss_gz_real - z_parent_center) / (dz/2e0);

  if(sub_gx_tmp > 1 || sub_gx_tmp < -1){ 
    cout << "x interpolation error..." ; 
    exit(1); 
  }
  if(sub_gy_tmp > 1 || sub_gy_tmp < -1){ 
    cout << "y interpolation error..." ; 
    exit(1);
  }
  if(sub_gz_tmp > 1 || sub_gz_tmp < -1){ 
    cout << "z interpolation error..." ; 
    exit(1);
  }

}

void FEM::subGaussWeight(double &weight, const int &ii, const int &jj, const int &kk, int &depth)
{
  int division = sub_div * pow(sub_div, depth-1);
  Gauss gauss(2);
  
  weight = (gauss.weight[ii] * gauss.weight[jj] * gauss.weight[kk]) / (double)(division * division * division);
}