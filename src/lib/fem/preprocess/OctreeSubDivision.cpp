#include "FEM.h"

void FEM::octreeSubDivision()
{
  sub_div_total = sub_div * sub_div * sub_div;

  for(int ic=0; ic<numOfElmGlobalFluid; ic++)
  { 
    if(phiEXFluid[ic] > 0.999) continue;

    if(elmFluid[ic]->getSubdomainId() == myId)
    {
      VDOUBLE2D x_parent;
      VDOUBLE1D sdf_parent;
      
      x_parent.resize(numOfNodeInElm, VDOUBLE1D(3, 0e0));
      sdf_parent.resize(numOfNodeInElm);
      
      for(int i=0; i<numOfNodeInElm; i++)
      {
        sdf_parent[i] = sdfFluid[elmFluid[ic]->nodeNumsPrevFluid[i]];

        for(int j=0;j<3;j++)
        {
          x_parent[i][j] = xFluid[elmFluid[ic]->nodeNumsPrevFluid[i]][j];
        }
      }
      int depth = 0;
      int tmp = 0;
      VDOUBLE1D x_center_parent(3,0e0);
      x_center_parent[0] = (x_parent[0][0] + x_parent[1][0]) / 2e0;
      x_center_parent[1] = (x_parent[0][1] + x_parent[3][1]) / 2e0;
      x_center_parent[2] = (x_parent[0][2] + x_parent[4][2]) / 2e0;
      
      gererateSubElms(sdf_parent, x_parent, x_center_parent, ic, tmp, depth);
    }
  }
}


void FEM::gererateSubElms(VDOUBLE1D &sdf_parent, VDOUBLE2D &x_sub, VDOUBLE1D &x_center_parent, const int &ic, int &tmp, int &depth)
{
  depth = depth + 1;

  for(int kk=0; kk<sub_div; kk++){
    for(int jj=0; jj<sub_div; jj++){
      for(int ii=0; ii<sub_div; ii++){
        
        VDOUBLE2D x_sub_sub;
        VDOUBLE1D sdf_sub_sub;
        
        x_sub_sub.resize(numOfNodeInElm, VDOUBLE1D(3, 0e0));
        sdf_sub_sub.resize(numOfNodeInElm, 0e0);
        getSubSubCoordinates(x_sub, x_sub_sub, sdf_sub_sub, ii, jj, kk);
        
        VDOUBLE1D N(numOfNodeInElm);
        
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
        
        //VDOUBLE1D().swap(sdf_sub_sub);
        //VDOUBLE2D().swap(x_sub_sub);
        //VDOUBLE1D().swap(N);
      }
    }
  }

}


void FEM::getSubSubCoordinates(VDOUBLE2D &x_sub, VDOUBLE2D &x_sub_sub, VDOUBLE1D &sdf_sub_sub, const int &ii, const int &jj, const int &kk)
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


void FEM::makeSubElmsData(VDOUBLE1D &sdf_sub_sub, VDOUBLE2D &x_sub_sub, VDOUBLE1D &x_center_parent, const int &ic, int &tmp, int &depth)
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

void FEM::subGaussPoint(VDOUBLE2D &x_sub_sub, double &sub_gx_tmp, double &sub_gy_tmp, double &sub_gz_tmp, VDOUBLE1D &x_center_parent, const int &ii, const int &jj, const int &kk)
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