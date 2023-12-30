#include "FEM.h"

void FEM::line_serch(vector<int> &sub_cross_num, vector<vector<double>> &sub_x_tmp, vector<double> &sdf_current, vector<vector<double>> &x_current, const int ic)
{
  int numOfLineInElm = 12;
  double s,t,u;
  int tmp = 0;

  for(int i=0; i<numOfLineInElm; i++){

    bool cross = is_cross(sdf_current, i);
    if(!cross) continue;

    tmp++;
    vector<double> N(numOfNodeInElm);
    
    int linePartition = 500;
    double sdf_min = 1e10;
    double s_min_point;
    vector<double> cross_point(3);
    
    for(int j=0; j<linePartition; j++){
      s_t_u(linePartition,i,j,s,t,u);
      ShapeFunction3D::C3D8_N(N,s,t,u);
    
      double sdf_line = 0e0;
      for(int p=0;p<numOfNodeInElm;p++){
        sdf_line += N[p] * sdf_current[p];
      } 
      if(fabs(sdf_line) < fabs(sdf_min)){
        sdf_min = sdf_line;
        s_min_point = s;
      }
    }
    //cout << sdf_min << endl;
    
    ShapeFunction3D::C3D8_N(N,s_min_point,t,u);
    for(int p=0;p<3;p++){
      for(int q=0; q<numOfNodeInElm; q++){
        cross_point[p] += N[q] * x_current[q][p];
      }
    }
    sub_x_tmp.push_back(cross_point);
  }
  sub_cross_num[ic] = tmp;
}


bool FEM::is_cross(vector<double> &sdf_current, const int i)
{
  if(i==0){
    if(sdf_current[0]*sdf_current[1] > 0) return false;
  }else if(i==1){
    if(sdf_current[1]*sdf_current[2] > 0) return false;
  }else if(i==2){
    if(sdf_current[2]*sdf_current[3] > 0) return false;
  }else if(i==3){
    if(sdf_current[3]*sdf_current[0] > 0) return false;
  }else if(i==4){
    if(sdf_current[4]*sdf_current[5] > 0) return false;
  }else if(i==5){
    if(sdf_current[5]*sdf_current[6] > 0) return false;
  }else if(i==6){
    if(sdf_current[6]*sdf_current[7] > 0) return false;
  }else if(i==7){
    if(sdf_current[7]*sdf_current[4] > 0) return false;
  }else if(i==8){
    if(sdf_current[0]*sdf_current[4] > 0) return false;
  }else if(i==9){
    if(sdf_current[1]*sdf_current[5] > 0) return false;
  }else if(i==10){
    if(sdf_current[2]*sdf_current[6] > 0) return false;
  }else if(i==11){
    if(sdf_current[3]*sdf_current[7] > 0) return false;
  }else{
    cout << "error" << endl;
  }
  return true;
}

void FEM::s_t_u(const int &linePartition, const int &i, const int &j, double &s, double &t, double &u)
{
  if(i==0){
    s = -1e0 + j*(2e0/linePartition); t = -1e0; u = -1e0;
  }else if(i==1){
    s = 1e0; t = -1e0 + j*(2e0/linePartition); u = -1e0;
  }else if(i==2){
    s = 1e0 - j*(2e0/linePartition); t = 1e0; u = -1e0;
  }else if(i==3){
    s = -1e0; t = 1e0 - j*(2e0/linePartition); u = -1e0;
  }else if(i==4){
    s = -1e0 + j*(2e0/linePartition); t = -1e0; u = 1e0;
  }else if(i==5){
    s = 1e0; t = -1e0 + j*(2e0/linePartition); u = 1e0;
  }else if(i==6){
    s = 1e0 - j*(2e0/linePartition); t = 1e0; u = 1e0;
  }else if(i==7){
    s = -1e0; t = 1e0 - j*(2e0/linePartition); u = 1e0;
  }else if(i==8){
    s = -1e0; t = -1e0; u = -1e0 + j*(2e0/linePartition);
  }else if(i==9){
    s = 1e0; t = -1e0; u = -1e0 + j*(2e0/linePartition);
  }else if(i==10){
    s = 1e0; t = 1e0; u = -1e0 + j*(2e0/linePartition);
  }else if(i==11){
    s = -1e0; t = 1e0; u = -1e0 + j*(2e0/linePartition);
  }else{
    cout << "error" << endl;
  }
}