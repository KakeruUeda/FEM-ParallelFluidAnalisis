#include "mathFEM.h"
#include "basicFunctions.h"

void MathFEM::calc_dxdr(double (&dxdr)[3][3],vector<vector<double>>  &dNdr,vector<vector<double>>  &x1,const int &numOfNodeInElm){
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      dxdr[i][j] = 0e0;
      for(int p=0;p<numOfNodeInElm;p++){
        dxdr[i][j] += dNdr[p][j] * x1[p][i];
      }
    }
  }
}
void MathFEM::calc_dNdx(vector<vector<double>> &dNdx,vector<vector<double>> &dNdr,const double (&dxdr)[3][3],const int &numOfNodeInElm){
  double drdx[3][3];

  BasicFunctions::calcInverseMatrix_3x3(drdx,dxdr);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      dNdx[p][i] = 0e0;
      for(int j=0;j<3;j++) dNdx[p][i] += dNdr[p][j] * drdx[j][i];
    }
  }
}