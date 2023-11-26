#include <iostream>
#include <vector>
using namespace std;

class ShapeFunction3D{
 public:

  static void C3D8_N(vector<double> &N,const double &g1,const double &g2,const double &g3)
  {
  N[0]= 1.25e-1 * (1e0-g1) * (1e0-g2) * (1e0-g3);
  N[1] = 1.25e-1 * (1e0+g1) * (1e0-g2) * (1e0-g3);
  N[2] = 1.25e-1 * (1e0+g1) * (1e0+g2) * (1e0-g3);
  N[3] = 1.25e-1 * (1e0-g1) * (1e0+g2) * (1e0-g3);
  N[4] = 1.25e-1 * (1e0-g1) * (1e0-g2) * (1e0+g3);
  N[5] = 1.25e-1 * (1e0+g1) * (1e0-g2) * (1e0+g3);
  N[6] = 1.25e-1 * (1e0+g1) * (1e0+g2) * (1e0+g3);
  N[7] = 1.25e-1 * (1e0-g1) * (1e0+g2) * (1e0+g3);
  }


  static void C3D8_dNdr(vector<vector<double>> &dNdr,const double &g1,const double &g2,const double &g3)
  {
    dNdr[0][0] = -1.25e-1 * (1e0-g2) * (1e0-g3);
    dNdr[0][1] = -1.25e-1 * (1e0-g1) * (1e0-g3);
    dNdr[0][2] = -1.25e-1 * (1e0-g1) * (1e0-g2);
    dNdr[1][0] =  1.25e-1 * (1e0-g2) * (1e0-g3);
    dNdr[1][1] = -1.25e-1 * (1e0+g1) * (1e0-g3);
    dNdr[1][2] = -1.25e-1 * (1e0+g1) * (1e0-g2);
    dNdr[2][0] =  1.25e-1 * (1e0+g2) * (1e0-g3);
    dNdr[2][1] =  1.25e-1 * (1e0+g1) * (1e0-g3);
    dNdr[2][2] = -1.25e-1 * (1e0+g1) * (1e0+g2);
    dNdr[3][0] = -1.25e-1 * (1e0+g2) * (1e0-g3);
    dNdr[3][1] =  1.25e-1 * (1e0-g1) * (1e0-g3);
    dNdr[3][2] = -1.25e-1 * (1e0-g1) * (1e0+g2);
    dNdr[4][0] = -1.25e-1 * (1e0-g2) * (1e0+g3);
    dNdr[4][1] = -1.25e-1 * (1e0-g1) * (1e0+g3);
    dNdr[4][2] =  1.25e-1 * (1e0-g1) * (1e0-g2);
    dNdr[5][0] =  1.25e-1 * (1e0-g2) * (1e0+g3);
    dNdr[5][1] = -1.25e-1 * (1e0+g1) * (1e0+g3);
    dNdr[5][2] =  1.25e-1 * (1e0+g1) * (1e0-g2);
    dNdr[6][0] =  1.25e-1 * (1e0+g2) * (1e0+g3);
    dNdr[6][1] =  1.25e-1 * (1e0+g1) * (1e0+g3);
    dNdr[6][2] =  1.25e-1 * (1e0+g1) * (1e0+g2);
    dNdr[7][0] = -1.25e-1 * (1e0+g2) * (1e0+g3);
    dNdr[7][1] =  1.25e-1 * (1e0-g1) * (1e0+g3);
    dNdr[7][2] =  1.25e-1 * (1e0-g1) * (1e0+g2);
  }
 
};