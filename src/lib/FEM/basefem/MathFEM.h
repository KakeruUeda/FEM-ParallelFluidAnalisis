#include <iostream>
#include <vector>
using namespace std;

class MathFEM{
  public:
    static void calc_dxdr(double (&dxdr)[3][3],vector<vector<double>> &dNdr,vector<vector<double>> &x1,const int &numOfNodeInElm);
    static void calc_dNdx(vector<vector<double>> &dNdx,vector<vector<double>> &dNdr,const double (&dxdr)[3][3],const int &numOfNodeInElm);
};

