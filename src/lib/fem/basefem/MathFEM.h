#include <define.h>
#include <iostream>
#include <vector>
using namespace std;

class MathFEM{
  public:
    static void calc_dxdr(double (&dxdr)[3][3],VDOUBLE2D &dNdr,VDOUBLE2D &x1,const int &numOfNodeInElm);
    static void calc_dNdx(VDOUBLE2D &dNdx,VDOUBLE2D &dNdr,const double (&dxdr)[3][3],const int &numOfNodeInElm);
};

