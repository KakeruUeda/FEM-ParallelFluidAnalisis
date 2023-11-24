#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <omp.h>

class FEM :public{
    public:
    TextParser tp;
    int OMPnumThreads;
    std::string outputDir,fileName;
    std::vector<ARRAY2D<double>> LHS_Diff,LHS_Kp;
    std::vector<ARRAY3D<double>> LHS_P;
    std::vector<ARRAY4D<double>> LHS_Adv,LHS_SUPG_u;
    std::vector<ARRAY3D<double>> LHS_SUPG_p;

    std::vector<ARRAY4D<double>> LHS_uu;
    std::vector<ARRAY3D<double>> LHS_up;
    std::vector<ARRAY3D<double>> LHS_pu;
    std::vector<ARRAY2D<double>> LHS_pp;

    std::vector<ARRAY2D<double>> Qu;
    std::vector<ARRAY1D<double>> Qp;

  double rho,mu,nu;
  double NRtolerance;
  double relaxationParam;
}