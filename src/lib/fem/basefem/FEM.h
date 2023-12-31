#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <omp.h>

#include "gauss.h"
#include "metis.h"
#include "TextParser.h"
#include "DomainFEM.h"
#include "PetscSolver.h"
#include "BasicFunctions.h"
#include "ShapeFunction.h"
#include "ElementBaseFEM.h"
#include "MathFEM.h"
using namespace std;


enum class SOLVER{
  STEADY_STOKES = 0,
  STEADY_NAVIERSTOKES = 1,
  UNSTEADY_NAVIERSTOKES = 2
};

enum class BOUNDARY{
  XFEM = 0,
  DARCY = 1
};


class FEM :public DomainFEM{
  public:

    SOLVER solver;
    BOUNDARY bd;

    SolutionData  SolnData;
    ElementBaseFEM **elm;
    PetscSolver  *solverPetsc;
    
    TextParser tp;
    PetscErrorCode  errpetsc;
    string outputDir,fileName;
    
    int numOfDofsNode=4;
    PetscInt  numOfId, myId, nNode_owned;
    PetscInt  node_start, node_end, elem_start, elem_end;
    PetscInt  row_start, row_end;
    PetscInt numOfDofsLocal, numOfDofsGlobal;

    vector<int>  assyForSoln, OutputNodes;
    vector<double> DirichletBCs; 
    vector<vector<double>> DirichletBCs_tmp; 

    vector<vector<int>>  nodeDofArrayBCsPrev, nodeDofArrayPrev, nodeDofArrayBCs, nodeDofArray;
    vector<vector<bool>>  nodeTypePrev, nodeType;

    int numOfOMP;

    double rho,mu,nu,U,D,Re;
    
    double NRtolerance;
    int NRitr_initial, NRitr; 
    double relaxationParam, relaxationParam_initial;

    double dt;
    int timeMax;
    
    vector<double> phi, phiEX, phiVOF;
    vector<double> sdf, sdf_node;

    vector<double> u;
    vector<double> v;
    vector<double> w;
    vector<double> p;

    /// XFEM ///
    int max_depth;
    int sub_div;
    int sub_div_total;
    
    /// Darcy ///
    double resistance;
    double alpha;

  ///// FLUID ONLY /////
  public:
    PetscInt numOfDofsLocalFluid, numOfDofsGlobalFluid;
    vector<int>  assyForSolnFluid;
    vector<double> DirichletBCsFluid; 
    vector<vector<double>> DirichletBCsFluid_tmp; 

    vector<vector<int>>  nodeDofArrayBCsPrevFluid, nodeDofArrayPrevFluid, nodeDofArrayBCsFluid, nodeDofArrayFluid;
    vector<vector<bool>>  nodeTypePrevFluid, nodeTypeFluid;

    SolutionData  SolnDataFluid;
    ElementBaseFEM **elmFluid;
    PetscSolver  *solverPetscFluid;
    
    vector<double> phiFluid;
    vector<double> phiEXFluid;
    vector<double> phiVOFFluid;
    vector<double> sdfFluid_node;
    vector<double> sdfFluid;
    vector<int> sortElm, sortNode;

    vector<double> uFluid;
    vector<double> vFluid;
    vector<double> wFluid;
    vector<double> pFluid;

    vector<vector<double>> uf;
    vector<vector<double>> vf;
    vector<vector<double>> wf;
    vector<vector<double>> pf;



  public:
    FEM();
    ~FEM();
    void initialize();
    void readInput();
    void readMPI();
    void readOutput();
    void readPysicalParam();
    void readBoundaryMethod();
    void readXFEMParam();
    void readDarcyParam();
    void readBTSubDivParam();
    void readNRParam();
    void readTimeParam();
    void readDomain();
    void readBoundary();
    void readImage();
    void setDomain();
    void prepare();
    void setBoundary();
    void setFluidDomain();

    void prepareMatrix();
    int divideMesh();
    int prepareForParallel();

    int deallocate();

    /// XFEM PARTITION ///
    void binaryTreeSubDivision();
    void gererateSubElms(vector<double> &sdf_parent, vector<vector<double>> &x_sub, vector<double> &x_center_parent, const int &ic, int &tmp, int &depth);
    void getSubSubCoordinates(vector<vector<double>> &x_sub, vector<vector<double>> &x_sub_sub, vector<double> &sdf_sub_sub, const int &ii, const int &jj, const int &kk);
    void makeSubElmsData(vector<double> &sdf_sub_sub, vector<vector<double>> &x_sub_sub, vector<double> &x_center_parent, const int &ic, int &tmp, int &depth);
    void subGaussPoint(vector<vector<double>> &x_sub_sub, double &sub_gx_tmp, double &sub_gy_tmp, double &sub_gz_tmp, vector<double> &x_center_parent, const int &ii, const int &jj, const int &kk);
    void subGaussWeight(double &weight, const int &ii, const int &jj, const int &kk, int &depth);


    void interfacePartition();
    void line_serch(vector<int> &sub_cross_num, vector<vector<double>> &sub_x_tmp, vector<double> &sdf_current, vector<vector<double>> &x_current, const int ic);
    void s_t_u(const int &linePartition, const int &i, const int &j, double &s, double &t, double &u);
    bool is_cross(vector<double> &sdf_current, const int i);
    void tetraPartition1(vector<vector<double>> &sub_x_tmp, vector<double> &sdf_current, vector<vector<double>> &x_current, const int ic);
    void tetraPartition2(vector<vector<double>> &sub_x_tmp, vector<double> &sdf_current, vector<vector<double>> &x_current, const int ic);
    void tetraPartition3(vector<vector<double>> &sub_x_tmp, vector<double> &sdf_current, vector<vector<double>> &x_current, const int ic);

    /// STEADY STOKES  ///
    void Stokes();
   
    void MatAssySTT(const int ic,MatrixXd &Klocal, VectorXd &Flocal);
    void XFEM_MatAssySTT(const int ic,MatrixXd &Klocal, VectorXd &Flocal);
    void XFEM_MatAssySTT2(const int ic,MatrixXd &Klocal, VectorXd &Flocal);
    void XFEM_MatAssySTT3(const int ic, MatrixXd &Klocal, VectorXd &Flocal);
    void SAWADA_XFEM_MatAssySTT(const int ic, MatrixXd &Klocal, VectorXd &Flocal);
    void Darcy_MatAssySTT(const int ic,MatrixXd &Klocal, VectorXd &Flocal);
    
    /// STEADY NAVIER STOKES  ///
    void SteadyNavierStokes();
   
    void MatAssySNS(const int ic,MatrixXd &Klocal, VectorXd &Flocal);
    void XFEM_MatAssySNS(const int ic,MatrixXd &Klocal, VectorXd &Flocal);
    
    void DiffusionInGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic);
    void PressureInGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<double> &N,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic);
    void PSPGInGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic);
    
    void localRefinement(const int n,std::vector<double> &b);
    void DiffusionInGaussIntegralXFEM(MatrixXd &Klocal, VectorXd &Flocal, vector<vector<double>> &dNPdr, vector<vector<double>> &dNVdr, vector<vector<double>> &x_current, const int numOfNodeInElm,const double weight,const int ic);
    void PressureInGaussIntegralXFEM(MatrixXd &Klocal, VectorXd &Flocal, vector<double> &NP, vector<vector<double>> &dNPdr, vector<vector<double>> &dNVdr, vector<vector<double>> &x_current, const int numOfNodeInElm, const double weight, const int ic);
    void PSPGInGaussIntegralXFEM(MatrixXd &Klocal, VectorXd &Flocal,vector<vector<double>> &dNPdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic);
   
    /// UNSTEADY NAVIER STOKES  ///
    void UnsteadyNavierStokes();
    
    void MatAssyUSNS(MatrixXd &Klocal, VectorXd &Flocal, const int ic, const int t_itr);
    void XFEM_MatAssyUSNS(MatrixXd &Klocal, VectorXd &Flocal, const int ic, const int t_itr);
    void Darcy_MatAssyUSNS(MatrixXd &Klocal, VectorXd &Flocal, const int ic, const int t_itr);
    void VelocityValue(double (&vel)[3], double (&advel)[3], double (&dvdx)[3][3],vector<double> &N, vector<vector<double>> &dNdx, const int ic, const int t_itr);
    
    void assignBCs();
    void assignPulsatileBCs(const double t_itr);
    void applyBCs();
    double calc_tau(const double (&dxdr)[3][3],const double (&vel)[3]);
    double calc_tau2(const double (&vel)[3]);
  
    void postCaluculation();
    void postCaluculation_itr(const int loop);
    void postCaluculation_timeItr(const int t_itr);
    
    void export_vti_metis(const string &file, vector<int> &node, vector<int> &element);
    void export_vti_node(const string &file, vector<double> &node);
    void export_vti_elm(const string &file, vector<double> &element);
    void export_vti_domain(const string &file);
    void export_vti_result(const std::string &file, vector<double> &u, vector<double> &v, vector<double> &w, vector<double> &p);
    void export_vti_result_2D(const std::string &file, vector<double> &u, vector<double> &v, vector<double> &p);
};

