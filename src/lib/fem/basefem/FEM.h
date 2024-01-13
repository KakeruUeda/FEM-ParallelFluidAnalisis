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

    int numOfDofsNode = 4;
    PetscInt  numOfId, myId, nNode_owned;
    PetscInt  node_start, node_end, elem_start, elem_end;
    PetscInt  row_start, row_end;
    PetscInt numOfDofsLocal, numOfDofsGlobal;

    VINT1D  assyForSoln, OutputNodes;
    VDOUBLE1D DirichletBCs; 
    VDOUBLE2D DirichletBCs_tmp; 

    VINT2D  nodeDofArrayBCsPrev, nodeDofArrayPrev, nodeDofArrayBCs, nodeDofArray;
    VBOOL2D  nodeTypePrev, nodeType;

    int numOfOMP;

    double rho, mu, nu, Re;
    
    double NRtolerance;
    int NRitr_initial, NRitr; 
    double relaxationParam, relaxationParam_initial;

    double dt;
    int timeMax;
    int pulsatile_flow;
    int pulse_begin_itr;
    double T;
    
    VDOUBLE1D phi, phiEX, phiVOF;
    VDOUBLE1D sdf, sdf_node;

    VDOUBLE1D u;
    VDOUBLE1D v;
    VDOUBLE1D w;
    VDOUBLE1D p;

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
    VINT1D  assyForSolnFluid;
    VDOUBLE1D DirichletBCsFluid; 
    VDOUBLE2D DirichletBCsFluid_tmp; 

    VINT2D  nodeDofArrayBCsPrevFluid, nodeDofArrayPrevFluid, nodeDofArrayBCsFluid, nodeDofArrayFluid;
    VBOOL2D  nodeTypePrevFluid, nodeTypeFluid;

    SolutionData  SolnDataFluid;
    ElementBaseFEM **elmFluid;
    PetscSolver  *solverPetscFluid;
    
    VDOUBLE1D phiFluid;
    VDOUBLE1D phiEXFluid;
    VDOUBLE1D phiVOFFluid;
    VDOUBLE1D sdfFluid_node;
    VDOUBLE1D sdfFluid;
    VINT1D sortElm, sortNode;

    VDOUBLE1D uFluid;
    VDOUBLE1D vFluid;
    VDOUBLE1D wFluid;
    VDOUBLE1D pFluid;

    VDOUBLE2D uf;
    VDOUBLE2D vf;
    VDOUBLE2D wf;
    VDOUBLE2D pf;

  public:
    FEM();
    ~FEM();
    
    void initialize();         void readInput();
    void readBase();           void readPysicalParam();
    void readBoundaryMethod(); void readXFEMParam();
    void readDarcyParam();     void readBTSubDivParam();
    void readNRParam();        void readTimeParam();
    void readDomain();         void readBoundary();
    void readImage();          void setDomain();
    void setBoundary();        void setFluidDomain();

    void prepareMatrix();
    void exportDomain();
    
    int divideMesh();
    int prepareForParallel();
    
    void allocateObj();
    void resizeVariables();

    int solverDeallocate();

    /// XFEM PARTITION ///
    void octreeSubDivision();
    void gererateSubElms(VDOUBLE1D &sdf_parent, VDOUBLE2D &x_sub, VDOUBLE1D &x_center_parent, const int &ic, int &tmp, int &depth);
    void getSubSubCoordinates(VDOUBLE2D &x_sub, VDOUBLE2D &x_sub_sub, VDOUBLE1D &sdf_sub_sub, const int &ii, const int &jj, const int &kk);
    void makeSubElmsData(VDOUBLE1D &sdf_sub_sub, VDOUBLE2D &x_sub_sub, VDOUBLE1D &x_center_parent, const int &ic, int &tmp, int &depth);
    void subGaussPoint(VDOUBLE2D &x_sub_sub, double &sub_gx_tmp, double &sub_gy_tmp, double &sub_gz_tmp, VDOUBLE1D &x_center_parent, const int &ii, const int &jj, const int &kk);
    void subGaussWeight(double &weight, const int &ii, const int &jj, const int &kk, int &depth);


    void interfacePartition();
    void line_serch(VINT1D &sub_cross_num, VDOUBLE2D &sub_x_tmp, VDOUBLE1D &sdf_current, VDOUBLE2D &x_current, const int ic);
    void s_t_u(const int &linePartition, const int &i, const int &j, double &s, double &t, double &u);
    bool is_cross(VDOUBLE1D &sdf_current, const int i);
    void tetraPartition1(VDOUBLE2D &sub_x_tmp, VDOUBLE1D &sdf_current, VDOUBLE2D &x_current, const int ic);
    void tetraPartition2(VDOUBLE2D &sub_x_tmp, VDOUBLE1D &sdf_current, VDOUBLE2D &x_current, const int ic);
    void tetraPartition3(VDOUBLE2D &sub_x_tmp, VDOUBLE1D &sdf_current, VDOUBLE2D &x_current, const int ic);

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
   
    void MatAssySNS(const int ic, MatrixXd &Klocal, VectorXd &Flocal);
    void XFEM_MatAssySNS(const int ic, MatrixXd &Klocal, VectorXd &Flocal);
    void Darcy_MatAssySNS(const int ic, MatrixXd &Klocal, VectorXd &Flocal);
    
    void DiffusionInGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,VDOUBLE2D &dNdr,VDOUBLE2D &x_current,const int numOfNodeInElm,const double weight,const int ic);
    void PressureInGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,VDOUBLE1D &N,VDOUBLE2D &dNdr,VDOUBLE2D &x_current,const int numOfNodeInElm,const double weight,const int ic);
    void PSPGInGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,VDOUBLE2D &dNdr,VDOUBLE2D &x_current,const int numOfNodeInElm,const double weight,const int ic);
    
    void LocalRefinement(const int n, VDOUBLE1D &b);
    void DiffusionInGaussIntegralXFEM(MatrixXd &Klocal, VectorXd &Flocal, VDOUBLE2D &dNPdr, VDOUBLE2D &dNVdr, VDOUBLE2D &x_current, const int numOfNodeInElm,const double weight,const int ic);
    void PressureInGaussIntegralXFEM(MatrixXd &Klocal, VectorXd &Flocal, VDOUBLE1D &NP, VDOUBLE2D &dNPdr, VDOUBLE2D &dNVdr, VDOUBLE2D &x_current, const int numOfNodeInElm, const double weight, const int ic);
    void PSPGInGaussIntegralXFEM(MatrixXd &Klocal, VectorXd &Flocal,VDOUBLE2D &dNPdr,VDOUBLE2D &x_current,const int numOfNodeInElm,const double weight,const int ic);
   
    /// UNSTEADY NAVIER STOKES  ///
    void UnsteadyNavierStokes();
    
    void MatAssyUSNS(MatrixXd &Klocal, VectorXd &Flocal, const int ic, const int t_itr);
    void XFEM_MatAssyUSNS(MatrixXd &Klocal, VectorXd &Flocal, const int ic, const int t_itr);
    void Darcy_MatAssyUSNS(MatrixXd &Klocal, VectorXd &Flocal, const int ic, const int t_itr);
    void velocityValue(double (&vel)[3], double (&advel)[3], double (&dvdx)[3][3],VDOUBLE1D &N, VDOUBLE2D &dNdx, const int ic, const int t_itr);
    
    void setMatAndVecZero();
    void setNRInitialValue();
    
    void assignBCs();
    void assignPulsatileBCs(const double t_itr);
    void applyBCs();
    double calc_tau(const double (&dxdr)[3][3],const double (&vel)[3]);
    double calc_tau2(const double (&vel)[3]);
    double calc_tau3(VDOUBLE2D &dNdx, double (&vel)[3]);
  
    void postCaluculation();
    void postCaluculation_itr(const int loop);
    void postCaluculation_timeItr(const int t_itr);
    
    void export_vti_metis(const string &file, VINT1D &node, VINT1D &element);
    void export_vti_node(const string &file, VDOUBLE1D &node);
    void export_vti_elm(const string &file, VDOUBLE1D &element);
    void export_vti_domain(const string &file);
    void export_vti_result(const std::string &file, VDOUBLE1D &u, VDOUBLE1D &v, VDOUBLE1D &w, VDOUBLE1D &p);
    void export_vti_result_2D(const std::string &file, VDOUBLE1D &u, VDOUBLE1D &v, VDOUBLE1D &p);
};

