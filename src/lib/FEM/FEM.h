#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <omp.h>

#include "gauss.h"
#include "metis.h"
#include "TextParser.h"
#include "DomainFEM.h"
#include "PetscSolver.h"
//#include "solutionData.h"
#include "ShapeFunction.h"
#include "ElementBaseFEM.h"
#include "MathFEM.h"
using namespace std;


class FEM :public DomainFEM{
  public:
    TextParser tp;
    PetscErrorCode  errpetsc;

    int numOfDofsNode=4;
    PetscInt  numOfId, myId, nNode_owned;
    PetscInt  node_start, node_end, elem_start, elem_end;
    PetscInt  row_start, row_end;
    //PetscInt numOfDofsLocal, numOfDofsGlobal;

    //vector<int>  assyForSoln, OutputNodes;
    //vector<double> DirichletBCs; 
    //vector<vector<double>> DirichletBCs_tmp; 

    //vector<vector<int>>  nodeDofArrayBCsPrev, nodeDofArrayPrev, nodeDofArrayBCs, nodeDofArray;
    //vector<vector<bool>>  nodeTypePrev, nodeType;
    
    //SolutionData  SolnData;
    //ElementBaseFEM **elm;
    //PetscSolver  *solverPetsc;

    double rho,mu,nu;
  
    vector<double> phi;
    vector<double> phiEX;
    vector<double> sdf;


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
    vector<double> sdfFluid;

    vector<int> sortElm, sortNode;

  public:
    FEM();
    ~FEM();
    void initialize();
    void readInput();
    void setDomain();
    void prepare();
    void setBoundary();
    void setFluidDomain();

    void prepareMatrix();
    int divideMesh();
    int prepareForParallel();

    int deallocate();

    void stokes();
   
    void calcStokesMatrix(const int ic,MatrixXd &Klocal, VectorXd &Flocal);
    void calcStokesMatrixXFEM(const int ic,MatrixXd &Klocal, VectorXd &Flocal);
    
    void DiffusionInGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic);
    void PressureInGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<double> &N,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic);
    void PSPGInGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic);
    
    void localRefinement(const int n,std::vector<double> &b);
    void DiffusionInGaussIntegralXFEM(MatrixXd &Klocal, VectorXd &Flocal, vector<vector<double>> &dNPdr, vector<vector<double>> &dNVdr, vector<vector<double>> &x_current, const int numOfNodeInElm,const double weight,const int ic);
    void PressureInGaussIntegralXFEM(MatrixXd &Klocal, VectorXd &Flocal, vector<double> &NP, vector<vector<double>> &dNPdr, vector<vector<double>> &dNVdr, vector<vector<double>> &x_current, const int numOfNodeInElm, const double weight, const int ic);
    void PSPGInGaussIntegralXFEM(MatrixXd &Klocal, VectorXd &Flocal,vector<vector<double>> &dNPdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic);
    
    void assignBoundaryConditions();
    void applyBoundaryConditions(VectorXd& Flocal, const int size, const int ic);
    void getSolution();
    

    void postCaluculation();
    
    void export_vti(const string &file, vector<int> &node, vector<int> &element);
    void export_vti_domain(const string &file);
    void export_vti_result(const std::string &file, vector<double> &u, vector<double> &v, vector<double> &w, vector<double> &p);
    void export_vti_result_2D(const std::string &file, vector<double> &u, vector<double> &v, vector<double> &p);

};

