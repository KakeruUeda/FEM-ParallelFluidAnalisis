#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>

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
    PetscInt  row_start, row_end, numOfDofsLocal, numOfDofsGlobal;


    vector<int>  assyForSoln, OutputNodes;
    vector<double> DirichletBCs; 
    vector<vector<double>> DirichletBCs_tmp; 

    vector<vector<int>>  nodeDofArrayBCsPrev, nodeDofArrayPrev, nodeDofArrayBCs, nodeDofArray;
    vector<vector<bool>>  nodeTypePrev, nodeType;
    
    SolutionData  SolnData;
    ElementBaseFEM **elm;
    PetscSolver  *solverPetsc;

    double rho,mu,nu;
  
  public:
    FEM();
    ~FEM();
    void initialize();
    void readInput();
    void setDomain();
    void prepare();
    void setBoundary();

    void prepareMatrix();
    int divideMesh();
    int prepareForParallel();

    int deallocate();

    void Stokes();
    void calcStokesMatrix(const int ic,MatrixXd &Klocal, VectorXd &Flocal);

    void Pressure_inGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<double> &N,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic);
    void Diffusion_inGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic);
    void PSPG_inGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic);
    
    void assignBoundaryConditions();
    void applyBoundaryConditions(VectorXd& Flocal, const int size, const int ic);
    void getSolution();
    
    void export_vti(const string &file, vector<int> &node, vector<int> &element);
    void export_vti_result(const std::string &file, vector<double> &u, vector<double> &v, vector<double> &w, vector<double> &p);
    void export_vti_result_2D(const std::string &file, vector<double> &u, vector<double> &v, vector<double> &p);
};

