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
#include "domainFEM.h"
#include "petscSolver.h"
//#include "solutionData.h"
#include "ShapeFunction.h"
#include "elementBase.h"
#include "mathFEM.h"
using namespace std;


class FEM :public DomainFEM{
  public:
    TextParser tp;
    PetscErrorCode  errpetsc;

    int ndof=4;
    PetscInt  n_mpi_procs, this_mpi_proc, nNode_owned;
    PetscInt  node_start, node_end, elem_start, elem_end;
    PetscInt  row_start, row_end, ntotdofs_local, ntotdofs_global;


    vector<int>  assyForSoln, OutputNodes;
    vector<vector<double>>  DirichletBCs; 

    vector<vector<int>>  NodeDofArrayOld, NodeDofArrayOld_withoutBd, NodeDofArrayNew;
    vector<vector<bool>>  NodeTypeOld, NodeTypeNew;
    
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
    int divideMesh();
    int prepareForParallel();
    
    void assignBoundaryConditions();
    void Stokes();
    void calcStokesMatrix(const int ic,MatrixXd &Klocal, VectorXd &Flocal);

    void Pressure_inGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<double> &N,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic);
    void Diffusion_inGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic);
    void PSPG_inGaussIntegral(MatrixXd &Klocal, VectorXd &Flocal,vector<vector<double>> &dNdr,vector<vector<double>> &x_current,const int numOfNodeInElm,const double weight,const int ic);

    void export_vti(const string &file, vector<int> &node, vector<int> &element);
    void export_vti_result(const std::string &file, vector<double> &u, vector<double> &v, vector<double> &w, vector<double> &p);
    void export_vti_result_2D(const std::string &file, vector<double> &u, vector<double> &v, vector<double> &p);
};

