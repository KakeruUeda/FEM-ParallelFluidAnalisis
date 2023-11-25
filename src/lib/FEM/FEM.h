#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "metis.h"
#include "TextParser.h"
#include "domainFEM.h"
#include "petscSolver.h"

using namespace std;

class FEM :public DomainFEM{
  public:
    TextParser tp;
    PetscErrorCode  errpetsc;

    PetscInt  n_mpi_procs, this_mpi_proc;
    
    ElementBaseFEM  **elm;
  
  public:
    FEM();
    ~FEM();
    void initialize();
    void readInput();
    void setDomain();
    void prepare();
    int partitionMesh();

    void export_vti(const string &file, vector<int> &node, vector<int> &element);
};

