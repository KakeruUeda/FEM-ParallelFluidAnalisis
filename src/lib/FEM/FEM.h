#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "FEMDomain.h"

class FEM :public FEMDomain{
  public:
   PetscInt  n_mpi_procs, this_mpi_proc;
  public:
    FEM();
    ~FEM();
    void initialize(string& domainFile);
    void readInput(string& domainFile);
    void setDomain();
    void partitionMesh();
};

