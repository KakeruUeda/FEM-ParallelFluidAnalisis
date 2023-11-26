#include <string>
#include "define.h"
#include "petscksp.h"
#include "petscmat.h"



using namespace std;

class DomainFEM{
  public:
    int nx,ny,nz;
    double dx,dy,dz,Lx,Ly,Lz;
    PetscInt numOfNodeGlobal,numOfElmGlobal;
    PetscInt numOfNodeLocal,numOfElmLocal;
    PetscInt numOfNodeInElm,numOfBd,numOfBdWall;
    PetscInt numOfBdNode;
    vector<vector<double>> x;
    vector<vector<int>> element;

    vector<double> bd_p;
    vector<int> bd_ip;
    vector<vector<double>> bd_u;
    vector<vector<int>> bd_iu;

    vector<int>  node_proc_id, node_map_get_old, node_map_get_new;

  private:
};

