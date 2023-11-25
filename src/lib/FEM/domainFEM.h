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
    vector<vector<double>> x;
    //vector<ElementType> element;
    vector<vector<int>> element;
    //vector<vector<int>> elmConnect; 
    vector<int>  node_proc_id, node_map_get_old, node_map_get_new;

  private:
};

class ElementBaseFEM{
  public:
    int subdomainId;

    void  setSubdomainId(int id)
    {  subdomainId = id; return;  }

  private:
};