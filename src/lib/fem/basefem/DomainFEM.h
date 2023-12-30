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
    vector<vector<int>> elementPrev;

    vector<double> bd_p;
    vector<int> bd_ip;
    vector<vector<double>> bd_u;
    vector<vector<int>> bd_iu;

    vector<int>  nodeId, nodeMapPrev, nodeMap;
    
    
     /// FLUID ONLY ///
  public:
    PetscInt numOfNodeGlobalFluid,numOfElmGlobalFluid;
    PetscInt numOfNodeLocalFluid,numOfElmLocalFluid;
    PetscInt numOfBdNodeFluid;
    
    vector<vector<double>> xFluid;
    vector<vector<int>> elementFluid;
    vector<vector<int>> elementFluidPrev;

    vector<double> bd_p_fluid;
    vector<int> bd_ip_fluid;
    vector<vector<double>> bd_u_fluid;
    vector<vector<int>> bd_iu_fluid;

    vector<int>  nodeIdFluid, nodeMapPrevFluid, nodeMapFluid;

  private:
};

