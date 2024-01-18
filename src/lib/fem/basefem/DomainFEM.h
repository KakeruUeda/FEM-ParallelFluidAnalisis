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
    VDOUBLE2D x;
    VINT2D element;
    VINT2D elementPrev;

    VDOUBLE1D bd_p;
    VINT1D bd_ip;
    VDOUBLE2D bd_u;
    VINT2D bd_iu;

    VINT1D  nodeId, nodeMapPrev, nodeMap;
    VINT1D  elmId;
    
    
     /// FLUID ONLY ///
  public:
    PetscInt numOfNodeGlobalFluid,numOfElmGlobalFluid;
    PetscInt numOfNodeLocalFluid,numOfElmLocalFluid;
    PetscInt numOfBdNodeFluid;
    
    VDOUBLE2D xFluid;
    VINT2D elementFluid;
    VINT2D elementFluidPrev;

    VDOUBLE1D bd_p_fluid;
    VINT1D bd_ip_fluid;
    VDOUBLE2D bd_u_fluid;
    VINT2D bd_iu_fluid;

    VINT1D  nodeIdFluid, nodeMapPrevFluid, nodeMapFluid;
    VINT1D  elmIdFluid;

  private:
};

