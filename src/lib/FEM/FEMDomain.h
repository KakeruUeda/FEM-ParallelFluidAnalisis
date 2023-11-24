#include <string>
#include "define.h"
#include "petscSolver.h"
using namespace std;

class FEMDomain{
  public:
    int nx,ny,nz;
    double dx,dy,dz,Lx,Ly,Lz;
    PetscInt numOfNode,numOfElm,numOfBd,numOfBdWall;

    vector<vector<double>> x;
    vector<ElementType> element;

  private:
};