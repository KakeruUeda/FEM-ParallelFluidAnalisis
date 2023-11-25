#include "eigen.h"
#include "petscksp.h"
#include "petscmat.h"

using namespace std;

class Petsc
{
  public:
    Vec  rhsVec, solnVec, solnVecPrev, reacVec;
    Mat  mtx; // linear system matrix
    KSP  ksp; // linear solver context
    PC   pc; // preconditioner context
};