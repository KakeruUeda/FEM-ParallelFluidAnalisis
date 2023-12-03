#include "Eigen.h"
#include "petscksp.h"
#include "petscmat.h"

using namespace std;

enum { SOLVER_EMPTY, PATTERN_OK, INIT_OK, ASSEMBLY_OK, FACTORISE_OK };

class PetscSolver
{
  public:
    Vec  rhsVec, solnVec, solnVecPrev, reacVec;
    Mat  mtx; // linear system matrix
    KSP  ksp; // linear solver context
    PC   pc; // preconditioner context
    
    PetscInt nRow, nCol, nnz;

    int currentStatus;

    bool  checkIO;

    PetscReal norm; // norm of solution error

    PetscErrorCode errpetsc;

    
    PetscSolver();
    virtual ~PetscSolver();

    virtual int initialise(int size_local, int size_global, int* diag_nnz, int* offdiag_nnz);
    virtual int zeroMtx();
    virtual int assembleMatrixAndVectorSerial(vector<int>& forAssyElem, vector<int>& forAssyElemRHS, MatrixXd& Klocal, VectorXd& Flocal);
    virtual int factorise();
    virtual int solve();
    virtual int factoriseAndSolve();
    virtual int free();
};