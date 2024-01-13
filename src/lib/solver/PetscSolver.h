#include "Eigen.h"
#include "petscksp.h"
#include "petscmat.h"
#include "define.h"


enum { SOLVER_EMPTY, PATTERN_OK, INIT_OK, ASSEMBLY_OK, FACTORISE_OK };

class PetscSolver
{
  public:
    Vec  rhsVec, solnVec, solnVecPrev;
    Mat  mtx; // linear system matrix
    KSP  ksp; // linear solver context
    PC   pc; // preconditioner context
    
    PetscInt nRow, nCol, nnz;

    int nnz_max_row;
    PetscInt  *diag_nnz, *offdiag_nnz;

    int currentStatus;

    bool  checkIO;

    PetscReal norm; // norm of solution error

    PetscErrorCode errpetsc;
     
    double vectorNorm(const int &nump, VectorXd &x);

    PetscSolver();
    virtual ~PetscSolver();

    int initialize(int size_local, int size_global);
    int initialAssembly();
    int setZeroValue();
    int setValue(VINT1D& forAssyElem, VINT1D& forAssyElemRHS, MatrixXd& Klocal, VectorXd& Flocal);
    int factorise();
    int solve();
    int factoriseAndSolve();
};