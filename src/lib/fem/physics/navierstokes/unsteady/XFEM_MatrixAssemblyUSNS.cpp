#include "FEM.h"

void FEM::XFEM_MatAssyUSNS(MatrixXd &Klocal, VectorXd &Flocal, const int ic, const int t_itr)
{
  PetscPrintf(MPI_COMM_WORLD, "XFEM not yet implemented");
  exit(1);
}