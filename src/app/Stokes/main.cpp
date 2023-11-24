#include "FEM.h"
 
int main(int argc, char* argv[])
{
  PetscInitialize(NULL, NULL,(char *)0, NULL);

  FEM Fem;

 string  domainFile    = argv[1];
 //string  controlfile = argv[2];
 //string  petscfile   = argv[3];

  Fem.initialize(domainFile);

  PetscFinalize();
  return 0;
}