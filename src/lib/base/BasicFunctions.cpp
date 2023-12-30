#include "BasicFunctions.h"

void BasicFunctions::calcInverseMatrix_3x3(double (&inv_a)[3][3],const double (&a)[3][3])
{
  double det;

  det = calcDeterminant_3x3(a);

  inv_a[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
  inv_a[0][1] = a[0][2]*a[2][1] - a[0][1]*a[2][2];
  inv_a[0][2] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
  inv_a[1][0] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
  inv_a[1][1] = a[0][0]*a[2][2] - a[0][2]*a[2][0];
  inv_a[1][2] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
  inv_a[2][0] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
  inv_a[2][1] = a[0][1]*a[2][0] - a[0][0]*a[2][1];
  inv_a[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) inv_a[i][j] = inv_a[i][j] / det;
  }
}

double BasicFunctions::calcDeterminant_3x3(const double (&a)[3][3])
{
  double det  = a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2] + a[2][0] * a[0][1] * a[1][2]
              - a[2][0] * a[1][1] * a[0][2] - a[1][0] * a[0][1] * a[2][2] - a[0][0] * a[2][1] * a[1][2];
  return det;
}