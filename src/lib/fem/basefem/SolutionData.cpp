#include "SolutionData.h"

SolutionData::SolutionData()
{
  td.resize(100);
  return;
}

void SolutionData::initialize(int size1)
{
  vecsize = size1;

  soln.resize(vecsize);
  soln.setZero();

  solnInit  = soln;
  solnPrev  = soln;
  solnPrev2 = soln;
  solnPrev3 = soln;
  solnPrev4 = soln;
  solnCur   = soln;

  solnDot     = soln;
  solnDotPrev = soln;
  solnDotCur  = soln;
  solnApplied = soln;

  return;
}

void SolutionData::setZero()
{
  soln.setZero();

  return;
}