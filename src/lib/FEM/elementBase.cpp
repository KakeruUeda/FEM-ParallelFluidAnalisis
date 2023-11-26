#include "elementBase.h"

using namespace std;

ElementBaseFEM::ElementBaseFEM()
{
  subdomainId = 10;
}


ElementBaseFEM::~ElementBaseFEM()
{
//  cout << "     ElementBase: destructor ...\n\n";
}

void ElementBaseFEM::prepareElemData(const int nsize)
{
  int ndof = 4;
  int numOfNodeInElm = 8;
  int ii;
  globalDOFnums.resize(nsize);
  for(ii=0; ii<numOfNodeInElm; ii++)
  {
    globalDOFnums[ii*ndof]   = nodeNums[ii]*ndof;
    globalDOFnums[ii*ndof+1] = nodeNums[ii]*ndof + 1;
    globalDOFnums[ii*ndof+2] = nodeNums[ii]*ndof + 2;
    globalDOFnums[ii*ndof+3] = nodeNums[ii]*ndof + 3;
  }
}