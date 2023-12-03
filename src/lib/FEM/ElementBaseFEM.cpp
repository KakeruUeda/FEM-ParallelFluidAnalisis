#include "ElementBaseFEM.h"

using namespace std;

ElementBaseFEM::ElementBaseFEM()
{
  subdomainId = 0;
  numOfDofsNode = 4;
  numOfNodeInElm = 8;
}


ElementBaseFEM::~ElementBaseFEM()
{
//  cout << "     ElementBase: destructor ...\n\n";
}

void ElementBaseFEM::prepareElemData(const int nsize)
{
  int ii;
  globalDOFnums.resize(nsize);
  for(ii=0; ii<numOfNodeInElm; ii++)
  {
    globalDOFnums[ii*numOfDofsNode]   = nodeNums[ii]*numOfDofsNode;
    globalDOFnums[ii*numOfDofsNode+1] = nodeNums[ii]*numOfDofsNode + 1;
    globalDOFnums[ii*numOfDofsNode+2] = nodeNums[ii]*numOfDofsNode + 2;
    globalDOFnums[ii*numOfDofsNode+3] = nodeNums[ii]*numOfDofsNode + 3;
  }
}