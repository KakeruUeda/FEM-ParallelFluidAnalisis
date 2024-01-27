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
  /* 
  int ii;
  globalDOFnums.resize(nsize);
  for(ii=0; ii<numOfNodeInElm; ii++)
  {
    globalDOFnums[ii*numOfDofsNode]   = nodeNums[ii]*numOfDofsNode;
    globalDOFnums[ii*numOfDofsNode+1] = nodeNums[ii]*numOfDofsNode + 1;
    globalDOFnums[ii*numOfDofsNode+2] = nodeNums[ii]*numOfDofsNode + 2;
    globalDOFnums[ii*numOfDofsNode+3] = nodeNums[ii]*numOfDofsNode + 3;
  }
  */

  int ii;
  globalDOFnumsFluid.resize(nsize);
  for(ii=0; ii<numOfNodeInElm; ii++)
  {
    globalDOFnumsFluid[ii*numOfDofsNode]   = nodeNumsFluid[ii]*numOfDofsNode;
    globalDOFnumsFluid[ii*numOfDofsNode+1] = nodeNumsFluid[ii]*numOfDofsNode + 1;
    globalDOFnumsFluid[ii*numOfDofsNode+2] = nodeNumsFluid[ii]*numOfDofsNode + 2;
    globalDOFnumsFluid[ii*numOfDofsNode+3] = nodeNumsFluid[ii]*numOfDofsNode + 3;
  }
}

void ElementBaseFEM::prepareAdjointElemDataBd(VINT2D &numOfDofsNodeInElementAdjointFluid, int nsize, const int ie)
{ 
  int tmp = 0;
  globalDOFnumsAdjointFluid.resize(nsize);
  for(int p=0; p<numOfNodeInElm; p++)
  {  
    if(numOfDofsNodeInElementAdjointFluid[ie][p] > numOfDofsNode){
      globalDOFnumsAdjointFluid[tmp]     = dofsNumsFluid[p];
      globalDOFnumsAdjointFluid[tmp + 1] = dofsNumsFluid[p] + 1;
      globalDOFnumsAdjointFluid[tmp + 2] = dofsNumsFluid[p] + 2;
      globalDOFnumsAdjointFluid[tmp + 3] = dofsNumsFluid[p] + 3;
      globalDOFnumsAdjointFluid[tmp + 4] = dofsNumsFluid[p] + 4;
      globalDOFnumsAdjointFluid[tmp + 5] = dofsNumsFluid[p] + 5;
      globalDOFnumsAdjointFluid[tmp + 6] = dofsNumsFluid[p] + 6; 
      tmp = tmp + numOfDofsNode + 3;
    }else{
      globalDOFnumsAdjointFluid[tmp]     = dofsNumsFluid[p];
      globalDOFnumsAdjointFluid[tmp + 1] = dofsNumsFluid[p] + 1;
      globalDOFnumsAdjointFluid[tmp + 2] = dofsNumsFluid[p] + 2;
      globalDOFnumsAdjointFluid[tmp + 3] = dofsNumsFluid[p] + 3;
      tmp = tmp + numOfDofsNode;
    }
  }

}

void ElementBaseFEM::prepareAdjointElemData(const int nsize)
{
  int ii;
  globalDOFnumsAdjointFluid.resize(nsize);
  for(ii=0; ii<numOfNodeInElm; ii++)
  {
    globalDOFnumsAdjointFluid[ii * numOfDofsNode]     = nodeNumsFluid[ii] * numOfDofsNode;
    globalDOFnumsAdjointFluid[ii * numOfDofsNode + 1] = nodeNumsFluid[ii] * numOfDofsNode + 1;
    globalDOFnumsAdjointFluid[ii * numOfDofsNode + 2] = nodeNumsFluid[ii] * numOfDofsNode + 2;
    globalDOFnumsAdjointFluid[ii * numOfDofsNode + 3] = nodeNumsFluid[ii] * numOfDofsNode + 3;
  }
}