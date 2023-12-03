#include "solutionData.h"

class ElementBaseFEM{
  public:
    int subdomainId,numOfDofsNode,numOfNodeInElm;
    vector<int>  nodeNums,nodeNumsPrev,nodeForAssyBCs,nodeForAssy,globalDOFnums;

    ElementBaseFEM();
    virtual ~ElementBaseFEM();
    
    SolutionData  *SolnData;

    void  setSubdomainId(int sid)
    {  subdomainId = sid; return;  }

    int getSubdomainId()
    {  return  subdomainId;  }
    
    void prepareElemData(const int nsize);

  private:
};