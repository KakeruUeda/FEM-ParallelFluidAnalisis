#include "solutionData.h"

class ElementBaseFEM{
  public:
    int subdomainId;
    vector<int>  nodeNums,forAssyVec,forAssyVec_withoutBd,globalDOFnums;

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