#include "solutionData.h"

class ElementBaseFEM{
  public:
    int subdomainId,numOfDofsNode,numOfNodeInElm,numOfSubdomain;
    //vector<int>  nodeNums,nodeNumsPrev,nodeForAssyBCs,nodeForAssy,globalDOFnums;

    ElementBaseFEM();
    virtual ~ElementBaseFEM();
    
    //SolutionData  *SolnData;

    void setSubdomainId(int sid)
    { subdomainId = sid; return; }

    int getSubdomainId()
    { return  subdomainId; }

    void setNumOfSubdomain(int i)
    { numOfSubdomain = i; return; }
    
    int getNumOfSubdomain()
    { return  numOfSubdomain; }
    
    void prepareElemData(const int nsize);


    /// ONLY FLUID ///
    vector<int>  nodeNumsFluid,nodeNumsPrevFluid;
    vector<int> nodeForAssyBCsFluid,nodeForAssyFluid;
    vector<int> globalDOFnumsFluid;
    vector<vector<vector<double>>> sub_x;
    
    SolutionData  *SolnDataFluid;


  private:
};