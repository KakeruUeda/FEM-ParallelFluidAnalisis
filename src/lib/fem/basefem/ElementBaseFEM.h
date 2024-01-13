#include "solutionData.h"

class SubProperty{
 public:
  VDOUBLE1D sub_elm_sdf;
  VDOUBLE1D sub_elm_x;
  VDOUBLE1D sub_gx;
  VDOUBLE1D sub_gy;
  VDOUBLE1D sub_gz;
  VDOUBLE1D sub_weight;
};

class ElementBaseFEM{
  public:
    int subdomainId,numOfDofsNode,numOfNodeInElm,numOfSubdomain;
    //VINT1D  nodeNums,nodeNumsPrev,nodeForAssyBCs,nodeForAssy,globalDOFnums;

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
    VINT1D nodeNumsFluid,nodeNumsPrevFluid;
    VINT1D nodeForAssyBCsFluid,nodeForAssyFluid;
    VINT1D globalDOFnumsFluid;
    vector<vector<vector<double>>> sub_x;
    
    SolutionData  *SolnDataFluid;

    vector<vector<SubProperty>> sub_elm_node;
    vector<SubProperty> sub_elm;

  private:
};


