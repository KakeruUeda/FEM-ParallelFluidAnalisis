#include "FEM.h"

class ObservedParam
{
  public:
    int nx, ny, nz;
    double Lx, Ly, Lz;
    double dx, dy, dz;
    
    VDOUBLE3D u, v, w;
    
    void inputVelocityMAP(const string file);
    void export_VelocityMAP(const string &file);
};

class CFDParam
{
  public:
    VDOUBLE3D u, v, w;
};


class VarDA :public FEM
{
  public: 
    VarDA();
    ~VarDA();

    ObservedParam obs;
    CFDParam cfd;
  
  public:
    int optMaxItr; int output_iter;
    
    double alpha_cf; 
    double beta1_cf; double beta2_cf;
    string controlBoundary;
    
    string mapFile;

    VDOUBLE1D u_adjoint; VDOUBLE1D v_adjoint;
    VDOUBLE1D w_adjoint; VDOUBLE1D p_adjoint;

    VDOUBLE3D ue, ve, we;

    VDOUBLE2D X;
    VDOUBLE2D grad;
    VDOUBLE2D grad_allNode;
    VDOUBLE2D feedbackForce;

    VINT2D control_elm_node;
    VINT1D control_node;

  public: 
    void initializeVarDA();
    void mainInitialize();
    void adjointInitialize();

    void readInputVarDA(); 
    void setControlBooundary();
    
    void calcFeedbackForce();
    double ArmijoCriteria(const double fk,double loop);
    void calcOptimalCondition();
    void calcCostFunction();
};

