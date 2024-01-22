#include "FEM.h"

class ObservedParam
{
  public:
    ExportFile obs_export_file;
    
  public:
    int nx, ny, nz;
    double Lx, Ly, Lz;
    double dx, dy, dz;
    
    VDOUBLE3D u, v, w;

  public: 
    void inputVelocityMAP(const string file);
    void export_VelocityMAP(const string &file);
};

class CFDParam
{
  public:
    VDOUBLE3D u, v, w;
};

class CostFunction
{  
  public:
    double term1, term2, term3, total;
    VDOUBLE1D history;
    void sum(){
      total = term1 + term2 + term3;
    }
};

enum class EstimatedVelocity
{
  CENTER = 0,
  AVERAGE = 1
};


class VarDA :public FEM
{
  public: 
    VarDA(): FEM() {}
    ~VarDA();

    ObservedParam obs;
    CFDParam cfd;
    CostFunction costFunction;
    EstimatedVelocity estVel;
  
  public:
    int optMaxItr; int output_itr;
    
    double alpha_cf; 
    double beta1_cf; double beta2_cf;
    string controlBoundary;
    
    string mapFile;

    VDOUBLE1D u_adjoint; VDOUBLE1D v_adjoint;
    VDOUBLE1D w_adjoint; VDOUBLE1D p_adjoint;

    VDOUBLE3D ue, ve, we;
    VDOUBLE3D ue_edge, ve_edge, we_edge;

    VDOUBLE2D X;
    VDOUBLE2D grad;
    VDOUBLE2D grad_allNode;
    VDOUBLE2D feedbackForce;

    VINT2D control_elm_node;
    VINT1D control_node;
    VINT1D bdface_dir;

    int nx_in_voxel, ny_in_voxel, nz_in_voxel;
    double volume;

  public: 

    // INITIALIZE
    void initializeVarDA();
    void mainInitialize();
    void adjointInitialize();

    void readInputVarDA(); 
    void setControlBooundary();

    void resizeDAVariables();
    
    // OPTIMIZATION
    
    double ArmijoCriteria(const double fk,double loop);
    void calcOptimalCondition();
    
    // COST FUNCTION
    void calcCostFunction();
    void getAveragedVelocityFromCFD();
    void getCenterVelocityFromCFD();
    void gaussIntegral(VDOUBLE1D &N, VDOUBLE2D &dNdr, VDOUBLE2D &x_current, VDOUBLE2D &u_current, VDOUBLE1D &vel_local_ave_in_voxel, const double weight);
    void ReguralizationTerm2_inGaussIntegral(VDOUBLE1D &N, VDOUBLE2D &dNdr, VDOUBLE2D &x_current, double &value, const double weight, const int ic);
    void ReguralizationTerm3_inGaussIntegral(VDOUBLE1D &N, VDOUBLE2D &dNdr, VDOUBLE2D &x_current, double &value, const double weight, const int ic);
    
    // FEEDBACK FORCE
    void calcFeedbackForce();
    void calcInterpolatedFeedback(VDOUBLE1D &N, VDOUBLE2D &x_current, double (&feedbackGaussPoint)[3], const int ic);
    void calcEdgeVale();
    void feedback_inGaussIntegral(VDOUBLE3D &feedbackIntegral, VDOUBLE1D &N, VDOUBLE2D &dNdr, double (&feedbackGaussPoint)[3], VDOUBLE2D &x_current, const double weight, const int ic);
};

