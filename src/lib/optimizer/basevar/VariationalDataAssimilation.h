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
    PetscSolver *adjointSolverPetscFluid;

  
  public:
    int optMaxItr; int output_itr;
    
    double alpha_cf; 
    double beta1_cf; double beta2_cf;
    string controlBoundary;
    
    string mapFile;

    PetscInt numOfDofsAdjointLocalFluid, numOfDofsAdjointGlobalFluid;
    int row_start_adjoint, row_end_adjoint;
    VINT1D  assyForSolnAdjointFluid;
    VINT2D  nodeDofArrayBCsPrevAdjointFluid, nodeDofArrayPrevAdjointFluid, nodeDofArrayBCsAdjointFluid, nodeDofArrayAdjointFluid;
    VBOOL2D  nodeTypePrevAdjointFluid, nodeTypeAdjointFluid;
    
    VINT1D numOfDofsNodePrevAdjointFluid;
    VINT2D numOfDofsNodeInElementPrevAdjointFluid;
    VINT1D numOfDofsNodeAdjointFluid;
    VINT2D numOfDofsNodeInElementAdjointFluid;

    VINT1D nodeDofsMapPrevFluid;
    VINT1D nodeDofsMapFluid;

    VDOUBLE1D DirichletBCsAdjointFluid;

    VDOUBLE1D u_adjoint; VDOUBLE1D v_adjoint;
    VDOUBLE1D w_adjoint; VDOUBLE1D p_adjoint;
    VDOUBLE1D lambda_u; VDOUBLE1D lambda_v; VDOUBLE1D lambda_w;

    VDOUBLE1D u_adjoint_fluid; VDOUBLE1D v_adjoint_fluid;
    VDOUBLE1D w_adjoint_fluid; VDOUBLE1D p_adjoint_fluid;
    VDOUBLE1D lambda_u_fluid; VDOUBLE1D lambda_v_fluid; VDOUBLE1D lambda_w_fluid;
    
    VDOUBLE3D ue, ve, we;
    VDOUBLE3D ue_edge, ve_edge, we_edge;

    VDOUBLE2D X;
    VDOUBLE2D grad;
    VDOUBLE2D grad_allNode;
    VDOUBLE2D feedbackForce;
    VDOUBLE2D feedbackForceFluid;
    
    VBOOL1D control_elm_prev_type;
    VBOOL1D control_node_prev_type;
    VBOOL2D control_elm_node_prev_type;
    VBOOL1D control_elm_type;
    VBOOL1D control_node_type;
    VINT2D control_elm_node;
    VINT1D control_node;

    VBOOL1D control_elm_prev_type_fluid;
    VBOOL1D control_node_prev_type_fluid;
    VBOOL2D control_elm_node_prev_type_fluid;
    VBOOL1D control_elm_type_fluid;
    VBOOL1D control_node_type_fluid;
    VINT2D  control_elm_node_fluid;
    VINT1D  control_node_fluid;
    
    VINT1D bdface_dir;
    VINT1D bdnodeInElm;

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
    void setPetscAdjointSolver();
    
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
    
    // ADJOINT EQUATION
    void adjointSteadyNavierStokes(VDOUBLE2D &externalForce);

    void assignBCsAdjoint();
    void applyBCsAdjoint();
    void AdjointMatAssySNS(const int ic, MatrixXd &Klocal, VectorXd &Flocal, VDOUBLE2D &externalForce);
    void Darcy_AdjointMatAssySNS(const int ic, MatrixXd &Klocal, VectorXd &Flocal, VDOUBLE2D &externalForce);

    void AdjointBdMatAssySNS(const int ic, MatrixXd &Klocal, VectorXd &Flocal, VDOUBLE2D &externalForce);
    void Darcy_AdjointBdMatAssySNS(const int ic, MatrixXd &Klocal, VectorXd &Flocal, VDOUBLE2D &externalForce);
   
    void setAdjointMatAndVecZero();
    void calcBoundaryIntegral(const int ic, MatrixXd &Klocal, VectorXd &Flocal);

    void scatterPysicalvVariables();
};

