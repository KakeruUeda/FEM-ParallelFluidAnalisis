#include "FEM.h"
#include "define.h"

using namespace std;

class ElementPropertyDA
{
  public:
    VINT1D node;
};

class PostProcess
{
   public:
     ElementPropertyDA **elmDA;
     ExportFile export_file_post;
     TextParser tp_post;

   public:
    
    // DA OBSERVED VELOCITY
    VINT2D numOfObsVoxels;

    double Lx_obs, Ly_obs, Lz_obs;
    int  nx_obs, ny_obs, nz_obs;
    double dx_obs, dy_obs, dz_obs;

    int nx_in_voxel, ny_in_voxel, nz_in_voxel;
    int numOfElmInVoxel, numOfNodeInVoxel;
   
    double x0, y0, z0;
    int voxelDivision;

    // DA REFERENCE VELOCITY && DOMAIN
    VDOUBLE3D sdf_opt;
    VDOUBLE3D phiVOF_opt;
    VDOUBLE4D x_opt;
    VDOUBLE4D vel_ref_opt;

    int nx_opt, ny_opt, nz_opt;
    double dx_opt, dy_opt, dz_opt;

  public:
    
    void prepareForDataAssimilation(FEM &fem, const int xNumOfVoxels, const int yNumOfVoxels, const int zNumOfVoxels);
    void readDAParam();
    void validateDADomain(FEM &fem, const int xNumOfVoxels, const int yNumOfVoxels, const int zNumOfVoxels);
    void makeVoxelAveragedVelocity(FEM &fem, VDOUBLE4D &vel_ave_voxel);
    void gaussIntegral(FEM &fem, VDOUBLE1D &N, VDOUBLE2D &dNdr, VDOUBLE2D &x_current, VDOUBLE2D &u_current, VDOUBLE1D &vel_local_ave_in_voxel, const double weight);
    void writeAveragedVelocityToFile(FEM &fem, VDOUBLE4D &vel_ave_voxel);
    void extractDomain(FEM &fem);
   
    double sdfToChi(FEM &fem, const double &sdf);
    void writeDADomainToFile(FEM &fem);
};