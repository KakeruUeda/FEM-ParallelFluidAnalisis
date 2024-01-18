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

    //VDOUBLE4D vel_ave_voxel;
    VINT2D numOfObsVoxels;

    double Lx_obs, Ly_obs, Lz_obs;
    int nx_obs, ny_obs, nz_obs;
    double dx_obs, dy_obs, dz_obs;

    int nx_in_voxel, ny_in_voxel, nz_in_voxel;
    int numOfElmInVoxel, numOfNodeInVoxel;
   
    double x0, y0, z0;
    int voxelDivision;

  public:
    
    void prepareForDataAssimilation(FEM &fem, const int xNumOfVoxels, const int yNumOfVoxels, const int zNumOfVoxels);
    void readDAParam();
    void makeVoxelAveragedVelocity(FEM &fem, VDOUBLE4D &vel_ave_voxel);
    void gaussIntegral(FEM &fem, VDOUBLE1D &N, VDOUBLE2D &dNdr, VDOUBLE2D &x_current, VDOUBLE2D &u_current, VDOUBLE1D &vel_local_ave_in_voxel, const double weight);
    void writeVelocityToFile(FEM &fem, VDOUBLE4D &vel_ave_voxel);

};