#include "VariationalDataAssimilation.h"

void ObservedParam::inputVelocityMAP(const std::string file)
{
  FILE *fp;
  if ((fp = fopen(file.c_str(), "r")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }

  fscanf(fp, "%d %d\n", &nx, &ny, &nz);
  fscanf(fp, "%lf %lf\n", &Lx, &Ly, &Lz);
  
  dx = Lx/nx;
  dy = Ly/ny;
  dz = Lz/nz;

  u.resize(nz, VDOUBLE2D(ny, VDOUBLE1D(nx, 0e0)));
  v.resize(nz, VDOUBLE2D(ny, VDOUBLE1D(nx, 0e0)));
  w.resize(nz, VDOUBLE2D(ny, VDOUBLE1D(nx, 0e0)));

  for(int k=0; k<nz; k++){
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
        int tmp1, tmp2, tmp3;
        fscanf(fp, "%d %d %lf %lf", &tmp1, &tmp2, &tmp3, &u[j][i], &v[j][i], &w[j][i]);
      }
    }
  } 

  //for debug
   string vti = "test.vti";
   //export_obs(vti);
}
