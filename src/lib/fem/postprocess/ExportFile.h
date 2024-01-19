#include <omp.h>
#include <string>
#include <iostream>
#include <cmath>
#include "define.h"

class ExportFile
{
  public:
    
    void export_vti_metis(const string &file, VINT1D &node, VINT1D &element, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
    void export_vti_node(const string &file, VDOUBLE1D &node, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
    void export_vti_elm(const string &file, VDOUBLE1D &element, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
    void export_vti_node_xyz(const string &file, VDOUBLE3D &node, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
    void export_vti_elm_xyz(const string &file, VDOUBLE3D &element, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
    void export_vti_domain(const string &file, VDOUBLE1D &node, VDOUBLE1D &phi, VDOUBLE1D &phiEX, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
    void export_vti_result(const string &file, VDOUBLE1D &u, VDOUBLE1D &v, VDOUBLE1D &w, VDOUBLE1D &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
    void export_vti_result_2D(const string &file, VDOUBLE1D &u, VDOUBLE1D &v, VDOUBLE1D &p, const int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
    void export_vti_velocity_cell(const string &file, VDOUBLE4D &vel, int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
    void export_vti_velocity_node(const string &file, VDOUBLE4D &vel, int nx, const int ny, const int nz, const double dx, const double dy, const double dz);
};