#include "VTK.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>

VTK::VTK(unsigned int nx, unsigned int ny, unsigned int nz, double dx)
    : nx_(nx), ny_(ny), nz_(nz), dx_(dx) {}

void VTK::writeVTKFrame(const std::string& filename,
                        int iteration,
                        const TDLBM& lbm) const
{
    std::ofstream vtkFile(filename);
    if (!vtkFile.is_open()) {
        std::cerr << "Error opening file " << filename << "\n";
        exit(EXIT_FAILURE);
    }

    vtkFile << "# vtk DataFile Version 2.0\n";
    vtkFile << "Driven Cavity Simulation, Iteration " << iteration << "\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_POINTS\n";
    vtkFile << "DIMENSIONS " << nx_ << " " << ny_ << " " << nz_ << "\n";
    vtkFile << "ORIGIN 0 0 0\n";
    vtkFile << "SPACING " << lbm.get_dx() << " "
            << lbm.get_dx() << " " << lbm.get_dx() << "\n";
    vtkFile << "POINT_DATA " << (nx_ * ny_ * nz_) << "\n";

    vtkFile << "SCALARS rho double 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int z = 0; z < nz_; ++z)
      for (int y = 0; y < ny_; ++y)
        for (int x = 0; x < nx_; ++x)
          vtkFile << std::setprecision(8)
                  << lbm.get_rho()[lbm.idx(x,y,z)] << "\n";

    vtkFile << "VECTORS velocity double\n";
    for (int z = 0; z < nz_; ++z)
      for (int y = 0; y < ny_; ++y)
        for (int x = 0; x < nx_; ++x) {
          auto v = lbm.get_u()[lbm.idx(x,y,z)];
          vtkFile << v.x << " " << v.y << " " << v.z << "\n";
        }

    vtkFile.close();
}
