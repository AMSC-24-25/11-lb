#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <filesystem>
#include <omp.h>

namespace fs = std::filesystem;


#include "3DLBM.hpp"

const int nx = 40;
const int ny = 40;
const int nz = 40;
const int Steps = 1000;             
const int ITERATIONS_PER_FRAME = 20;  
const int ITERATIONS_PER_PROGRESS_UPDATE = 10;
const int Re = 1000; 
double u_lid = 0.1;


void writeVTKFrame(const std::string &filename, int iteration,
    const TDLBM& lbm) {
std::ofstream vtkFile(filename);
if (!vtkFile.is_open()) {
std::cerr << "Errore nell'apertura del file " << filename << "\n";
exit(EXIT_FAILURE);
}

vtkFile << "# vtk DataFile Version 2.0\n";
vtkFile << "Driven Cavity Simulation, Iteration " << iteration << "\n";
vtkFile << "ASCII\n";
vtkFile << "DATASET STRUCTURED_POINTS\n";
vtkFile << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
vtkFile << "ORIGIN 0 0 0\n";
vtkFile << "SPACING " << lbm.get_dx() << " " << lbm.get_dx()<< " " << lbm.get_dx() << "\n";
vtkFile << "POINT_DATA " << (nx * ny * nz) << "\n";

vtkFile << "SCALARS rho double 1\n";
vtkFile << "LOOKUP_TABLE default\n";
for (int z = 0; z < nz; ++z) {
for (int y = 0; y < ny; ++y) {
for (int x = 0; x < nx; ++x) {
 int index = lbm.idx(x, y, z);
 vtkFile << std::setprecision(8) << lbm.get_rho()[index] << "\n";
}
}
}

vtkFile << "VECTORS velocity double\n";
for (int z = 0; z < nz; ++z) {
for (int y = 0; y < ny; ++y) {
for (int x = 0; x < nx; ++x) {
    int index = lbm.idx(x, y, z);
 vtkFile << lbm.get_u()[index].x << " " << lbm.get_u()[index].y << " " << lbm.get_u()[index].z << "\n";
}
}
}
vtkFile.close();
}

int main(){
    TDLBM Cavity(nx, ny, nz, u_lid, Re);
    std::cout << "### Lattice-Boltzmann D3Q19 per Driven Cavity (TRT) ###\n";
    std::cout << "Numero di Reynolds target: " << Re << "\n";
    std::cout << "Viscosità calcolata: " << Cavity.get_nu() << "\n";
    std::cout << "τ (tau) = " << Cavity.get_tau() << "\n";

    

    
    std::vector<std::pair<int, std::string>> frames;
    
    std::string outputDir = "output_" + std::to_string(Re) + "_" +
                            std::to_string(nx) + "x" + std::to_string(ny) + "x" + std::to_string(nz);

    try {
        fs::create_directory(outputDir);
    } catch (const fs::filesystem_error& err) {
        std::cerr << "Errore nella creazione della directory: " << err.what() << "\n";
        return 1;
    }

    
    auto startTime = std::chrono::high_resolution_clock::now();
    int frameCounter = 0;
    for (int n = 1; n <= Steps; ++n) {
        Cavity.evolution();

        if (n == 1 || n % ITERATIONS_PER_FRAME == 0) {
            std::string filename = outputDir + "/frame_" + std::to_string(frameCounter) + ".vtk";
            writeVTKFrame(filename, n, Cavity);
            frames.push_back({n, filename});
            ++frameCounter;
        }

        if (n % ITERATIONS_PER_PROGRESS_UPDATE == 0 || n == Steps) {
            float progress = static_cast<float>(n) / Steps;
            auto currentTime = std::chrono::high_resolution_clock::now();
            auto elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();
            double estimatedTotalTime = (elapsedTime / progress);
            int remainingTime = static_cast<int>(estimatedTotalTime - elapsedTime);
            progress *= 100;
            std::cout << "\rProgress: " << std::fixed << std::setprecision(2) << progress 
                      << "% | Elapsed: " << elapsedTime << "s, Remaining (est.): " 
                      << remainingTime << "s" << std::flush;
        }

        
    }

    std::cout << "\nSimulazione completata.\n";
    return 0;



}