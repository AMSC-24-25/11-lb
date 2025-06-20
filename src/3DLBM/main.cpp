#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <filesystem>
#include "include/3DLBM.hpp"
#include "include/VTK.hpp"

#ifdef USE_OPENMP
  #include <omp.h>
#endif

namespace fs = std::filesystem;

const int ITERATIONS_PER_PROGRESS_UPDATE = 10;
int nx, ny, nz;

int main(int argc, char* argv[]) {
    if (argc !=7) {
        std::cerr << "Usage: " << argv[0]
                  << " <mesh_size> <time_steps> <reynolds> \n";
        return EXIT_FAILURE;
    }

    int mesh  = std::stoi(argv[1]);
    int Steps = std::stoi(argv[2]);
    int Re    = std::stoi(argv[3]);
    int ITERATIONS_PER_FRAME=std::stoi(argv[5]);
    double u_lid = std::stod(argv[6]);;
    nx = mesh; ny = mesh; nz = mesh;

    #ifdef USE_OPENMP
    #pragma omp parallel
    {
        int n_threads = omp_get_num_threads();

        #pragma omp single
        std::cout << ">>> OpenMP active with " 
                  << n_threads 
                  << " thread" <<"\n";

    }
    #else
    std::cout << "OpenMP not active\n";
    #endif

    TDLBM Cavity(nx, ny, nz, u_lid, Re);
    VTK writer(nx, ny, nz, Cavity.get_dx());

    std::cout << "### Lattice-Boltzmann D3Q19 for Driven Cavity (TRT) ###"<< "\n";
    std::cout << "Target Reynolds number: " << Re << "\n" ;
    std::cout << "Calculated viscosity: " << Cavity.get_nu() << "\n";
    std::cout << "τ (tau) = " << Cavity.get_tau() << "\n";

    if(Cavity.get_tau()<0.6){
        std::cout << "\n-----WARNING: TAU IS LOW-----" << "\n"; 
        std::cout << "SIMULATION CAN DIVERGE\n";
    }

    std::vector<std::pair<int, std::string>> frames;
    std::string outputDir;
    std::string ArgDir=argv[4];

   
       
        outputDir = ArgDir + "/output_" + std::to_string(Re) + "_"
                    + std::to_string(nx) + "x" + std::to_string(ny) + "x" + std::to_string(nz);
    

    try {
        fs::create_directory(outputDir);
    } catch (const fs::filesystem_error& err) {
        std::cerr << "Error creating directory: " << err.what() << "\n";
        return EXIT_FAILURE;
    }

    auto startTime = std::chrono::high_resolution_clock::now();
    int frameCounter = 0;

    for (int n = 1; n <= Steps; ++n) {
        Cavity.simulate();

        if (n == 1 || n % ITERATIONS_PER_FRAME == 0) {
            std::string filename = outputDir + "/frame_" + std::to_string(frameCounter) + ".vtk";
            writer.writeVTKFrame(filename, n, Cavity);
            frames.emplace_back(n, filename);
            ++frameCounter;
        }

        if (n % ITERATIONS_PER_PROGRESS_UPDATE == 0 || n == Steps) {
            float progress = static_cast<float>(n) / Steps;
            auto currentTime = std::chrono::high_resolution_clock::now();
            auto elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();
            double estimatedTotalTime = elapsedTime / progress;
            int remainingTime = static_cast<int>(estimatedTotalTime - elapsedTime);
            progress *= 100.0f;
            std::cout << "\rProgress: "
                      << std::fixed << std::setprecision(2) << progress
                      << "% | Elapsed: " << elapsedTime << "s, Remaining (est.): "
                      << remainingTime << "s" << std::flush;
        }
    }

    std::cout << "\nSimulation completed.\n";
    return EXIT_SUCCESS;
}
