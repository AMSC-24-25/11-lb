#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip> 


#include "LBM.hpp"



const int maxSteps = 10000; // Maximum number of time steps
const double Re = 10000;
const double u_lid = 0.5;

const int ITERATIONS_PER_FRAME = 10;
const int ITERATIONS_PER_PROGRESS_UPDATE = 100;

const int NX = 100; // Dimension in the x-direction
const int NY = 100; // Dimension in the y-direction

int main() {
    // Create the output file for velocity
    std::ofstream file("vel_data.txt");
    if (!file.is_open()) {
        std::cerr << "could not opene/create 'vel_data.txt'.\n";
        return 1;
    }
    file << NX << "\n" << NY << "\n";

    LBM lbm(NX, NY, u_lid, Re);

    auto startTime = std::chrono::high_resolution_clock::now();

    for (int n = 1; n <= maxSteps; n++) {
        lbm.evolution(); // System evolution

        //Every ITERATIONS_PER_FRAME steps, save velocity data
        if (n % ITERATIONS_PER_FRAME == 0) {
            for (int i = 0; i < lbm.NX; ++i) {
                for (int j = 0; j < lbm.NY; ++j) {
                    double vx = lbm.get_vel(i,j,0); 
                    double vy = lbm.get_vel(i,j,1);
                    double v = sqrt(vx*vx + vy*vy); 
                    file << v << "\n";
                }
            }
            
        }
        
        // Update the progress bar
        if (n % ITERATIONS_PER_PROGRESS_UPDATE == 0 || n == maxSteps - 1) {
            float progress = (static_cast<float>(n) / maxSteps);
            auto currentTime = std::chrono::high_resolution_clock::now();
            auto elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();
            
            double estimatedTotalTime = elapsedTime / progress;
            int remainingTime = estimatedTotalTime - elapsedTime;

            progress *= 100;
            std::cout << "\rProgress: " << std::fixed << std::setprecision(2) << progress << "% completed "
                      << "| Elapsed Time: " << elapsedTime << "s, "
                      << "Remaining Time (estimated): " << static_cast<int>(remainingTime) << "s"
                      << std::flush;
                      
        }
    }

    file.close();

    std::cout << "\n";

    return 0;
}
