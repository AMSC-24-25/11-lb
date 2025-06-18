#include "LDRIVEN.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <iomanip>
#ifdef USE_OPENMP
  #include <omp.h>
#endif

LDRIVEN::LDRIVEN(unsigned int nx, unsigned int ny, double u_lid, double Re,std::string outdir) : NX(nx), NY(ny), u_lid(u_lid), Re(Re), outdir(outdir) {
    dx = 1.0;         // Spatial step
    dy = 1.0;         // Spatial step
    Lx = dx * double(NY); // Domain length in y
    Ly = dy * double(NX); // Domain length in x
    dt = dx;          // Time step
    c = dx / dt;      // Speed of sound (in lattice units)
    rho0 = 1.0;       // Initial density
    nu = u_lid * Lx / Re; // Kinematic viscosity
    tau_f = 3.0 * nu + 0.5; // Relaxation time

    rho = (double*)malloc(sizeof(double)*(NX)*(NY));
    rho2 = (double*)malloc(sizeof(double)*(NX)*(NY));
    u = (double*)malloc(sizeof(double)*(NX)*(NY)*D);
    u2 = (double*)malloc(sizeof(double)*(NX)*(NY)*D);
    f = (double*)malloc(sizeof(double)*(NX)*(NY)*Q);
    F = (double*)malloc(sizeof(double)*(NX)*(NY)*Q);



    // Initialize variables in the domain
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for(int k = 0; k < D; k++) {
                velocity(i,j,k) = 0.0; // Initial velocity
                velocity_2(i,j,k) = 0.0;
            }

            density(i,j) = rho0; // Initial density
            density_2(i,j) = rho0;
			// velocity(i, NY-1, 0) = u_lid;  // Velocity imposed on the upper boundary
            // velocity_2(i, NY-1, 0) = u_lid;

            for (int k = 0; k < Q; k++) {
                field(i,j,k) = feq(k, i, j); // Initial equilibrium function
            }
        }
    }
    
}

LDRIVEN::~LDRIVEN() {
    // freeing the allocated memory
    free(rho);
    free(rho2);
    free(u);
    free(u2);
    free(f);
    free(F);
}

double LDRIVEN::feq(unsigned int k, unsigned int x, unsigned int y) {  
    double eu = direction(k,0) * velocity(x,y,0) + direction(k,1) * velocity(x,y,1);  
    double uv = velocity(x,y,0) * velocity(x,y,0) + velocity(x,y,1) * velocity(x,y,1);

    // Compute equilibrium distribution function value
    return w[k] * density(x,y) * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv);
}

void LDRIVEN::evolution() {
    compute();      // Collision and Update macroscopic density and velocities
    apply_boundary_conditions(); // Apply boundary conditions
}


void LDRIVEN::simulate(unsigned int iterations,int iter_per_frame) {
    maxSteps=iterations;
    ITERATIONS_PER_FRAME=iter_per_frame;
    // Apertura file output
    std::ofstream file(outdir+"/vel_data.txt");
    if (!file.is_open()) {
        std::cerr << "Could not open/create 'vel_data.txt'.\n";
        return;
    }
    file << NX << "\n" << NY << "\n";

    auto startTime = std::chrono::high_resolution_clock::now();
    
    std::ofstream param_file(outdir+"/parameters.txt", std::ios::trunc);    
    if (param_file.is_open()) {
        param_file << "# ---- LBM Simulation Parameters ----\n\n";
        param_file << "# Domain Size\n";
        param_file << "NX = " << NX << "\n";
        param_file << "NY = " << NY << "\n\n";

        param_file << "# Physical Parameters\n";
        param_file << "Re = " << Re << "\n";
        param_file << "Lid velocity = " << u_lid << "\n";

        param_file << "\n# Derived Parameters\n";
        param_file << "nu (kinematic viscosity) = " << nu << "\n";
        param_file << "tau (relaxation time) = " << tau_f << "\n";

        param_file << "\n# Output files:\n";
        param_file << "vel_data.txt (velocity field)\n";
    } 
    else {
        std::cerr << "Could not open/create 'parameters.txt'\n";
    }
    std::cout << "\nParameters saved in 'parameters.txt'" << std::endl;
    param_file.close();

    for(int iter = 0; iter <= maxSteps; iter++) {
        compute();      // Collision and Update macroscopic density and velocities
        apply_boundary_conditions(); // Apply boundary conditions
        if(iter % ITERATIONS_PER_FRAME == 0) {
                for (int j = 0; j < NY; ++j) {
                    for (int i = 0; i < NX; ++i) {
                        double vx = u[((NY)*i+j)*D + 0];
                        double vy = u[((NY)*i+j)*D + 1];
                        double v  = std::sqrt(vx * vx + vy * vy);
                        file << v << "\n";
                    }
                }
            }

            if(iter % ITERATIONS_PER_PROGRESS_UPDATE == 0 || iter == maxSteps - 1) {
            float progress = (static_cast<float>(iter) / maxSteps);
            auto currentTime = std::chrono::high_resolution_clock::now();
            auto elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();

            std::cout << "\r" << std::string(100, ' ') << "\r"; // clean row
            
            if (progress > 0.001f) { // evita divisioni assurde
                double estimatedTotalTime = elapsedTime / progress;
                int remainingTime = static_cast<int>(estimatedTotalTime - elapsedTime);
                
                progress *= 100;
                std::cout << "\rProgress: " << std::fixed << std::setprecision(2) << progress << "% "<<"completed "
                          << "| Elapsed Time: " << elapsedTime << "s, "
                          << "Remaining Time (estimated): " << remainingTime << "s"
                          << std::flush;
            } else {
                std::cout << "\rProgress: " << std::fixed << std::setprecision(2) << progress * 100 << "% "<<"completed "
                          << "| Elapsed Time: " << elapsedTime << "s, "
                          << "Remaining Time (estimated): calculating..."
                          << std::flush;
            }      
    }
}
file.close();
    std::cout << "\nSimulation completed.\n";
}


void LDRIVEN::compute() {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            density_2(i,j) = 0; // Initialize updated density

            velocity_2(i,j,0) = 0; // Initialize updated velocity x-component
            velocity_2(i,j,1) = 0; // Initialize updated velocity y-component

            for (int k = 0; k < Q; k++) {
                int ip = (i - direction(k,0) + NX) % NX; // Node from which the contribution comes
                int jp = (j - direction(k,1) + NY) % NY;

                // Compute collision step
                field_2(i,j,k) = field(ip,jp,k) +(feq(k, ip, jp) - field(ip,jp,k)) / tau_f; // Collision

                density_2(i,j) += field_2(i,j,k);
            
                velocity_2(i,j,0) += direction(k,0) * field_2(i,j,k);
                velocity_2(i,j,1) += direction(k,1) * field_2(i,j,k);
            }

            // Normalize velocity components
            velocity_2(i,j,0) /= density_2(i,j);
            velocity_2(i,j,1) /= density_2(i,j);
        }
    }

    // Swap pointers to update fields
    double *t;

    t = u;
    u = u2;
    u2 = t;

    t = rho;
    rho = rho2;
    rho2 = t;

    t = f;
    f = F;
    F = t;
}


void LDRIVEN::apply_boundary_conditions() {
    #pragma omp parallel
    {
        // Left e Right boundaries
        #pragma omp for nowait
        for (int j = 1; j < NY - 1; j++) {
            // Left boundary (x = 0)
            velocity(0, j, 0) = 0.0;
            velocity(0, j, 1) = 0.0;
            velocity(NX-1, j, 0) = 0.0;
            velocity(NX-1, j, 1) = 0.0;

            density(0,j) = density(1,j);
            density(NX-1,j) = density(NX-2,j);
            for (int k = 0; k < Q; k++) {
                field(0, j, k) = feq(k, 0, j) + field(1, j, k) - feq(k, 1, j);
                field(NX - 1, j, k) = feq(k, NX - 1, j) + field(NX - 2, j, k) - feq(k, NX - 2, j);
            }
        }

        // Top e Bottom boundaries
        #pragma omp for
        for (int i = 0; i < NX; i++) {
            // Top boundary (y = NY-1)
            velocity(i, NY - 1, 0) = u_lid;
            velocity(i, NY - 1, 1) = 0.0;
            velocity(i, 0, 0) = 0.0;
            velocity(i, 0, 1) = 0.0;

            density(i,NY-1) = density(i,NY - 2);
            density(i,0) = density(i,1);
            for (int k = 0; k < Q; k++) {
                field(i, NY - 1, k) = feq(k, i, NY - 1) + field(i, NY - 2, k) - feq(k, i, NY - 2);
                field(i, 0, k) = feq(k, i, 0) + field(i, 1, k) - feq(k, i, 1);
            }
        }

    }
}
