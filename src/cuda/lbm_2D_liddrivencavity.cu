#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>

// Definition of simulation parameters
constexpr int Q = 9;      // Number of directions (D2Q9)
constexpr int D = 2;      // Spatial dimension

// Inline functions to translate indices to linear indices
__device__ inline int idx_density(int i, int j, int NY) {
    return i * NY + j;
}
__device__ inline int idx_field(int i, int j, int k, int NY) {
    return (i * NY + j) * Q + k;
}
__device__ inline int idx_velocity(int i, int j, int d,int NY) {
    return (i * NY + j) * D + d;
}

__device__ inline double feq_func(int k, double rho, double ux, double uy, const double* w) {
    // Definition of D2Q9 directions
    int cx[Q] = { 0, 1, 0, -1,  0, 1, -1, -1, 1 };
    int cy[Q] = { 0, 0, 1,  0, -1, 1,  1, -1, -1 };

    double eu = cx[k] * ux + cy[k] * uy;
    double uv = ux * ux + uy * uy;
    return w[k] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uv);
}

// Kernel for collision and streaming step (compute)
__global__ void kernel_compute(double* f, double* f2, double* rho, double* rho2, double* u, double* u2, double tau_f, const double* w, int nx, int ny) {
    // Identify the cell managed by this thread
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    // Check if the thread is inside the domain (skipping borders)
    if (i > 0 && i < nx - 1 && j > 0 && j < ny - 1) {
        int index = idx_density(i, j, ny);
        double local_rho = 0.0;
        double local_ux = 0.0;
        double local_uy = 0.0;

        int cx[Q] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
        int cy[Q] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };

        // Loop over the Q directions to perform collision and streaming
        for (int k = 0; k < Q; k++) {
            // Determine the index of the source node (backward streaming)
            int ip = i - cx[k];
            int jp = j - cy[k];

            // Compute the indices for velocity, density, and distribution field in the source node
            int idx_comp = idx_density(ip, jp,ny);
            int idx_field_comp = idx_field(ip, jp, k,ny);

            // Load density and velocity values from the source node
            double rho_comp = rho[idx_comp];
            double u_x_comp = u[idx_velocity(ip, jp, 0,ny)];
            double u_y_comp = u[idx_velocity(ip, jp, 1,ny)];

            // Compute equilibrium value (feq)
            double feq = feq_func(k, rho_comp, u_x_comp, u_y_comp, w);

            // Apply collision
            double new_val = f[idx_field_comp] + (feq - f[idx_field_comp]) / tau_f;

            // Update distribution field for the current cell
            int idx_field_current = idx_field(i, j, k,ny);
            f2[idx_field_current] = new_val;

            // Accumulate contributions to density and velocity
            local_rho += new_val;
            local_ux += cx[k] * new_val;
            local_uy += cy[k] * new_val;
        }

        // Update density and velocity (normalized)
        rho2[index] = local_rho;
        u2[idx_velocity(i, j, 0,ny)] = local_ux / local_rho;
        u2[idx_velocity(i, j, 1,ny)] = local_uy / local_rho;
    }
}



// Kernel to apply boundary conditions for both left/right and top/bottom borders
__global__ void kernel_apply_boundary(double* f, double* rho, double* u, double u_lid, int nx, int ny, const double* w) {
    // Compute (i,j) indices based on 2D thread block configuration
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    // Linear index for density and velocity (each cell has D components)
    int idx = i * ny + j;

    // LEFT AND RIGHT BORDERS: for j in [1, NY-2]
    if (j >= 1 && j < ny - 1) {
        if (i == 0 || i == nx - 1) {
            // Set zero velocity at lateral borders
            u[idx * D + 0] = 0.0;
            u[idx * D + 1] = 0.0;

            // Set density equal to the inner neighbor's
            if (i == 0)
                rho[idx] = rho[1 * ny + j];       // copy from (1,j)
            else // i == NX-1
                rho[idx] = rho[(nx - 2) * ny + j];  // copy from (NX-2,j)

            // Update distribution field for each direction
            for (int k = 0; k < Q; k++) {
                // For left border (i==0), the inner node is (1,j)
                // For right border (i==NX-1), the inner node is (NX-2,j)
                int i_comp = (i == 0) ? 1 : (nx - 2);
                int idx_comp = i_comp * ny + j;
                // Get density and velocity of inner node
                double rho_comp = rho[idx_comp];
                double ux_comp = u[idx_comp * D + 0];
                double uy_comp = u[idx_comp * D + 1];

                // Compute feq for the current node (with zero velocity) and for the neighbor
                double feq_current = feq_func(k, rho[idx], 0.0, 0.0, w);
                double feq_comp = feq_func(k, rho_comp, ux_comp, uy_comp, w);

                // Updates are performed using the formula:
                // f(boundary, j, k) = feq(boundary, j, k) + f(inner node, j, k) - feq(inner node, j, k)
                int idx_field_boundary = (i * ny + j) * Q + k;
                int idx_field_comp = (i_comp * ny + j) * Q + k;
                f[idx_field_boundary] = feq_current + f[idx_field_comp] - feq_comp;
            }
        }
    }

    // TOP AND BOTTOM BORDERS: for any i, j==0 or j==NY-1
    if (i < nx) {
        if (j == 0 || j == ny - 1) {
            if (j == 0) {
                // Bottom border: zero velocity
                u[idx * D + 0] = 0.0;
                u[idx * D + 1] = 0.0;
                // Copy density from inner neighbor (node (i,1))
                rho[idx] = rho[i * ny + 1];
            }
            else {
                // Top border: set velocity to lid value
                u[idx * D + 0] = u_lid;
                u[idx * D + 1] = 0.0;
                // Copy density from inner neighbor (node (i,NY-2))
                rho[idx] = rho[i * ny + (ny - 2)];
            }

            // Update distribution field for each direction
            for (int k = 0; k < Q; k++) {
                int idx_field = (i * ny + j) * Q + k;
                if (j == 0) {
                    int idx_comp = (i * ny + 1);
                    double rho_comp = rho[idx_comp];
                    double ux_comp = u[idx_comp * D + 0];
                    double uy_comp = u[idx_comp * D + 1];
                    double feq_boundary = feq_func(k, rho[idx], 0.0, 0.0, w);
                    double feq_comp = feq_func(k, rho_comp, ux_comp, uy_comp, w);
                    int idx_field_comp = (i * ny + 1) * Q + k;
                    f[idx_field] = feq_boundary + f[idx_field_comp] - feq_comp;
                }
                else {  // j == NY-1
                    int idx_comp = (i * ny + (ny - 2));
                    double rho_comp = rho[idx_comp];
                    double ux_comp = u[idx_comp * D + 0];
                    double uy_comp = u[idx_comp * D + 1];
                    double feq_boundary = feq_func(k, rho[idx], u_lid, 0.0, w);
                    double feq_comp = feq_func(k, rho_comp, ux_comp, uy_comp, w);
                    int idx_field_comp = (i * ny + (ny - 2)) * Q + k;
                    f[idx_field] = feq_boundary + f[idx_field_comp] - feq_comp;
                }
            }
        }
    }
}

__global__ void kernel_init(double* rho, double* rho2,
    double* u, double* u2,
    double* f, double* F,
    int nx, int ny, double rho0, double u_lid, const double* w) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < nx && j < ny) {
        int idx = idx_density(i, j,ny);

        // Initialize density and velocity fields
        rho[idx] = rho0;
        rho2[idx] = rho0;

        u[idx_velocity(i, j, 0,ny)] = 0.0;
        u[idx_velocity(i, j, 1,ny)] = 0.0;

        u2[idx_velocity(i, j, 0,ny)] = 0.0;
        u2[idx_velocity(i, j, 1,ny)] = 0.0;

        // Initialize distribution function with equilibrium values
        for (int k = 0; k < Q; k++) {
            double feq = feq_func(k, rho0, 0.0, 0.0, w);
            f[idx_field(i, j, k,ny)] = feq;
            F[idx_field(i, j, k,ny)] = feq;
        }
    }
}


constexpr int ITERATIONS_PER_PROGRESS_UPDATE = 100;
#include <chrono>
#include <iomanip> 

int main(int argc, char* argv[]) {

    if (argc != 7) {
        std::cerr << "Usage: " << argv[0]
                  << " <mesh_size> <time_steps> <reynolds> [output_dir]\n";
        return EXIT_FAILURE;
    }
    // Simulation parameters
    int NX = std::atoi(argv[1]);
    int NY = std::atoi(argv[2]);
    int MAX_STEPS = std::atoi(argv[3]);
    double Re = std::atof(argv[4]);
    int ITER_PER_FRAME = std::atoi(argv[5]);
    std::string out_dir = argv[6];
    const double u_lid = 0.5;
    int ITERATIONS_PER_FRAME = 200;
    
    const double dx = 1.0;
    double Lx = NY * dx;
    double nu = u_lid * Lx / Re;
    double tau_f = 3.0 * nu + 0.5;
    double rho0 = 1.0;

    // D2Q9 weights array
    double h_w[Q] = { 4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9,
                     1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36 };

    // Compute field sizes
    size_t size_density = NX * NY * sizeof(double);
    size_t size_field = NX * NY * Q * sizeof(double);
    size_t size_velocity = NX * NY * D * sizeof(double);

    // Device memory allocation
    double* d_f, * d_f2, * d_rho, * d_rho2, * d_u, * d_u2, * d_w;
    cudaMalloc(&d_f, size_field);
    cudaMalloc(&d_f2, size_field);
    cudaMalloc(&d_rho, size_density);
    cudaMalloc(&d_rho2, size_density);
    cudaMalloc(&d_u, size_velocity);
    cudaMalloc(&d_u2, size_velocity);
    cudaMalloc(&d_w, Q * sizeof(double));
    cudaMemcpy(d_w, h_w, Q * sizeof(double), cudaMemcpyHostToDevice);

    // Define block size and number
    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((NX + threadsPerBlock.x - 1) / threadsPerBlock.x,
        (NY + threadsPerBlock.y - 1) / threadsPerBlock.y);

    // Initialization of fields
    kernel_init << <numBlocks, threadsPerBlock >> > (d_rho, d_rho2, d_u, d_u2, d_f, d_f2, NX, NY, rho0, u_lid, d_w);
    cudaDeviceSynchronize();

    std::ofstream file_velocity(out_dir+"vel_data_cuda.txt");
    if (!file_velocity.is_open()) {
        std::cerr << "Error opening the file for velocity.\n";
        return 1;
    }
    file_velocity << NX << "\n" << NY << "\n";

    auto startTime = std::chrono::high_resolution_clock::now();

    for (int step = 1; step <= MAX_STEPS; step++) {
        // Collision and streaming step
        kernel_compute << <numBlocks, threadsPerBlock >> > (d_f, d_f2, d_rho, d_rho2, d_u, d_u2,
            tau_f, d_w, NX, NY);
        cudaDeviceSynchronize();

        // Apply boundary conditions
        kernel_apply_boundary << <numBlocks, threadsPerBlock >> > (d_f, d_rho, d_u,
            u_lid, NX, NY, d_w);
        cudaDeviceSynchronize();

        if (step % ITERATIONS_PER_FRAME == 0 || step == 1 || step == MAX_STEPS) {
            // Transfer results to host for saving to file 
            double* temp_vel = new double[NX * NY * D];
            cudaMemcpy(temp_vel, d_u, size_velocity, cudaMemcpyDeviceToHost);

            for (int j = 0; j < NY; ++j) {
                for (int i = 0; i < NX; ++i) {
                    double vx = temp_vel[(i * NY + j) * D + 0];
                    double vy = temp_vel[(i * NY + j) * D + 1];
                    double v = sqrt(vx * vx + vy * vy);
                    file_velocity << v << "\n";
                }
            }
            delete[] temp_vel;
        }

        // Swap pointers to prepare for next step
        double* tmp;
        tmp = d_f; d_f = d_f2; d_f2 = tmp;
        tmp = d_rho; d_rho = d_rho2; d_rho2 = tmp;
        tmp = d_u; d_u = d_u2; d_u2 = tmp;

        // Update the progress bar
        if (step % ITERATIONS_PER_PROGRESS_UPDATE == 0 || step == MAX_STEPS) {
            float progress = (static_cast<float>(step) / MAX_STEPS);
            auto currentTime = std::chrono::high_resolution_clock::now();
            auto elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();

            double estimatedTotalTime = elapsedTime / progress;
            int remainingTime = estimatedTotalTime - elapsedTime;

            progress *= 100;
            std::cout << "\rProgress: " << std::fixed << std::setprecision(2) << progress << "% completed "
                << "| Elapsed Time: " << elapsedTime << "s, "
                << "Remaining Time (estimated): " << static_cast<int>(remainingTime) << "s" << "        "
                << std::flush;

        }
    }

    file_velocity.close();

    // Free GPU and host memory
    cudaFree(d_f);
    cudaFree(d_f2);
    cudaFree(d_rho);
    cudaFree(d_rho2);
    cudaFree(d_u);
    cudaFree(d_u2);
    cudaFree(d_w);

    return 0;
}
