#include "Lattice.hpp"
#include "Masks.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip> 
#include <cmath>
#include <string>
#ifdef USE_OPENMP
  #include <omp.h>
#endif

std::vector <int> evaluateBoundary(const std::vector<int>& indices, const Matrix<bool> &obstacleMatrix)
{
    const int &col = indices.at(0);
    const int &row = indices.at(1);
    const int &NX = obstacleMatrix.shape().at(0);
    const int &NY = obstacleMatrix.shape().at(1);
    std::vector<int> boundary_here(4);

    // Evaluate boundary array, considering near objects
    // Horizontal 
    if( col>0 && obstacleMatrix.getCopy(col-1,row))
    {
        boundary_here.at(0) = -1;
    }
    else if(col<(NX-1) && obstacleMatrix.getCopy(col+1,row)) 
    {
        boundary_here.at(0) = 1;
    }
    else 
    {
        boundary_here.at(0) = 0;
    }

    // Vertical 
    if (row>0 && obstacleMatrix.getCopy(col,row-1))
    {
        boundary_here.at(1) = -1;
    }
    else if (row<(NY-1) && obstacleMatrix.getCopy(col,row+1))
    {
        boundary_here.at(1) = 1;
    }
    else 
    {
        boundary_here.at(1) = 0;
    }

    // First Diag
    if (col>0 && row>0 && obstacleMatrix.getCopy(col-1,row-1))
    {
        boundary_here.at(2) = -1;
    }
    else if (col<(NX-1) && row<(NY-1) && obstacleMatrix.getCopy(col+1, row+1))
    {
        boundary_here.at(2) = 1;
    }
    else 
    {
        boundary_here.at(2) = 0;
    }

    // Second Diag
    if (col>0 && row<(NY-1) && obstacleMatrix.getCopy(col-1,row+1))
    {
        boundary_here.at(3) = 1;
    }
    else if (col<(NX-1) && row>0 && obstacleMatrix.getCopy(col+1, row-1))
    {
        boundary_here.at(3) = -1;
    }
    else 
    {
        boundary_here.at(3) = 0;
    }  
    return boundary_here;
}

Lattice::Lattice(unsigned int nx, unsigned int ny, double u_lid, double Re, std::string outdir): NX(nx), NY(ny), u_lid(u_lid), Re(Re), outdir(outdir)
{
    /*  TO SET MANUALLY  */
    // Hard coded variables

    ITERATIONS_PER_PROGRESS_UPDATE = 10;
    boundary_velocity.resize(2);
    boundary_velocity.at(0) = 0.1; // Vx
    boundary_velocity.at(1) = 0.0; // Vy

    // Calculating other parameters
    const double nu = 2.0/3.0*std::max(NX,NY)*boundary_velocity.at(0)/Re;
    const double tau = 0.5 + nu * 3.0;
    if (tau<0.6)
    {
        std::cout << "\n-----WARNING: TAU IS LOW-----" << std::endl; 
        std::cout << "Tau: " << tau << std::endl;
    }
    sigma = 10.0*std::min(NY,NX);
    const double lambda_trt = 1.0 / 4.0;
    const double tau_minus = lambda_trt / (tau - 0.5) + 0.5;
    omega_P = 1.0 / tau;
    omega_M = 1.0 / tau_minus;

    // Initialize node and obstacles matrices
    node_matrix = Matrix<Node> ({NX,NY});
    obstacles = Matrix<bool> ({NX,NY});

    for (int i=0; i<NX; i++)
        for (int j=0; j<NY; j++)
        {   
            obstacles.set({i,j}, false);
        }
    
    std::ofstream param_file(outdir+"/parameters.txt", std::ios::trunc);    
    if (param_file.is_open()) {
        param_file << "# ---- LBM Simulation Parameters ----\n\n";
        param_file << "# Domain Size\n";
        param_file << "NX = " << NX << "\n";
        param_file << "NY = " << NY << "\n\n";

        param_file << "# Physical Parameters\n";
        param_file << "Re = " << Re << "\n";
        param_file << "boundary_velocity:\n";
        param_file << "  Vx = " << boundary_velocity.at(0) << "\n";
        param_file << "  Vy = " << boundary_velocity.at(1) << "\n\n";

        param_file << "# Derived Parameters\n";
        param_file << "nu (kinematic viscosity) = " << nu << "\n";
        param_file << "tau (relaxation time) = " << tau << "\n";
        param_file << "omega_P = " << omega_P << "\n";
        param_file << "lambda_trt = " << lambda_trt << "\n";
        param_file << "tau_minus = " << tau_minus << "\n";
        param_file << "omega_M = " << omega_M << "\n\n";

        param_file << "# Sigma for force calculation\n";
        param_file << "sigma = " << sigma << "\n\n";

        param_file << "# Output files:\n";
        param_file << "vel_data.txt (velocity field)\n";
        param_file << "lift_drag.txt (aerodynamic forces)\n\n";

       
    } 
    else {
        std::cerr << "Could not open/create 'parameters.txt'\n";
    }

    std::string user_input;
    size = 1; // Default in case no mask is added
    
    while (true) {
        std::cout << "\nDo you want to add an obstacle? (y/n): ";
        std::cin >> user_input;

        if (user_input == "n" || user_input == "N") {
            break;
        }
        param_file << "# Obstacles:\n";
        std::string mask_type;
        std::cout << "Which type of obstacle? (circle / rect / airfoil): ";
        std::cin >> mask_type;

        if (mask_type == "circle") {
            int radius, x_center, y_center;
            std::cout << "Enter radius: ";
            std::cin >> radius;
            std::cout << "Enter center x-coordinate: ";
            std::cin >> x_center;
            std::cout << "Enter center y-coordinate: ";
            std::cin >> y_center;
            size = create_circular_mask(radius, x_center, y_center, obstacles);
            object_count++;
     
            param_file << "  - type: circle\n";
            param_file << "    radius: " << radius << "\n";
            param_file << "    center_x: " << x_center << "\n";
            param_file << "    center_y: " << y_center << "\n";
            param_file << "\n";
        }
        else if (mask_type == "rect") {
            double rot;
            int width, height, x_pos, y_pos;
            std::cout << "Enter width: ";
            std::cin >> width;
            std::cout << "Enter height: ";
            std::cin >> height;
            std::cout << "Enter bottom-left x-coordinate: ";
            std::cin >> x_pos;
            std::cout << "Enter bottom-left y-coordinate: ";
            std::cin >> y_pos;
            std::cout << "Enter rotation [deg] (default 0.0°): ";
            std::cin >> rot;
            size = create_rectangular_mask(width, height, x_pos, y_pos, obstacles, rot);
            object_count++;

            param_file << "  - type: rect\n";
            param_file << "    width: " << width << "\n";
            param_file << "    height: " << height << "\n";
            param_file << "    x_pos: " << x_pos << "\n";
            param_file << "    y_pos: " << y_pos << "\n";
            param_file << "    rotation_deg: " << rot << "\n\n";
        }
        else if (mask_type == "airfoil") {
            int x_attack, y_attack, chord_length;
            double t;
            std::cout << "Enter x attack-coordinate of the airfoil: ";
            std::cin >> x_attack;
            std::cout << "Enter y attack-coordinate of the airfoil: ";
            std::cin >> y_attack;
            std::cout << "Enter chord length: ";
            std::cin >> chord_length;
            std::cout << "Enter profile Thickness (chord fraction, i.e. 0.12): ";
            std::cin >> t;
            size = create_airfoil_mask(chord_length, x_attack, y_attack, t, obstacles);
            object_count++;

            param_file << "  - type: airfoil\n";
            param_file << "    x_attack: " << x_attack << "\n";
            param_file << "    y_attack: " << y_attack << "\n";
            param_file << "    chord_length: " << chord_length << "\n";
            param_file << "    thickness: " << t << "\n\n";
        }
        else {
            std::cerr << "Unrecognized mask type. Please enter one of: circle, rect, airfoil." << std::endl;
        }
    }
    std::cout << "\nParameters saved in 'parameters.txt'" << std::endl;
    param_file.close();

    std::vector<int> boundary;
    bool isObstacle;
    #pragma omp parallel
    {
        #pragma omp for collapse(2)
        for (int i=0; i<NX; i++)
            for (int j=0; j<NY; j++)
            {
                #pragma omp critical
                {
                isObstacle = obstacles.getCopy(i,j);
                boundary = evaluateBoundary( {i,j} , obstacles);
                Node node(boundary, isObstacle, {i,j});
                node_matrix.set( {i,j}, node );
                node_matrix(i,j).initializeEquilibrium();
                }
                node_matrix(i,j).updateMacroscopic();
            }
    }
    // Create the output file for velocity
    std::ofstream file_velocity(outdir+"/vel_data.txt",std::ios::trunc);
    if (!file_velocity.is_open()) {
        std::cerr << "could not opene/create 'vel_data.txt'.\n";
        return;
    }
    file_velocity << NX << "\n" << NY << "\n";

    std::ofstream file_lift_drag(outdir+"/lift_drag.txt", std::ios::trunc);
        if (!file_lift_drag.is_open()) {
            std::cerr << "could not opene/create 'lift_drag.txt'.\n";
            return;
    }

}

void Lattice::simulate(int max_steps, int iter_per_frame)
{
    maxSteps = max_steps;
    ITERATIONS_PER_FRAME = iter_per_frame;

    std::ofstream param_file(outdir+"/parameters.txt", std::ios::app);    
    if (param_file.is_open()) {
    param_file << "# Time-Stepping\n";
        param_file << "maxSteps = " << maxSteps << "\n";
        param_file << "ITERATIONS_PER_FRAME = " << ITERATIONS_PER_FRAME << "\n";
        param_file << "ITERATIONS_PER_PROGRESS_UPDATE = " << ITERATIONS_PER_PROGRESS_UPDATE << "\n\n";
    }    
    param_file.close();

    const double t1 = 2.0*sigma*sigma;
    const double halfOmegaSum = (omega_P+omega_M)/2.0;
    const double halfOmegaSub = (omega_P-omega_M)/2.0;
    auto startTime = std::chrono::high_resolution_clock::now();

    for(int currentStep = 0; currentStep < maxSteps; currentStep++)
    {
        const double Vx = boundary_velocity.at(0)*(1.0 - std::exp(-static_cast<double>(currentStep*currentStep)/t1));
        const double Vy = boundary_velocity.at(1)*(1.0 - std::exp(-static_cast<double>(currentStep*currentStep)/t1));
        double Cd = 0.0, Cl = 0.0;

        #pragma omp parallel num_threads(64)
        {
            #pragma omp for collapse(2)
            for( int i=0; i<node_matrix.shape().at(0); i++)
                for (int j=0; j<node_matrix.shape().at(1); j++)
            {
                    node_matrix(i,j).applyInletBoundary(*this, {Vx,Vy});
                    node_matrix(i,j).applyZouHeBoundary(*this);
                    node_matrix(i,j).updateMacroscopic();
                    node_matrix(i,j).equilibriumCollision(omega_P, halfOmegaSum, halfOmegaSub);
                    node_matrix(i,j).bounceBack();
            }

            // Streaming
            #pragma omp for collapse(2)
            for( int i=0; i<node_matrix.shape().at(0); i++)
                for (int j=0; j<node_matrix.shape().at(1); j++)
                    node_matrix(i,j).streaming(*this);
        }

        if (object_count == 1)
        {
            for( int i=0; i<node_matrix.shape().at(0); i++)
                for (int j=0; j<node_matrix.shape().at(1); j++)
                    node_matrix(i,j).computeDragAndLift(Cd, Cl, 1.0, size , Vx);
        }
        std::ofstream file_velocity(outdir+"/vel_data.txt", std::ios::app);
        std::ofstream file_lift_drag(outdir+"/lift_drag.txt", std::ios::app);

        //Every ITERATIONS_PER_FRAME steps, save velocity data
        if (currentStep==1 || currentStep % ITERATIONS_PER_FRAME == 0) {
            for (int j = 0; j < node_matrix.shape().at(1); j++) {
                for (int i = 0; i < node_matrix.shape().at(0); i++) {
                    double vx = node_matrix(i,j).getVelocity().at(0); 
                    double vy = node_matrix(i,j).getVelocity().at(1); 
                    double v = sqrt(vx*vx + vy*vy);

                    if (node_matrix(i,j).isObstacle()==false && v>1){
                        std::cerr << "\nSimulation has diverged at #iter ... " << currentStep << std::endl;
                        std::cerr << "Problem happened at cell\n (i,j): (" << i << "," << j << ")" << std::endl;
                        return;
                    }

                    file_velocity << v << "\n";
                    if (object_count==1)
                        file_lift_drag << Cl << " " << Cd << "\n";
                }
            }
        }

        // Update the progress bar
        if (currentStep % ITERATIONS_PER_PROGRESS_UPDATE == 0 || currentStep == maxSteps - 1) {
            float progress = (static_cast<float>(currentStep) / maxSteps);
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
    std::cout << std::endl;
}

