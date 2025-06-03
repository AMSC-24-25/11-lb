#include "Lattice.hpp"
#include "Masks.cpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip> 
#include <omp.h>

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

Lattice::Lattice(unsigned NX_, unsigned steps_, double Re_,
                 bool useMask, const std::string &maskType,
                 double maskSize, const std::string &outDir)
  : NX(NX_), NY(NX_), steps(steps_), Re(Re_),
    outDir_(outDir), useMask_(useMask),
    maskType_(maskType), maskSize_(maskSize),
{
    const unsigned int NY = NX_;
    ITERATIONS_PER_FRAME = 25;
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
    #pragma omp parallel
    {
    #pragma omp parallel for collapse(2)
    for (int i=0; i<NX; i++)
        for (int j=0; j<NY; j++)
        {   
            obstacles.set({i,j}, false);
        }
    
    // CREATE MASK
    // Circular
    //size = create_circular_mask(30, 50, 150, obstacles);

    // Step
    /*for (int i=0;i<35;i++){
        size = create_rectangular_mask (1,9+i,1+i,0,obstacles);
        create_rectangular_mask(1,9+i,1+i, 99-9-i, obstacles);
    }
    for (int i=0;i<37;i++){
        size = create_rectangular_mask (1,9+i,70+i,0,obstacles);
        create_rectangular_mask(1,9+i,70+i, 99-9-i, obstacles);
    }*/
    // NACA 00xx Airfoil
    create_airfoil_mask( 100, 100, obstacles );

    std::vector<int> boundary;
    bool isObstacle;
    #pragma omp parallel for collapse(2)
    for (int i=0; i<NX; i++)
        for (int j=0; j<NY; j++)
        {
            isObstacle = obstacles.getCopy(i,j);
            boundary = evaluateBoundary( {i,j} , obstacles);
            Node node(boundary, isObstacle, {i,j});
            node_matrix.set( {i,j}, node );
            node_matrix(i,j).initializeEquilibrium();
            node_matrix(i,j).updateMacroscopic();
        }
    }
    // Create the output file for velocity
    std::ofstream file_velocity("vel_data.txt");
    if (!file_velocity.is_open()) {
        std::cerr << "could not opene/create 'vel_data.txt'.\n";
        return;
    }
    file_velocity << NX << "\n" << NY << "\n";
    std::ofstream file_lift_drag("lift_drag.txt");
    if (!file_lift_drag.is_open()) {
        std::cerr << "could not opene/create 'lift_drag.txt'.\n";
        return;
    }
}

void Lattice::simulate()
{
    const double t1 = 2.0*sigma*sigma;
    const double halfOmegaSum = (omega_P+omega_M)/2.0;
    const double halfOmegaSub = (omega_P-omega_M)/2.0;
    auto startTime = std::chrono::high_resolution_clock::now();

    while (currentStep<=maxSteps)
    {
        const double Vx = boundary_velocity.at(0)*(1.0 - std::exp(-static_cast<double>(currentStep*currentStep)/t1));
        const double Vy = boundary_velocity.at(1)*(1.0 - std::exp(-static_cast<double>(currentStep*currentStep)/t1));
        double Cd = 0.0, Cl = 0.0;
        #pragma omp parallel
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
        #pragma omp parallel for collapse(2)
        for( int i=0; i<node_matrix.shape().at(0); i++)
            for (int j=0; j<node_matrix.shape().at(1); j++)
                node_matrix(i,j).streaming(*this);

        #pragma omp parallel for collapse(2)
        for( int i=0; i<node_matrix.shape().at(0); i++)
            for (int j=0; j<node_matrix.shape().at(1); j++)
                node_matrix(i,j).computeDragAndLift(Cd, Cl, 1.0, size , Vx);
        }

        std::ofstream file_velocity("vel_data.txt", std::ios::app);
        std::ofstream file_lift_drag("lift_drag.txt", std::ios::app);

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
                std::cout << "\rProgress: " << std::fixed << std::setprecision(2) << progress << "% completed "
                          << "| Elapsed Time: " << elapsedTime << "s, "
                          << "Remaining Time (estimated): " << remainingTime << "s"
                          << std::flush;
            } else {
                std::cout << "\rProgress: " << std::fixed << std::setprecision(2) << progress * 100 << "% completed "
                          << "| Elapsed Time: " << elapsedTime << "s, "
                          << "Remaining Time (estimated): calculating..."
                          << std::flush;
            }
        }
        currentStep += 1;
    }
    // std::cout << "\a" << std::flush; 
}