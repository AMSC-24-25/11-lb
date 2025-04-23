#include "Lattice.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip> 

Lattice::Lattice()
{
    /*  TO SET MANUALLY  */
    // Hard coding variables
    const unsigned int NX = 100;
    const unsigned int NY = 100;
    const double Re = 1000;
    maxSteps = 3000;
    ITERATIONS_PER_FRAME = 20;
    ITERATIONS_PER_PROGRESS_UPDATE = 10;
    boundary_velocity.at(0) = 0.1; // Vx
    boundary_velocity.at(1) = 0.1; // Vy

    // Calculating other parameters
    const double nu = 2/3*NX*boundary_velocity.at(0)/Re;
    const double tau = 0.5 + nu * 3.0;
    sigma = 10.0*NX;
    const double lambda_trt = 1.0 / 4.0;
    const double tau_minus = lambda_trt / (tau - 0.5) + 0.5;
    omega_P = 1.0 / tau;
    omega_M = 1.0 / tau_minus;

    // Initialize node and obstacles matrices
    node_matrix = Matrix<Node> ({NX,NY});
    obstacles = Matrix<bool> ({NX,NY});
    for (int i=0; i<NX; i++)
        for (int j=0; j<NY; j++)
        {   obstacles(i,j) = false;
            if (i<=40 && i>=20 && j<=60 && j>=40)
            {
                obstacles(i,j) = true;
            }
        }
    std::vector<int> boundary;
    bool isObstacle;
    for (int i=0; i<NX; i++)
        for (int j=0; j<NY; j++)
        {
            isObstacle = obstacles.getCopy(i,j);
            boundary = evaluateBoundary( {i,j} , obstacles);
            node_matrix.set( {i,j}, Node(boundary, isObstacle, {i,j}) );
            node_matrix(i,j).initializeEquilibrium();
            node_matrix(i,j).initializeEquilibrium();
        }    

    // Create the output file for velocity
    std::ofstream file_velocity("vel_data.txt");
    if (!file_velocity.is_open()) {
        std::cerr << "could not opene/create 'vel_data.txt'.\n";
        return;
    }
    file_velocity << NX << "\n" << NY << "\n";
}

void Lattice::simulate()
{
    const double t1 = 2.0*sigma*sigma;
    const double halfOmegaSum = (omega_P+omega_M)/2;
    const double halfOmegaSub = (omega_P-omega_M)/2;
    double Cd = 0.0, Cl = 0.0;
    double size = 20; //has to be set to the side of the obstacle 

    while (currentStep<=maxSteps)
    {
        const double Vx = boundary_velocity.at(0)*(1.0 - std::exp(-static_cast<double>(currentStep*currentStep)/t1));
        const double Vy = boundary_velocity.at(1)*(1.0 - std::exp(-static_cast<double>(currentStep*currentStep)/t1));

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
        for( int i=0; i<node_matrix.shape().at(0); i++)
            for (int j=0; j<node_matrix.shape().at(1); j++)
                node_matrix(i,j).streaming(*this);

        for( int i=0; i<node_matrix.shape().at(0); i++)
            for (int j=0; j<node_matrix.shape().at(1); j++)
                node_matrix(i,j).computeDragAndLift(Cd, Cl, 1.0, size , boundary_velocity.at(0));

        auto startTime = std::chrono::high_resolution_clock::now();

        //Every ITERATIONS_PER_FRAME steps, save velocity data
        if (currentStep=1 || currentStep % ITERATIONS_PER_FRAME == 0) {
            for (int j = 0; j < node_matrix.shape().at(1); j++) {
                for (int i = 0; i < node_matrix.shape().at(0); i++) {
                    double vx = node_matrix(i,j).getVelocity().at(0); 
                    double vy = node_matrix(i,j).getVelocity().at(1); 
                    double v = sqrt(vx*vx + vy*vy);

                    std::ofstream file_velocity("vel_data.txt");
                    if (!file_velocity.is_open()) {
                        std::cerr << "could not opene/create 'vel_data.txt'.\n";
                        return;
                    } 
                    file_velocity << v << "\n";

                    std::ofstream file_lift_drag("lift_drag.txt");
                    if (!file_lift_drag.is_open()) {
                        std::cerr << "could not opene/create 'lift_drag.txt'.\n";
                        return;
                    } 
                    file_lift_drag << Cl << " " << Cd << "\n";
                }
            }
        }

        // Update the progress bar
        if (currentStep % ITERATIONS_PER_PROGRESS_UPDATE == 0 || currentStep == maxSteps - 1) {
            float progress = (static_cast<float>(currentStep) / maxSteps);
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
        currentStep++;
    }
}