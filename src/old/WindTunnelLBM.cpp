#include "WindTunnelLBM.hpp"
#include <iostream>

WindTunnelLBM::WindTunnelLBM(unsigned int nx_, unsigned int ny_, double inlet_velocity_x_,
                                double inlet_velocity_y_, double Re_)
    : LBM(nx_, ny_, 0.0, Re_), inlet_velocity_x(inlet_velocity_x_), inlet_velocity_y(inlet_velocity_y_),
                                is_solid(nx_, std::vector<bool>(ny_, false)) {

        // Redefinitions of nu and tau_f because they depend on u_lid in LBM.cpp
        nu = inlet_velocity_x*Lx/Re; // Kinematic viscosity
        tau_f = nu*3 + 0.5; // Relaxation time

        std::cout << tau_f << std::endl;

        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                is_solid[i][j] = false;
            }
        }
}


void WindTunnelLBM::create_rectangular_mask( const unsigned int b, const unsigned int h, 
                                                const unsigned int x, const unsigned int y)
{
    for (unsigned int i=0; i<b+1; i++)
        for (unsigned int j=0; j<b+1; j++)
            is_solid[x+i][y+j] = true;
}


double WindTunnelLBM::naca_airfoil( double x, double chord )
{
    const double t = 0.15; // profile thickness (expressed as chord%)
    return (5*t*chord) * (0.2969*sqrt(x/chord) - 0.1260*x/chord - 0.3516*x*x/chord/chord + 
            0.2843 *x*x*x/chord/chord/chord - 0.1016*x*x*x*x/chord/chord/chord/chord);
}

void WindTunnelLBM::create_airfoil_mask( double chord, double x){

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            //check if position is well placed
            if (i >= x && i<(x+chord) ) {
                // if so evaluate airfoil shape
                double y_airfoil = abs(naca_airfoil((i-x), chord));
                // check if (i, j) is outside or belongs to the airfoil
                if (j <= NY/2+y_airfoil && j>=NY/2-y_airfoil) {
                    // if belongs
                    is_solid[i][j] = true;
                }
            }
        }
    }
}

void WindTunnelLBM::apply_boundary_conditions() {
    #pragma omp parallel
    {
        // inlet conditions (left boundary)
        #pragma omp for nowait
        for (int j = 0; j < NY; j++) {
            double u_in_x = inlet_velocity_x;
            double u_in_y = inlet_velocity_y; 
            velocity(0, j, 0) = u_in_x;
            velocity(0, j, 1) = u_in_y;
        
            density(0,j) = 1.0;

        for (int k = 0; k < Q; k++) {
                double eu = direction(k, 0) * u_in_x + direction(k, 1) * u_in_y;
                field(0, j, k) = w[k] * rho0 * (1 + 3 * eu + 4.5 * eu * eu - 1.5 * (u_in_x * u_in_x + u_in_y * u_in_y));
            }
        }

        // right boundary open outlet
        #pragma omp for nowait
        for (int j = 0; j < NY; j++) {
                field(NX - 1, j, 6) = field(NX - 2, j, 6);
                field(NX - 1, j, 7) = field(NX - 2, j, 7);
                field(NX - 1, j, 8) = field(NX - 2, j, 8);
        }

        // conditions on upper and lower boundaries
        for (int i = 1; i < NX - 1; i++) {
                field(i, 0, 2) = field(i, 1, 2);       // upp
                field(i, 0, 3) = field(i, 1, 3);       // upp
                field(i, 0, 4) = field(i, 1, 4);       // upp

                field(i, NY-1, 2) = field(i, NY-2, 2); // low
                field(i, NY-1, 3) = field(i, NY-2, 3); // low
                field(i, NY-1, 4) = field(i, NY-2, 4); // low
        }

        // solid boundaries (noslip linear distribution)
        #pragma omp for
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                if (is_solid[i][j]) {
                    // zero velocity on solid body
                    velocity(i, j, 0) = 0.0;
                    velocity(i, j, 1) = 0.0;
                    for (int k = 0; k < Q; k++) {
                        field(i, j, k) = field(i, j, Q-k-1);
                        }
                }
            }
        }
    }
}

void WindTunnelLBM::evolution() {
    compute();      // Collision and Update macroscopic density and velocities
    computeForces();
    apply_boundary_conditions(); // Apply boundary conditionsù
}

void WindTunnelLBM::computeForces() {
    double drag = 0.0;
    double lift = 0.0;
    
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            if (is_solid[i][j]) { // If node is a solid body
                for (int k = 0; k < Q; k++) {
                    // Drag (X)
                    drag += field(i, j, k) * direction(k, 0);

                    // Lift (Y)
                    lift += field(i, j, k) * direction(k, 1);
                }
            }
        }
    }
    //std::cout << "Drag: " << drag << " Lift: " << lift << std::endl;
}

