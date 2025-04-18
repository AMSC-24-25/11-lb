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
    const double t = 0.20; // profile thickness (expressed as chord%)
    return (5*t*chord) * (0.2969*sqrt(x/chord) - 0.1260*x/chord - 0.3516*x*x/chord/chord + 
            0.2843 *x*x*x/chord/chord/chord - 0.1016*x*x*x*x/chord/chord/chord/chord);
}

void WindTunnelLBM::create_airfoil_mask( double chord){

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            //check if position is well plqced
            if (i >= 0.5*(NX-chord) && i<0.5*(NX+chord) ) {
                // if so evaluate airfoil shape
                double y_airfoil = abs(naca_airfoil((i-0.5*(NX-chord)), chord));
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

        // open outlet (right boundary)
        #pragma omp for nowait
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < Q; k++) {
                velocity(NX-1, j, 0) = velocity(NX-2, j, 0);
                velocity(NX-1, j, 1) = velocity(NX-2, j, 1);
                field(NX-1, j, k) = field(NX-2, j, k); // zero-gradient
            }
        }

        // conditions on upper and lower boundaries
        #pragma omp for
        for (int i = 1; i < NX; i++) {
            velocity(i, NY - 1, 0) = 0.0;
            velocity(i, NY - 1, 1) = 0.0;
            velocity(i, 0, 0) = 0.0;
            velocity(i, 0, 1) = 0.0;
            for (int k = 0; k < Q; k++) {
                // bottom boundary 
                field(i, 0, k) = feq(k, i, 0) + field(i, 1, k) - feq(k, i, 1);
                // upper boundary (open outlet)
                field(i, NY - 1, k) = feq(k, i, NY - 1) + field(i, NY - 2, k) - feq(k, i, NY - 2);
            }
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
    apply_boundary_conditions(); // Apply boundary conditionsù

}
