#include "Auxiliary.cpp"
#include "Node.hpp"
#include <cmath>
#include <iostream>

/* FILE WITH FUNCTION TO CREATE VARIOUS TYPES OF MASKS   */
/* EVERY MASK WILL GENERATE AN OBJECT IN THE SIMULATION  */

double naca_airfoil( double x, double chord )
{
    const double t = 0.15; // profile thickness (expressed as chord%)
    return (5*t*chord) * (0.2969*sqrt(x/chord) - 0.1260*x/chord - 0.3516*x*x/chord/chord + 
            0.2843 *x*x*x/chord/chord/chord - 0.1016*x*x*x*x/chord/chord/chord/chord);
}

double create_airfoil_mask( double chord, double x, Matrix<bool>& obstacles){

    const unsigned int NX = obstacles.shape().at(0);
    const unsigned int NY = obstacles.shape().at(1);

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            //check if position is well placed
            if (i >= x && i<(x+chord) ) {
                // if so evaluate airfoil shape
                double y_airfoil = abs(naca_airfoil((i-x), chord));
                // check if (i, j) is outside or belongs to the airfoil
                if (j <= NY/2+y_airfoil && j>=NY/2-y_airfoil) {
                    // if belongs
                    obstacles.set({i,j}, true);
                }
            }
        }
    }
    return chord;
}

#include <cmath> // for cos, sin, M_PI

double create_rectangular_mask(
    const unsigned int b, const unsigned int h, 
    const int x, const int y, 
    Matrix<bool>& obstacles,
    double angle_degrees = 0.0) // default no rotation
{
    double theta = angle_degrees * M_PI / 180.0; // degrees to radians
    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);

    // Center of the rectangle
    double cx = b / 2.0;
    double cy = h / 2.0;

    for (int i = 0; i <= b; i++) {
        for (int j = 0; j <= h; j++) {
            // Rotate around center
            double rel_x = i - cx;
            double rel_y = j - cy;

            double rot_x = cos_theta * rel_x - sin_theta * rel_y + cx;
            double rot_y = sin_theta * rel_x + cos_theta * rel_y + cy;

            // Map back to integer grid
            int final_x = std::round(rot_x) + x;
            int final_y = std::round(rot_y) + y;

            // Optional: Check if final_x and final_y are inside obstacles bounds
            obstacles.set({final_x, final_y}, true);
        }
    }
    return std::max(b,h);
}

double create_circular_mask(
    const unsigned int radius,
    const int center_x,
    const int center_y,
    Matrix<bool>& obstacles)
{
    const int r_squared = radius * radius;

    for (int i = -static_cast<int>(radius); i <= static_cast<int>(radius); i++) {
        for (int j = -static_cast<int>(radius); j <= static_cast<int>(radius); j++) {
            // Check if (i,j) is inside the circle
            if (i*i + j*j <= r_squared) {
                int final_x = center_x + i;
                int final_y = center_y + j;

                // check matrix bounds here!
                obstacles.set({final_x, final_y}, true);
            }
        }
    }
    return radius;
}

