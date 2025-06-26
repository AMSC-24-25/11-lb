#pragma once
#include "Auxiliary.hpp"   // serve Matrix<bool>

double naca_airfoil(double x, double chord, double t);

double create_airfoil_mask(double chord, double x_attack, double y_attack,
                           double t, Matrix<bool>& obstacles);

double create_rectangular_mask(unsigned int b, unsigned int h,
                               int x, int y, Matrix<bool>& obstacles,
                               double angle_degrees = 0.0);

double create_circular_mask(unsigned int radius, int cx, int cy,
                            Matrix<bool>& obstacles);
