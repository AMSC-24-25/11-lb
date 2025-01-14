#ifndef WINDTUNNEL_LBM_HPP
#define WINDTUNNEL_LBM_HPP

#include "LBM.hpp"


/**
 * @class WindTunnelLBM
 * @brief Specialized class for simulating a 2D wind tunnel using the Lattice Boltzmann method.
 */
class WindTunnelLBM: public LBM {
private:
    double inlet_velocity_x; ///< Velocity component in the x-direction.
    double inlet_velocity_y; ///< Velocity component in the y-direction.
    std::vector<std::vector<bool>> is_solid; // 1->solid, 0->fluid

    /**
     * @brief Applies boundary conditions.
     */
    void apply_boundary_conditions() override;

    /** 
     * @brief Generate a NACA 00tt profile 
     * @param x Coordinate x in which the output is evaluated
     * @return Coordinate y associated to x and the profile 
     */
    double naca_airfoil(double x);
    
public:
    double chord;

    /**
     * @brief Constructor of the LBM class.
     * @param nx Domain size in the x-direction.
     * @param ny Domain size in the y-direction.
     * @param inlet_velocity_x x component velocity at the top boundary.
     * @param inlet_velocity_y y component velocity at the top boundary.
     * @param Re Reynolds number.
     */
    WindTunnelLBM(unsigned int nx, unsigned int ny, double inlet_velocity_x, double inlet_velocity_y, double Re);

    /**
     * @brief Generate NACA airfoil mask
     */
    void create_airfoil_mask();

    /**
     * @brief Evolves the system for a single iteration.
     */
    void evolution() override;
};






#endif