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

    /**
     * @brief Applies boundary conditions.
     */
    void apply_boundary_conditions() override;

    /** 
     * @brief Generate a NACA 00tt profile 
     * @param x Coordinate x in which the output is evaluated
     * @param chord lenght of the profile
     * @return Coordinate y associated to x and the profile 
     */
    double naca_airfoil(double x, double chord);
    
public:
    std::vector<std::vector<bool>> is_solid; // 1->solid, 0->fluid

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
     * @param chord lenght of the profile
     * @param x x position from where airfoil profile will start
     */
    void create_airfoil_mask( double chord, double x);

    /**
     * @brief Generate a rectangular mask
     * @param b base of the rectangular
     * @param h height of the rectangular 
     * @param x x-coordinate of the bottom-left corner
     * @param y y-coordinate of the bottom-left corner 
     */
    void create_rectangular_mask( const unsigned int b, const unsigned int h, const unsigned int x, const unsigned int y);

    /**
     * @brief Generate a circular mask
     * @param r ray of the circle
     * @param x x-coordinate of center
     * @param y y-coordinate of center
     */
    void create_circular_mask( const unsigned int r, const unsigned int x, const unsigned int y);

    /**
     * @brief Evolves the system for a single iteration.
     */
    void evolution() override;

    /**
     * @brief Compute drag and lift on a solid body
     */
    void computeForces();

};

#endif