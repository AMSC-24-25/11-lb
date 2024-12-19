#ifndef LBM_HPP
#define LBM_HPP


/**
 * @class LBM
 * @brief Implementation of the Lattice Boltzmann method (LBM) for the D2Q9 model.
 *
 * The LBM class represents a two-dimensional model with nine discrete velocities (D2Q9)
 * for fluid dynamics simulations based on the lattice Boltzmann approach.
 */
class LBM {
private:
    static constexpr int D = 2; ///< Number of dimensions (2D).
    static constexpr int Q = 9; ///< Number of discrete velocities (D2Q9).
    
    /**
     * @brief Discrete velocities of the D2Q9 model.
     *
     * The array contains the direction of velocities for each node.
     * Size: Q * D.
     * The order for the elements in this (and the next) array is not random.
     * The elements are ordered so that, when executing an iteration, the amound of
     * cache-miss is reduced. Ordering these arrays lead to a observable speedup.
     */
    static constexpr int e[Q * D] = {1,1,  1,0,  1,-1,  0,1,  0,0,  0,-1,  -1,1,  -1,0,  -1,-1};

    /**
     * @brief Equilibrium weights for each direction.
     *
     * Used in the equilibrium function.
     */
    static constexpr double w[Q] = {
        1.0/36.0,
        1.0/9.0,
        1.0/36.0,
        1.0/9.0,
        4.0/9.0,
        1.0/9.0,
        1.0/36.0,
        1.0/9.0,
        1.0/36.0
    };

    const double u_lid; ///< Velocity imposed on the top boundary.

    // Main variables
    double *rho;   ///< Current macroscopic density.
    double *rho2;  ///< Updated density.
    double *u;     ///< Current macroscopic velocity.
    double *u2;    ///< Updated velocity.
    double *f;     ///< Current distribution function.
    double *F;     ///< Updated distribution function.

    // Simulation parameters
    double c;     ///< Speed of sound (lattice units).
    double Re;    ///< Reynolds number.
    double dx, dy; ///< Spatial steps.
    double Lx, Ly; ///< Domain dimensions.
    double dt;    ///< Time step.
    double rho0;  ///< Initial density.
    double tau_f; ///< Relaxation time.
    double nu;    ///< Kinematic viscosity.

    /**
     * @brief Helper function to access directions.
     * @param d Discrete direction.
     * @param v Dimension (0 = x, 1 = y).
     * @return Velocity value.
     */
    static inline const int& direction(unsigned int d, unsigned int v){
        return e[d*D + v];
    }

    /** 
    * @brief Returns the density at a specific location.
    * @param x Position x.
    * @param y Position y.
    * @return Reference to the updated density.
    */
    inline double& density(unsigned int x, unsigned int y){
        return rho[(NY)*x+y];
    }

    /**
    * @brief Returns the density at a specific location (updated).
    * @param x Position x.
    * @param y Position y.
    * @return Reference to the updated density.
    */
    inline double& density_2(unsigned int x, unsigned int y){
        return rho2[(NY)*x+y];
    }

    /**
    * @brief Returns the velocity in a specific direction and position.
    * @param x Position x.
    * @param y Position y.
    * @param d Dimension (0 = x, 1 = y).
    * @return Reference to the velocity.
    */
    inline double& velocity(unsigned int x, unsigned int y, unsigned int d){
        return u[((NY)*x+y)*D + d];
    }

    /**
    * @brief Returns the updated velocity in a specific direction and position.
    * @param x Position x.
    * @param y Position y.
    * @param d Dimension (0 = x, 1 = y).
    * @return Reference to the updated velocity.
    */
    inline double& velocity_2(unsigned int x, unsigned int y, unsigned int d){
        return u2[((NY)*x+y)*D + d];
    }

    /**
    * @brief Returns the field value at a specific location and direction.
    * @param x Position x.
    * @param y Position y.
    * @param d Discrete direction.
    * @return Reference to the field value.
    */
    inline double& field(unsigned int x, unsigned int y,unsigned int d){
        return f[((NY)*x+y)*Q + d];
    }

    /**
    * @brief Returns the updated field value at a specific location and direction.
    * @param x Position x.
    * @param y Position y.
    * @param d Discrete direction.
    * @return Reference to the updated field value.
    */
    inline double& field_2(unsigned int x, unsigned int y,unsigned int d){
        return F[((NY)*x+y)*Q + d];
    }

    /**
     * @brief Computes the equilibrium function.
     * @param k Discrete direction.
     * @param x x-position.
     * @param y y-position.
     * @return Value of the equilibrium function.
     */
    double feq(unsigned int k, unsigned int x, unsigned int y);

    /**
     * @brief Performs collisions and updates velocities and densities.
     */
    void compute();  

    /**
     * @brief Applies boundary conditions.
     */
    void apply_boundary_conditions();

public : 
    const int NX; ///< Domain size in the x-direction.
    const int NY; ///< Domain size in the y-direction.


    /**
     * @brief Constructor of the LBM class.
     * @param nx Domain size in the x-direction.
     * @param ny Domain size in the y-direction.
     * @param u_lid Velocity at the top boundary.
     * @param Re Reynolds number.
     */
    LBM(unsigned int nx, unsigned int ny, double u_lid, double Re);
    ~LBM();

    /**
     * @brief Evolves the system for a single iteration.
     */
    void evolution();

    /**
     * @brief Evolves the system for a specified number of iterations.
     * @param iterations Number of iterations.
     */
    void evolution(unsigned int iterations);

    /**
     * @brief Returns the velocity at a specific location.
     * @param x x-position.
     * @param y y-position.
     * @param d Dimension (0 = x, 1 = y).
     * @return Velocity value.
     */
    inline const double& get_vel(unsigned int x, unsigned int y, unsigned int d) const {
        return u[((NY)*x+y)*D + d];
    }

}; 
#endif