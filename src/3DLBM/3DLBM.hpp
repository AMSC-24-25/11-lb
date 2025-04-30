#ifndef TD_LBM_HPP
#define TD_LBM_HPP
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <filesystem>
#include <omp.h>

/**
 * @class 3DLBM
 * @brief Implementation of the Lattice Boltzmann method (LBM) for the D3Q19 model.
 *
 * The LBM class represents a three-dimensional model with nineteen discrete velocities (D2Q9)
 * for fluid dynamics simulations based on the lattice Boltzmann approach.
*/


struct TDouble { double x, y, z; };

/** 
    * @brief Compute TDouble.
    * @param x direction x.
    * @param y direction y.
    * @param z direction z.
    * @return TDouble
    */
   inline TDouble Set_TDouble(double x, double y, double z){
    return {x, y, z};
}

class TDLBM {
protected:  
    static constexpr int Q = 19; ///< Number of discrete velocities (D3Q19).

    //Main variables
     
    int nx;
    int ny;
    int nz;
    int ncells;
    double dx;
    double dt;
    double rho0;
    double U_lid;
    double L;
    double nu;
    double tau;
    int Re;
    int Steps;
    int ITERATIONS_PER_FRAME;
    int ITERATIONS_PER_PROGRESS_UPDATE;
    std::vector<TDouble> e;
    std::vector<double> f;
    std::vector<double> f_new;
    std::vector<double> rho;
    std::vector<TDouble> u;
    std::vector<TDouble> u_prev;
    

    

    /** 
    * @brief Compute a dot.
    * @param a  TDouble.
    * @param b  Tdouble.
    */
    inline double dot(const TDouble& a, const TDouble& b){
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    /** 
    * @brief Compute viscosity nu.
    * @param U velocity.
    * @param L length .
    * @return nu
    */
    inline double compute_nu(double U, double L, int Re){
        return (U * L) / Re;
    }

    /** 
    * @brief Compute viscosity nu.
    * @param nu viscosity.
    * @return tau
    */
    inline double compute_tau(double nu){
        return (6.0 * nu + 1.0) / 2.0;
    }

    
   /** 
    * @brief Equilibrium weights
    * @param i index.
    */
    double weight(int i){
        if (i == 0)
            return 1.0/3.0;
        else if (i >= 1 && i < 7)
            return 1.0/18.0;
        else
            return 1.0/36.0;
    }

    /** 
    * @brief Set the 19 direction vector of D3Q19
    * @param direction array.
    */
     void set_e_values(std::vector<TDouble>& e);

    /** 
    * @brief function to obtain the opposite direction.
    * @param i index direction.
    */
    int get_opposite_direction(int i){
        switch(i) {
            case 0: return 0;  // riposo
            case 1: return 2;  // +x -> -x
            case 2: return 1;  // -x -> +x
            case 3: return 4;  // +y -> -y
            case 4: return 3;  // -y -> +y
            case 5: return 6;  // +z -> -z
            case 6: return 5;  // -z -> +z
            case 7: return 8;  // +x,+y -> -x,-y
            case 8: return 7;  // -x,-y -> +x,+y
            case 9: return 10; // -x,+y -> +x,-y
            case 10: return 9; // +x,-y -> -x,+y
            case 11: return 12; // +x,+z -> -x,-z
            case 12: return 11; // -x,-z -> +x,+z
            case 13: return 14; // +y,+z -> -y,-z
            case 14: return 13; // -y,-z -> +y,+z
            case 15: return 16; // -x,+z -> +x,-z
            case 16: return 15; // +x,-z -> -x,+z
            case 17: return 18; // -y,+z -> +y,-z
            case 18: return 17; // +y,-z -> -y,+z
            default: return 0;
        }
    }

    /** 
    * @brief compute the equilibrium function.
    * @param rho.
    * @param w weight.
    * @param e direction array.
    * @param u .
    * @return equilibrium function
    */
    double feq(double rho, double w, const TDouble& e, const TDouble& u);

    /** 
    * @brief init the density and velocity in every cells.
    * @param rho.
    * @param u.
    */
    void init_rho_v(std::vector<double>& rho, std::vector<TDouble>& u);

    /** 
    * @brief init function f.
    * @param f.
    * @param f_new new function. 
    * @param rho. 
    * @param u.
    * @param e array directions.
    */
    void init_f(std::vector<double>& f, std::vector<double>& f_new,
    const std::vector<double>& rho, const std::vector<TDouble>& u,
    const std::vector<TDouble>& e);

    /** 
    * @brief calculate rho value.
    * @param f function.
    * @param x direction x. 
    * @param y direction y. 
    * @param z direction z.
    */
    double compute_rho(const std::vector<double>& f, int x, int y, int z);

    /** 
    * @brief calculate u value for each direction.
    * @param f function.
    * @param rho. 
    * @param e direction array;
    * @param x direction x. 
    * @param y direction y. 
    * @param z direction z.
    */
    TDouble compute_u(const std::vector<double>& f, double rho,
    const std::vector<TDouble>& e, int x, int y, int z);

    /** 
    * @brief collision operator TRT.
    * @param f function.
    * @param rho. 
    * @param u. 
    * @param e.
    */
    void collide_trt(std::vector<double>& f, std::vector<double>& rho,
    std::vector<TDouble>& u, const std::vector<TDouble>& e);

    /** 
    * @brief streaming operatore with bounce-back for node out of the domain.
    * @param f function.
    * @param f_new new function.  
    * @param e.
    */
    void stream(const std::vector<double>& f, std::vector<double>& f_new, const std::vector<TDouble>& e);

    /** 
    * @brief apply the boundary condition for the driven cavity: Zou-he and no slip condition.
    * @param f.
    * @param e directions array.
    */
    void apply_boundary_conditions(std::vector<double>& f, const std::vector<TDouble>& e);



public:
    
    
    
    /**
     * @brief Constructor of the TDLBM class.
     * @param nx Domain size in the x-direction.
     * @param ny Domain size in the y-direction.
     * @param nz Domain size in the y-direction.
     * @param u_lid Velocity at the top boundary.
     * @param Re Reynolds number.
     */
    TDLBM(unsigned int nx, unsigned int ny, unsigned int nz, double u_lid, double Re);

    /** 
    * @brief helper to find index.
    * @param x direction x.
    * @param y direction y.
    * @param z direction z.
    * @return index
    */
    inline int idx(int x, int y, int z)const{
    return x + nx * y + nx * ny * z;
}

/** 
* @brief helper to find distribution vector index in a specific direction.
* @param x direction x.
* @param y direction y.
* @param z direction z.
* @param i index.
* @return index for distribution vector
*/
  inline int idxi(int x, int y, int z, int i){
    return x + (((y + z * ny) * nx) + (nx * ny * nz * i));
}


    /**
     * @brief Evolves the system for a specified number of iterations.
     * @param iterations Number of iterations.
     */
    void evolution();

    /**
     * @brief getter method.
     */
    const std::vector<double>& get_rho() const { return rho; }
    /**
     * @brief getter method.
     */
    const std::vector<TDouble>& get_u() const { return u; }
    /**
     * @brief getter method.
     */
    double get_dx() const { return dx; }
    /**
     * @brief getter method.
     */
    double get_tau() const { return tau; }

    /**
     * @brief getter method.
     */
    double get_nu() const { return nu; }
    
};
#endif