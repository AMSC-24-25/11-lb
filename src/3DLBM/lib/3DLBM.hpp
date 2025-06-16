#ifndef TD_LBM_HPP
#define TD_LBM_HPP

#include <vector>
#include <cmath>
#include <filesystem>
#ifdef USE_OPENMP
  #include <omp.h>
#endif

/**
 * @struct TDouble
 * @brief 3-component double vector for velocity, etc.
 */
struct TDouble {
    double x; ///< x-component
    double y; ///< y-component
    double z; ///< z-component
};

/**
 * @brief Construct a TDouble.
 * @param x x-component.
 * @param y y-component.
 * @param z z-component.
 * @return TDouble struct.
 */
inline TDouble Set_TDouble(double x, double y, double z) {
    return { x, y, z };
}

/**
 * @brief Dot product of two TDouble.
 * @param a first vector.
 * @param b second vector.
 * @return scalar product.
 */
inline double dot(const TDouble& a, const TDouble& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/**
 * @class TDLBM
 * @brief Implements the 3D D3Q19 Lattice–Boltzmann method with TRT collision model.
 */
class TDLBM {
protected:
    static constexpr int Q = 19; ///< Number of discrete directions.

    // Domain and simulation parameters
    int nx, ny, nz;      ///< grid dimensions
    double U_lid;        ///< lid velocity
    double Re;           ///< target Reynolds number

    // Numerical parameters
    double dx, dt;       ///< spatial and temporal steps
    double rho0;         ///< initial density
    double L;            ///< characteristic length
    double nu;           ///< kinematic viscosity
    double tau;          ///< relaxation time

    // Simulation data
    int ncells;                     ///< total cells
    std::vector<TDouble> e;         ///< direction vectors
    std::vector<double> f, f_new;   ///< distribution functions
    std::vector<double> rho;        ///< density field
    std::vector<TDouble> u, u_prev; ///< velocity fields

public:
    /**
     * @brief Constructor.
     * @param nx number of cells in x.
     * @param ny number of cells in y.
     * @param nz number of cells in z.
     * @param u_lid lid velocity.
     * @param Re target Reynolds number.
     */
    TDLBM(unsigned int nx, unsigned int ny, unsigned int nz, double u_lid, double Re);

    /**
     * @brief Equilibrium weight for a given direction.
     * @param i direction index.
     * @return equilibrium weight.
     */
    double weight(int i) const;

    /**
     * @brief Get the opposite direction index.
     * @param i direction index.
     * @return opposite direction index.
     */
    int get_opposite_direction(int i) const;

    /**
     * @brief Compute linear cell index.
     * @param x x-coordinate.
     * @param y y-coordinate.
     * @param z z-coordinate.
     * @return linear index.
     */
    int idx(int x, int y, int z) const;

    /**
     * @brief Compute linear distribution index.
     * @param x x-coordinate.
     * @param y y-coordinate.
     * @param z z-coordinate.
     * @param i direction index.
     * @return linear distribution index.
     */
    int idxi(int x, int y, int z, int i) const;

    /**
     * @brief Initialize density and velocity fields.
     * @param rho density vector.
     * @param u velocity vector.
     */
    void init_rho_v(std::vector<double>& rho, std::vector<TDouble>& u);

    /**
     * @brief Initialize distribution functions.
     * @param f distribution vector.
     * @param f_new temp distribution vector.
     * @param rho density vector.
     * @param u velocity vector.
     * @param e direction vectors.
     */
    void init_f(std::vector<double>& f,
                std::vector<double>& f_new,
                const std::vector<double>& rho,
                const std::vector<TDouble>& u,
                const std::vector<TDouble>& e);

    /**
     * @brief Set the D3Q19 direction vectors.
     * @param e vector of direction vectors.
     */
    void set_e_values(std::vector<TDouble>& e);

    /**
     * @brief TRT collision operator.
     * @param f distribution functions.
     * @param rho density field.
     * @param u velocity field.
     * @param e direction vectors.
     */
    void collide_trt(std::vector<double>& f,
                     std::vector<double>& rho,
                     std::vector<TDouble>& u,
                     const std::vector<TDouble>& e);

    /**
     * @brief Streaming step of the distributions.
     * @param f current distributions.
     * @param f_new post-streaming distributions.
     * @param e direction vectors.
     */
    void stream(const std::vector<double>& f,
                std::vector<double>& f_new,
                const std::vector<TDouble>& e);

    /**
     * @brief Apply boundary conditions (lid-driven cavity).
     * @param f distribution functions.
     * @param e direction vectors.
     */
    void apply_boundary_conditions(std::vector<double>& f,
                                   const std::vector<TDouble>& e);

    /**
     * @brief Perform one time step evolution.
     */
    void simulate();

    /**
     * @brief Equilibrium distribution function.
     * @param rho cell density.
     * @param w weight.
     * @param e direction vector.
     * @param u cell velocity.
     * @return equilibrium f.
     */
    double feq(double rho, double w, const TDouble& e, const TDouble& u) const;

    /**
     * @brief Compute density in a cell (read-only).
     * @param f distribution vector.
     * @param x x-coordinate.
     * @param y y-coordinate.
     * @param z z-coordinate.
     * @return cell density.
     */
    double compute_rho(const std::vector<double>& f, int x, int y, int z) const;

    /**
     * @brief Compute velocity in a cell (read-only).
     * @param f distribution vector.
     * @param rho cell density.
     * @param e direction vectors.
     * @param x x-coordinate.
     * @param y y-coordinate.
     * @param z z-coordinate.
     * @return cell velocity.
     */
    TDouble compute_u(const std::vector<double>& f,
                      double rho,
                      const std::vector<TDouble>& e,
                      int x, int y, int z) const;

    /**
     * @brief Get spatial step dx.
     * @return dx.
     */
    double get_dx() const { return dx; }

    /**
     * @brief Get relaxation time τ.
     * @return τ.
     */
    double get_tau() const { return tau; }

    /**
     * @brief Get kinematic viscosity ν.
     * @return ν.
     */
    double get_nu() const { return nu; }

    /**
     * @brief Get velocity field.
     * @return const reference to velocity vector.
     */
    const std::vector<TDouble>& get_u() const { return u; }

    /**
     * @brief Get density field.
     * @return const reference to density vector.
     */
    const std::vector<double>& get_rho() const { return rho; }
};

#endif // TD_LBM_HPP
