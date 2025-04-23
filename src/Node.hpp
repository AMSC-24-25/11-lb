#ifndef NODE_HPP
#define NODE_HPP

#include <vector>

// Forward declaration of Lattice class
class Lattice;

/**
 * @brief Represents a single node (cell) in a 2D Lattice Boltzmann grid using the D2Q9 model.
 */
class Node {
public:
    /**
     * @brief Constructs a Node with given boundary conditions, obstacle status, and position.
     * 
     * @param boundary Boundary condition flags for each direction.
     * @param isObstacle Whether the node represents a solid obstacle.
     * @param position Position of the node in the lattice grid.
     */
    Node(const std::vector<int>& boundary, const bool isObstacle, const std::vector<int>& position);

    /// Default constructor
    Node() = default;

    // --- Core simulation methods ---

    /**
     * @brief Updates macroscopic properties (density and velocity) from the distribution functions.
     */
    void updateMacroscopic();

    /**
     * @brief Performs the equilibrium-based collision using TRT model.
     * 
     * @param omegaP Primary relaxation parameter.
     * @param halfOmegaSum 0.5 * (omegaP + omegaM)
     * @param halfOmegaSub 0.5 * (omegaP - omegaM)
     */
    void equilibriumCollision(double omegaP, double halfOmegaSum, double halfOmegaSub);

    /**
     * @brief Initializes the distribution functions to equilibrium based on the problem type.
     * 
     */
    void initializeEquilibrium();

    /**
     * @brief Streams the distribution functions from neighboring nodes.
     * 
     * @param lattice Reference to the global lattice.
     */
    void streaming(Lattice& lattice);

    /**
     * @brief Applies inlet boundary condition (e.g., moving lid).
     * 
     * @param lattice Reference to the global lattice.
     * @param boundary_velocity Current lid velocity.
     */
    void applyInletBoundary(Lattice& lattice, std::vector<double> boundary_velocity);

    /**
     * @brief Applies Zou/He boundary condition for walls.
     * 
     * @param lattice Reference to the global lattice.
     */
    void applyZouHeBoundary(Lattice& lattice);

    /**
     * @brief Performs bounce-back for solid obstacle cells.
     */
    void bounceBack();

    /**
     * @brief Accumulates drag and lift contributions for this node.
     * 
     * @param Cx Drag Coefficient.
     * @param Cy Lift coefficient.
     */
    void computeDragAndLift(double &Cx, double &Cy, double dRef, double LRef, double URef);

    // --- Getters ---
    const double& getDensity() const;
    const std::vector<double>& getVelocity() const;
    const std::vector<double>& getF(bool age) const;
    const std::vector<int>& getBoundary() const;
    bool isObstacle() const;

    // --- Setters ---
    void setFAtIndex(int index, const double& value);

private:
    std::vector<double> f;         ///< Distribution function (9 directions for D2Q9)
    std::vector<double> newF;      ///< Post-streaming values

    std::vector<double> velocity;  ///< Macroscopic velocity (ux, uy)
    double density = 1.0;          ///< Macroscopic density

    std::vector<int> boundary;     ///< Boundary flags for each relevant direction
    bool obstacle = false;         ///< True if this node is a solid obstacle
    std::vector<int> position;     ///< Position in the lattice grid

    std::vector<double> w = {4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
                                    1.0 / 36.0};
                                    ///< Weights
    std::vector<std::vector<double>> dir = {{0, 0}, {1, 0}, {0, -1}, {-1, 0}, {0, 1}, {1, -1}, {-1, -1}, {-1, 1}, {1, 1}};
    std::vector<std::vector<double>> dir_by_dim = {{0, 1, 0, -1, 0, 1, -1, -1, 1}, {0, 0, -1, 0, 1, -1, -1, 1, 1}};
    std::vector<unsigned int> opp = {0, 3, 4, 1, 2, 7, 8, 5, 6};

};

#endif // NODE_HPP
