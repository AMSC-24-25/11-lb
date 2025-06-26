#pragma once

#include <string>
#include "3DLBM.hpp"

/**
 * @class VTK
 * @brief Class for writing VTK Structured Points files for density and velocity fields.
 */
class VTK {
public:
    /**
     * @brief Constructor.
     * @param nx number of grid points in x.
     * @param ny number of grid points in y.
     * @param nz number of grid points in z.
     * @param dx grid spacing.
     */
    VTK(unsigned int nx, unsigned int ny, unsigned int nz, double dx);

    /**
     * @brief Write one VTK frame.
     * @param filename output file name.
     * @param iteration current iteration number.
     * @param lbm reference to LBM solver.
     */
    void writeVTKFrame(const std::string& filename,
                       int iteration,
                       const TDLBM& lbm) const;

private:
    unsigned int nx_, ny_, nz_; ///< grid dimensions
    double dx_;                 ///< grid spacing
};
