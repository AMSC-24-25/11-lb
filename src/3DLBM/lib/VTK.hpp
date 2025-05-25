#pragma once

#include <string>
#include "3DLBM.hpp"

/**
 * @class VTKWriter
 * @brief Class for writing VTK Structured Points files for density and velocity fields.
 *
 * Provides a method writeVTKFrame whose implementation matches the original
 * standalone writeVTKFrame exactly, adapted as a class method.
 */
class VTK {
public:
    /**
     * @brief Constructs a VTKWriter with grid dimensions and spacing.
     * @param nx Number of grid points in the x-direction.
     * @param ny Number of grid points in the y-direction.
     * @param nz Number of grid points in the z-direction.
     * @param dx Grid spacing (uniform in all directions).
     */
    VTK(unsigned int nx, unsigned int ny, unsigned int nz, double dx);

    /**
     * @brief Writes density and velocity fields to a VTK Structured Points file.
     * @param filename Full path of the output VTK file (including .vtk extension).
     * @param iteration Iteration number to include in the file header.
     * @param lbm Reference to the TDLBM simulation object providing data.
     *
     * Implementation matches the original writeVTKFrame function exactly,
     * modified only to operate as a class method.
     */
    void writeVTKFrame(const std::string& filename,
                       int iteration,
                       const TDLBM& lbm) const;

private:
    unsigned int nx_;  ///< number of grid points in x
    unsigned int ny_;  ///< number of grid points in y
    unsigned int nz_;  ///< number of grid points in z
    double dx_;        ///< grid spacing in all directions
};