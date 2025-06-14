#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "Node.hpp"
#include "Auxiliary.hpp"
#include <cmath>
#include <string>


class Lattice
{
    public:
        Lattice(unsigned int nx, unsigned int ny, double u_lid, double Re, std::string outdir);
        void simulate(int max_steps, int iter_per_frame);
        Matrix<Node> node_matrix;
        unsigned int object_count = 0;

    private:
        int sigma;
        double omega_P, omega_M, size;
        std::vector<double> boundary_velocity;
        Matrix<bool> obstacles;
        int maxSteps;
        int ITERATIONS_PER_FRAME;
        int ITERATIONS_PER_PROGRESS_UPDATE;
        int NX;
        int NY;
        double Re;
        double u_lid;
        std::string outdir;

};

#endif
