#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "Node.hpp"
#include "Auxiliary.cpp"
#include <cmath>
#include <omp.h>

class Lattice
{
    public:
        Lattice();
        void simulate();
        Matrix<Node> node_matrix;
        unsigned int object_count = 0;

    private:
        int sigma;
        double omega_P, omega_M, size;
        std::vector<double> boundary_velocity;
        Matrix<bool> obstacles;
        int currentStep = 0;
        int maxSteps;
        int ITERATIONS_PER_FRAME;
        int ITERATIONS_PER_PROGRESS_UPDATE;

};

#endif
