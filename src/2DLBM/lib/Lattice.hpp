#ifndef LATTICE_HPP
#define LATTICE_HPP

#include "Node.hpp"
#include "Auxiliary.cpp"
#include <cmath>

class Lattice
{
    public:
        Lattice(unsigned NX, unsigned steps, double Re,
        bool useMask, const std::string &maskType,
        double maskSize, const std::string &outDir);
        void simulate();
        Matrix<Node> node_matrix;

    private:
        int sigma;
        double omega_P, omega_M, size;
        std::vector<double> boundary_velocity;
        Matrix<bool> obstacles;
        int currentStep = 0;
        int maxSteps;
        int ITERATIONS_PER_FRAME;
        int ITERATIONS_PER_PROGRESS_UPDATE;

        std::string outDir_;    
        bool        useMask_;   
        std::string maskType_;  
        double      maskSize_;  


        

};


#endif
