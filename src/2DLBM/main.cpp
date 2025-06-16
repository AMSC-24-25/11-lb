#include "lib/Lattice.hpp"
#include "lib/LDRIVEN.hpp"
#include <iostream>
#include <string>
#include <filesystem>

#ifdef USE_OPENMP
  #include <omp.h>
#endif

int main(int argc, char* argv[]){

    if (argc < 3 || argc > 10) {
        std::cerr << "Usage: " << argv[0]
                  << " <mesh_size> <time_steps> <reynolds> [output_dir]\n";
        return EXIT_FAILURE;
    }

    #ifdef USE_OPENMP
    #pragma omp parallel
    {
        int n_threads = omp_get_num_threads();

        #pragma omp single
        std::cout << ">>> OpenMP active with " 
                  << n_threads 
                  << " thread" <<"\n";

    }
    #else
    std::cout << "OpenMP do not work";
    #endif

    unsigned int nx  = std::stoi(argv[1]);
    unsigned int ny  = std::stoi(argv[2]);
    int Steps = std::stoi(argv[3]);
    double Re  = std::stod(argv[4]);
    int ITERATIONS_PER_FRAME = std::stoi(argv[5]);
    bool use_tunnel= false;
    std::string ArgDir=argv[6];     
    double u = std::stod(argv[7]);      
    if(argc == 9 ){
        use_tunnel = true;
    }

    if(use_tunnel == true){
        Lattice lattice(nx,ny,u,Re,ArgDir);
        lattice.simulate(Steps,ITERATIONS_PER_FRAME);
    }else{
         LDRIVEN cavity(nx,ny,u,Re,ArgDir);
         cavity.simulate(Steps,ITERATIONS_PER_FRAME);
    }
    
    

    return 0;
}