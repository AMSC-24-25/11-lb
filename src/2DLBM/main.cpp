#include "lib/Lattice.hpp"
#include "lib/LDRIVEN.hpp"
#include <iostream>
#include <string>
#include <filesystem>

#ifdef USE_OPENMP
  #include <omp.h>
#endif

int main(int argc, char* argv[]){

    if (argc < 3 || argc > 8) {
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
    const double u_lid = 0.1;
                      

    if(argc == 7 ){
        use_tunnel = true;
    }

    if(use_tunnel == true){
        Lattice lattice(nx,ny,u_lid,Re);
        lattice.simulate();
    }else{
         LDRIVEN cavity(nx,ny,u_lid,Re);
         cavity.evolution(Steps);
    }
    
    

    return 0;
}