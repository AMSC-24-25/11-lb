#include "lib/Lattice.hpp"
#include "lib/LDRIVEN.hpp"
#include <iostream>
#include <string>
#include <filesystem>

unsigned int nx, ny;

int main(int argc, char* argv[]){

    if (argc < 4 || argc > 7) {
        std::cerr << "Usage: " << argv[0]
                  << " <mesh_size> <time_steps> <reynolds> [output_dir]\n";
        return EXIT_FAILURE;
    }

    unsigned int mesh  = std::stoi(argv[1]);
    int Steps = std::stoi(argv[2]);
    double Re    = std::stod(argv[3]);
    int ITERATIONS_PER_FRAME = std::stoi(argv[4]);
    bool use_tunnel= false;
    const double u_lid = 0.1;
    nx = mesh; ny = mesh;                   

    if(argc == 6 ){
        use_tunnel = true;
    }

    if(use_tunnel == true){
        Lattice lattice;
        lattice.simulate();
    }else{
         LDRIVEN cavity(nx,ny,u_lid,Re);
         cavity.evolution(Steps);
    }
    
    

    return 0;
}