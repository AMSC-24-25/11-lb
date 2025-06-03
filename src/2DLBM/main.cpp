#include "Lattice.hpp"
#include <iostream>
#include <string>
#include <filesystem>

int main(int argc, char* argv[]){

    if (argc < 4 || argc > 7) {
        std::cerr << "Usage: " << argv[0]
                  << " NX Nsteps Re [maskType maskSize] [outputDir]\n";
        return EXIT_FAILURE;
    }

    unsigned int mesh = std::stoul(argv[1]);              
    unsigned int Steps = std::stoul(argv[2]);              
    double Re = std::stod(argv[3]);               
   
    bool        useMask   = false;                            
    std::string maskType  = "";                             
    double      maskSize  = 0.0;                              
    std::string outDir    = "./output";                     

    if (argc == 5) {
       
        outDir = argv[4];
    } else if (argc == 6) {
        
        useMask  = true;
        maskType = argv[4];
        maskSize = std::stod(argv[5]);
    } else if (argc == 7) {
        
        useMask  = true;
        maskType = argv[4];
        maskSize = std::stod(argv[5]);
        outDir   = argv[6];
    }

    // Ensure output directory exists
    //std::filesystem::create_directories(outDir);
    
    Lattice lattice(mesh, Steps, Re,
                    useMask, maskType, maskSize,
                    outDir);
    lattice.simulate();

    return 0;
}