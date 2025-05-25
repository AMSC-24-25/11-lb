cmake_minimum_required(VERSION 3.14)

project(11-LB)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-march=native")

option(WITH_OPENMP "Enable OpenMP parallelization" OFF)

if(WITH_OPENMP)
set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
set(OpenMP_CXX_LIB_NAMES "omp")
set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
link_directories("/usr/local/opt/libomp/lib")
add_compile_definitions(USE_OPENMP)
endif()

set(SOURCES
   src/3DLBM/main.cpp
    src/3DLBM/lib/3DLBM.cpp
     src/3DLBM/lib/VTK.cpp
)

set(HEADERS    
    src/3DLBM/lib/3DLBM.hpp 
    src/3DLBM/lib/VTK.hpp 
)

include_directories(src)

add_executable(11-LB ${SOURCES})

if(WITH_OPENMP)
    target_link_libraries(11-LB PRIVATE ${OpenMP_omp_LIBRARY})
endif()