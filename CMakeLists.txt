cmake_minimum_required(VERSION 3.14)

project(11-LB)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

option(WITH_OPENMP "Enable OpenMP parallelization" OFF)
option(ENABLE_3D "Build 3-D LBM solver" OFF)
option(ENABLE_2D "Build 2-D LBM solver" OFF)

if(NOT ENABLE_3D AND NOT ENABLE_2D)
    message(FATAL_ERROR "Enable at least one solver")
endif()

if(WITH_OPENMP)
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    link_directories("/usr/local/opt/libomp/lib")
endif()

include_directories(src)
include_directories(src/common)

if(ENABLE_3D)
    file(GLOB SOURCES_3D
         src/3DLBM/*.cpp
         src/3DLBM/lib/*.cpp
         src/common/*.cpp)
    add_executable(11-LB ${SOURCES_3D})
    if(WITH_OPENMP)
        target_link_libraries(11-LB PRIVATE ${OpenMP_omp_LIBRARY})
        target_compile_definitions(11-LB PRIVATE USE_OPENMP)
    endif()
endif()

if(ENABLE_2D)
    file(GLOB SOURCES_2D
         src/2DLBM/*.cpp
         src/2DLBM/lib/*.cpp
         src/common/*.cpp)
    add_executable(11-LB ${SOURCES_2D})
    if(WITH_OPENMP)
        target_link_libraries(11-LB PRIVATE ${OpenMP_omp_LIBRARY})
        target_compile_definitions(11-LB PRIVATE USE_OPENMP)
    endif()
endif()
