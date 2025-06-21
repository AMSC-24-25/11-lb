cmake_minimum_required(VERSION 3.14)

project(11-LB)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

option(WITH_OPENMP "Enable OpenMP parallelization" OFF)
option(ENABLE_3D "Build 3-D LBM solver" OFF)
option(ENABLE_2D "Build 2-D LBM solver" OFF)
option(ENABLE_CUDA   "Build CUDA solver" OFF) 

if(NOT ENABLE_3D AND NOT ENABLE_2D AND NOT ENABLE_CUDA)
    message(FATAL_ERROR "Enable at least one solver")
endif()

if(WITH_OPENMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_CXX_FOUND)
        message(STATUS "Found OpenMP")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

if(ENABLE_CUDA)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    set(CMAKE_CUDA_ARCHITECTURES AUTO)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

include_directories(src)
include_directories(src/common)

if(ENABLE_3D)
    file(GLOB SOURCES_3D
         src/3DLBM/*.cpp
         src/3DLBM/include/*.cpp
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
         src/2DLBM/include/*.cpp
         src/common/*.cpp)
    add_executable(11-LB ${SOURCES_2D})
    if(WITH_OPENMP)
        target_link_libraries(11-LB PRIVATE ${OpenMP_omp_LIBRARY})
        target_compile_definitions(11-LB PRIVATE USE_OPENMP)
    endif()
endif()

if(ENABLE_CUDA)

   add_executable(11-LB ${CMAKE_SOURCE_DIR}/src/cuda/lbm_2D_liddrivencavity.cu)
   target_compile_options(11-LB PRIVATE
   $<$<COMPILE_LANGUAGE:CUDA>:-O3 -use_fast_math -arch=sm_75>)
   set_target_properties(11-LB PROPERTIES
   CUDA_SEPARABLE_COMPILATION ON)
   target_link_libraries(11-LB PRIVATE CUDA::cudart)
   target_compile_definitions(11-LB PRIVATE USE_CUDA)
endif()
