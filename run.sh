#!/usr/bin/env bash
set -euo pipefail

print_help() {
  cat <<EOF
Usage: $0 [OPTIONS]

Algorithm (choose one):
  --3dLbm                     Run 3-D LBM solver
  --2dLbm                     Run 2-D LBM solver

Common options:
  -m,  --mesh NX        Mesh dimension (single value for cubic mesh)
  -s,  --steps N              Number of time steps to simulate
  -r,  --re RE                Reynolds number
  -d,  --dir PATH             Output directory (default: ./output)
  -omp,--Openmp               Enable OpenMP
  -h,  --help                 Show this help
  -itf,--iters-per-frame N    Iterations per frame for output (default: 10)

  2D specific
   -tunnel  use the wind tunnel instead of lid driven cavity
   -my      mesh dimensio on y 


EOF
}

# Defaults
ALGORITHM=""
MESH_SIZE=""
MESH_SIZE_Y=""
TIME_STEPS=""
REYNOLDS_NUMBER=""
ITER_PER_FRAME="10"
OUTPUT_DIR=""
USE_CUSTOM_DIR=false
USE_OPENMP=false
USE_TUNNEL=false




# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --3dLbm)   ALGORITHM="3dLbm"; shift ;;
    --2dLbm)   ALGORITHM="2dLbm"; shift ;;
    -tunnel) USE_TUNNEL=true; shift 1;;
    -my) MESH_SIZE_Y="$2"; shift 2;;
    -m|--mesh) MESH_SIZE="$2"; shift 2 ;;
    -s|--steps) TIME_STEPS="$2"; shift 2 ;;
    -r|--re|--reynolds) REYNOLDS_NUMBER="$2"; shift 2 ;;
    -d|--dir)  OUTPUT_DIR="$2"; USE_CUSTOM_DIR=true; shift 2 ;;
    -itf|--iters-per-frame|--frameite) ITER_PER_FRAME="$2"; shift 2 ;;
    -omp|--Openmp) USE_OPENMP=true; shift ;;
    
    -h|--help) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

MESH_SIZE_Y="${MESH_SIZE_Y:-$MESH_SIZE}"

# Validate common
if [[ -z $ALGORITHM || -z $MESH_SIZE || -z $TIME_STEPS || -z $REYNOLDS_NUMBER ]]; then
  echo "Error: must specify algorithm, mesh, steps, and Reynolds number." >&2
  print_help; exit 1
fi


# Build
BUILD_DIR="build"
mkdir -p "$BUILD_DIR"
pushd "$BUILD_DIR" > /dev/null
CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS=-fsanitize=address"
[[ $USE_OPENMP == true ]] && CMAKE_FLAGS+=" -DWITH_OPENMP=ON" || CMAKE_FLAGS+=" -DWITH_OPENMP=OFF"
if [[ $ALGORITHM == "3dLbm" ]]; then
  CMAKE_FLAGS+=" -DENABLE_3D=ON -DENABLE_2D=OFF"
elif [[ $ALGORITHM == "2dLbm" ]]; then
  CMAKE_FLAGS+=" -DENABLE_3D=OFF -DENABLE_2D=ON"
else
  echo "Choose --3dLbm or --2dLbm." >&2; exit 1
fi
cmake .. ${CMAKE_FLAGS}
make -j"$(nproc)"
popd > /dev/null

# Run
if [[ $ALGORITHM == "3dLbm" ]]; then

  if [[ $USE_CUSTOM_DIR == true ]]; then
    exec "${BUILD_DIR}/11-LB" "$MESH_SIZE" "$TIME_STEPS" "$REYNOLDS_NUMBER" "$OUTPUT_DIR" "$ITER_PER_FRAME"
  else
    exec "${BUILD_DIR}/11-LB" "$MESH_SIZE" "$TIME_STEPS" "$REYNOLDS_NUMBER" "$ITER_PER_FRAME"
  fi

elif [[ $ALGORITHM == "2dLbm" ]]; then

  if [[ $USE_TUNNEL == true ]]; then
    exec "${BUILD_DIR}/11-LB" "$MESH_SIZE" "$MESH_SIZE_Y" "$TIME_STEPS" "$REYNOLDS_NUMBER" "$ITER_PER_FRAME" "$USE_TUNNEL"
  else
    exec "${BUILD_DIR}/11-LB" "$MESH_SIZE" "$MESH_SIZE_Y" "$TIME_STEPS" "$REYNOLDS_NUMBER" "$ITER_PER_FRAME" 
  fi  
else
  echo "Unsupported algorithm: $ALGORITHM" >&2
  print_help
  exit 1
fi
