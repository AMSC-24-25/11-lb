#!/usr/bin/env bash
set -euo pipefail

print_help() {
  cat <<EOF
Usage: $0 [OPTIONS]

Algorithm (choose one):
  --3dLbm                     Run 3-D LBM solver
  --2dLbm                     Run 2-D LBM solver

Common options:
  -m,  --mesh NX,NY,NZ        Mesh dimension (single value for cubic mesh)
  -s,  --steps N              Number of time steps to simulate
  -r,  --re RE                Reynolds number
  -d,  --dir PATH             Output directory (default: ./output)
  -omp,--Openmp               Enable OpenMP
  -h,  --help                 Show this help

  2-D specific:
  --mask TYPE                 Enable obstacle mask of TYPE {circle|rect|airfoil}
  --mask-size VAL             Characteristic size of the mask (required with --mask)

  3-D specific:
  -itf,--iters-per-frame N    Iterations per frame for output (default: 10)


EOF
}

# Defaults
ALGORITHM=""
MESH_SIZE=""
TIME_STEPS=""
REYNOLDS_NUMBER=""
ITER_PER_FRAME="10"
OUTPUT_DIR=""
USE_CUSTOM_DIR=false
USE_OPENMP=false


USE_MASK=false
MASK_TYPE=""
MASK_SIZE=""

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --3dLbm)   ALGORITHM="3dLbm"; shift ;;
    --2dLbm)   ALGORITHM="2dLbm"; shift ;;
    -m|--mesh) MESH_SIZE="$2"; shift 2 ;;
    -s|--steps) TIME_STEPS="$2"; shift 2 ;;
    -r|--re|--reynolds) REYNOLDS_NUMBER="$2"; shift 2 ;;
    -d|--dir)  OUTPUT_DIR="$2"; USE_CUSTOM_DIR=true; shift 2 ;;
    -itf|--iters-per-frame|--frameite) ITER_PER_FRAME="$2"; shift 2 ;;
    -omp|--Openmp) USE_OPENMP=true; shift ;;
    --mask)  USE_MASK=true; MASK_TYPE="$2"; shift 2 ;;
    --mask-size) MASK_SIZE="$2"; shift 2 ;;
    -h|--help) print_help; exit 0 ;;
    *) echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

# Validate common
if [[ -z $ALGORITHM || -z $MESH_SIZE || -z $TIME_STEPS || -z $REYNOLDS_NUMBER ]]; then
  echo "Error: must specify algorithm, mesh, steps, and Reynolds number." >&2
  print_help; exit 1
fi

# Validate mask
if [[ $ALGORITHM == "2dLbm" && $USE_MASK == true && -z $MASK_SIZE ]]; then
  echo "Error: --mask-size is required when --mask is specified." >&2; exit 1
fi

# Build
BUILD_DIR="build"
mkdir -p "$BUILD_DIR"
pushd "$BUILD_DIR" > /dev/null
CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=Release"
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
  if [[ $USE_MASK == true ]]; then
    if [[ $USE_CUSTOM_DIR == true ]]; then
      exec "${BUILD_DIR}/11-LB" \
           "$MESH_SIZE" "$TIME_STEPS" "$REYNOLDS_NUMBER" \
           "$MASK_TYPE" "$MASK_SIZE" "$OUTPUT_DIR"
    else
      exec "${BUILD_DIR}/11-LB" \
           "$MESH_SIZE" "$TIME_STEPS" "$REYNOLDS_NUMBER"  \
           "$MASK_TYPE" "$MASK_SIZE"
    fi
  else
    if [[ $USE_CUSTOM_DIR == true ]]; then
      exec "${BUILD_DIR}/11-LB" \ 
           "$MESH_SIZE" "$TIME_STEPS" "$REYNOLDS_NUMBER" \ 
           "$OUTPUT_DIR"
    else
      exec "${BUILD_DIR}/11-LB" \ 
           "$MESH_SIZE" "$TIME_STEPS" "$REYNOLDS_NUMBER" 
    fi
  fi
else
  echo "Unsupported algorithm: $ALGORITHM" >&2
  print_help
  exit 1
fi
