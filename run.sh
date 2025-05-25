#!/usr/bin/env bash
set -euo pipefail

print_help() {
  cat <<EOF
Usage: $0 [OPTIONS]

OPTIONS:
  --3dLbm                      Select the 3D LBM algorithm
  -m, --mesh NX                Mesh dimension (single value for cubic mesh)
  -s, --steps N                Number of time steps to simulate
  -r, --re RE                  Reynolds number
  -d, --dir PATH               Output directory for simulation results (default: ./output)
  -h, --help                   Display this help message
  -itf, --iters-per-frame N    iteration per frame
  -omp, --Openmp               enable openmp
EOF
}

# Default settings
OUTPUT_DIR=""
USE_CUSTOM_DIR=false
ALGORITHM=""
MESH_SIZE=""
TIME_STEPS=""
REYNOLDS_NUMBER=""
ITER_PER_FRAME="10"
USE_OPENMP=false


# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --3dLbm)
      ALGORITHM="3dLbm"; shift ;;
      -omp|Openmp)
      USE_OPENMP=true; shift ;;
    -m|--mesh)
      MESH_SIZE="$2"; shift 2 ;;
    -s|--steps)
      TIME_STEPS="$2"; shift 2 ;;
    -r|--re|--reynolds)
      REYNOLDS_NUMBER="$2"; shift 2 ;;
    -d|--dir)
      OUTPUT_DIR="$2"
      USE_CUSTOM_DIR=true
      shift 2 ;;
      -itf|--frameite)
      ITERATION_PER_FRAME="$2"
      shift 2;;
    -h|--help)
      print_help; exit 0 ;;
    *)
      echo "Unknown option: $1"; print_help; exit 1 ;;
  esac
done

# Validate required parameters
if [[ -z $ALGORITHM || -z $MESH_SIZE || -z $TIME_STEPS || -z $REYNOLDS_NUMBER ]]; then
  echo "Error: algorithm, mesh size, time steps, and Reynolds number must be specified." >&2
  print_help
  exit 1
fi


# Build directory
BUILD_DIR="./build"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"
if [[ "$USE_OPENMP" == "true" ]]; then
  cmake .. -DCMAKE_BUILD_TYPE=Release -DWITH_OPENMP=ON
 else
  cmake .. -DCMAKE_BUILD_TYPE=Release -DWITH_OPENMP=OFF
fi
make -j
cd -

# Execute simulation
case "$ALGORITHM" in
  3dLbm)
    echo "Running simulation..."
    if [[ "$USE_CUSTOM_DIR" == true ]]; then
      "$BUILD_DIR/11-LB" "$MESH_SIZE" "$TIME_STEPS" "$REYNOLDS_NUMBER" "$OUTPUT_DIR" "$ITER_PER_FRAME"
    else
      "$BUILD_DIR/11-LB" "$MESH_SIZE" "$TIME_STEPS" "$REYNOLDS_NUMBER" "$ITER_PER_FRAME"
    fi
    ;;
  *)
    echo "Unsupported algorithm: $ALGORITHM" >&2
    print_help
    exit 1
    ;;
 esac
