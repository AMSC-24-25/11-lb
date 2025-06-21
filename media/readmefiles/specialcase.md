## 0. **Special Case** (in case of problem)

### macOS and OpenMP:
Running the code natively on macOS requires the LLVM‐based Clang toolchain bundled with Xcode and the **libomp** runtime installed via **Homebrew**:

```bash
brew install libomp      # installs the OpenMP runtime
```

Because Apple’s Clang does not automatically locate Homebrew’s libomp, you must replace the standard OpenMP discovery block in `CMakeLists.txt`:

```cmake
Default discovery (may fail on macOS)
if (WITH_OPENMP)
    find_package(OpenMP REQUIRED)
    if (OpenMP_CXX_FOUND)
        message(STATUS "Found OpenMP")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()
```

with the explicit linkage shown below:

```cmake
Explicit linkage for Homebrew‐installed libomp
if (WITH_OPENMP)
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY   "/usr/local/opt/libomp/lib/libomp.dylib")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    link_directories("/usr/local/opt/libomp/lib")
endif()
```

> **Why do the paths sometimes differ?**  
> Homebrew installs libraries under the prefix returned by `brew --prefix`.  
> - On Intel Macs this is typically `/usr/local/opt/...`.  
> - On Apple Silicon it is `/opt/homebrew/opt/...`.  
> Adjust the include and library paths accordingly.
---
<br>

### CUDA:
To execute the CUDA implementation you need an **NVIDIA GPU** and the **CUDA Toolkit** already installed. Occasionally, certain CUDA versions fail to auto-detect the correct compute capability, so you must set it manually.

1. **Identify your GPU and toolkit version**

```bash
nvidia-smi        # prints GPU model (e.g., Tesla T4 → CC 7.5)
nvcc --version    # prints CUDA Toolkit version
```


 2. **Edit `CMakeLists.txt`**

Replace the automatic architecture selection

```cmake
set(CMAKE_CUDA_ARCHITECTURES AUTO)
```

with the exact compute capability of your GPU, e.g. for a Tesla T4:

```cmake
set(CMAKE_CUDA_ARCHITECTURES 75)   # CC 7.5
```

> If multiple GPUs are present, list each architecture separated by semicolons, e.g. `set(CMAKE_CUDA_ARCHITECTURES 75;86)`.
> 

---
<br>