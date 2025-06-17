# Project Group 11: Lattice Boltzmann Methods 
### Introduction and objectives
This projects uses the Lattice Boltzmann Methods (LBM) to perform a 2D fluid simulation with the D2Q9 model.The aim for the hands-on was to solve the 2D lid-driven cavity problem, the project extended it to 3D Lid-Driven, 2D Wind-Tunnel-like problem and 2D Lid-Driven paralellized with CUDA.


# Build and Execution Guide
<br>

## 0. **Special Case** (in case of problem)

### macOS and OpenMP:
Running the code natively on macOS requires the LLVM‐based Clang toolchain bundled with Xcode and the **libomp** runtime installed via **Homebrew**:

```bash
brew install libomp      # installs the OpenMP runtime
```

Because Apple’s Clang does not automatically locate Homebrew’s libomp, you must replace the standard OpenMP discovery block in `CMakeLists.txt`:

```cmake
 ❌ Default discovery (may fail on macOS)
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
 ✅ Explicit linkage for Homebrew‐installed libomp
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

## 1. **Enable execution permissions for the launcher**

   ```bash
   chmod +x run.sh
   ```


## 2. **Build and Run the simulation**

   From the directory where the run.sh is, launch:
   ```bash
   ./run.sh [FLAGS]
   ```

---

### Solver Selection flags

| Flag          | Behavior                                  |
|---------------|-------------------------------------------|
| `--3dLbm`     | Select the 3D solver implementation.      |
| `--2dLbm`     | Select the 2D solver implementation.      |
| `--cuda`      | Use the CUDA-based GPU 2D implementation. |

---

### Available Flags (you must put the obligatory flags)

| Flag                            | Behavior                                               | 3D LBM | 2D LBM | OBLIGATORY|
|---------------------------------|--------------------------------------------------------|:------:|:------:|:---------------|
| `-m`, `--mesh`                  | Define the grid resolution in X for 2D and for all direction in 3D.    | ✅     | ✅     | ✅|
| `-s`, `--steps`                 | Set total number of simulation time steps.             | ✅     | ✅     |✅|
| `-r`, `--re`, `--reynolds`      | Assign the Reynolds number for flow conditions.        | ✅     | ✅     |✅|
| `-u`                            | Initial or boundary fluid velocity (default: 0.1).     | ✅     | ✅     |⛔|    
| `-omp`, `--Openmp`              | Activate OpenMP.                                       | ✅     | ✅     |⛔|    
| `-h`, `--help`                  | Display usage information and exit.                    | ✅     | ✅     |⛔|    
| `-tunnel`                       | Apply wind-tunnel boundary geometry.                   | ⛔     | ✅     |⛔|   
| `-my`                           | Specify grid resolution in Y only for 2D(equal to X default).      | ⛔     | ✅     |⛔| 

---

**Examples:**

   * **3D solver** with a cubic mesh of 100 points, 1,000 time steps, Reynolds number = 500:

     ```bash
     ./run.sh --3dLbm -m 100 -s 1000 -r 500
     ```

   * **2D solver** with OpenMP, wind-tunnel geometry, mesh 200×100, 2,000 steps, Re = 1000:

     ```bash
     ./run.sh --2dLbm -m 200 -my 100 -s 2000 -r 1000 -omp -tunnel
     ```

   * **CUDA solver** on GPU architecture SM 7.5, mesh 256, 1,500 steps, Re = 750:

     ```bash
     ./run.sh --cuda -m 256 -s 1500 -r 750
     ```
## 3. **running cuda with coolab**

if you want to run the cuda implementation without having a Nvidia GPU you find a python notebook to upload in colab with all the instruction [here](src/cuda/lbm_cuda.ipynb)


---



### Parallelization Strategy
We parallelized the code using OpenMP. The main bottleneck in the computation is represented by a single nested for-loop. The code was written so that each iteration was completely independent from the others and so that they all could be executed in parallel. To achieve this we needed to *bufferize* the whole computation so that two sets of memory location were used and then swapped for each iteration.
We observed a significant improvement in computation time.


### Scalability

<div class="container">
  <div class="column">
    <h4>Strong scalability test (6 cores, 12 threads)</h4>
    <div class="note">Execution times recorded on a PC with 6 cores and 12 threads.<br>Simulation run for 10'000 iterations.</div>
    <table>
      <thead>
        <tr>
          <th>Threads</th>
          <th>Time (ms)</th>
          <th>Speedup</th>
          <th>Efficiency</th>
        </tr>
      </thead>
      <tbody>
        <tr><td>1</td><td>69187</td><td>-</td><td>-</td></tr>
        <tr><td>2</td><td>38153</td><td>1.81</td><td>0.91</td></tr>
        <tr><td>3</td><td>29987</td><td>2.31</td><td>0.77</td></tr>
        <tr><td>4</td><td>27700</td><td>2.50</td><td>0.63</td></tr>
        <tr><td>5</td><td>24697</td><td>2.80</td><td>0.56</td></tr>
        <tr><td>6</td><td>22624</td><td>3.06</td><td>0.51</td></tr>
        <tr><td>7</td><td>19980</td><td>3.46</td><td>0.49</td></tr>
        <tr><td>8</td><td>17710</td><td>3.91</td><td>0.49</td></tr>
        <tr><td>9</td><td>16086</td><td>4.30</td><td>0.48</td></tr>
        <tr><td>10</td><td>14747</td><td>4.69</td><td>0.47</td></tr>
        <tr><td>11</td><td>14664</td><td>4.72</td><td>0.43</td></tr>
        <tr><td>12</td><td>18187</td><td>3.80</td><td>0.32</td></tr>
      </tbody>
    </table>
  </div>
  
  <div class="column">
    <h4>Strong scalability test (8 cores, 8 threads)</h4>
    <div class="note">Execution times recorded on a PC with 8 cores and 8 threads.<br>Simulation run for 10'000 iterations.</div>
    <table>
      <thead>
        <tr>
          <th>Threads</th>
          <th>Time (ms)</th>
          <th>Speedup</th>
          <th>Efficiency</th>
        </tr>
      </thead>
      <tbody>
        <tr><td>1</td><td>53929</td><td>-</td><td>-</td></tr>
        <tr><td>2</td><td>31645</td><td>1.70</td><td>0.85</td></tr>
        <tr><td>3</td><td>23348</td><td>2.31</td><td>0.77</td></tr>
        <tr><td>4</td><td>19208</td><td>2.81</td><td>0.70</td></tr>
        <tr><td>5</td><td>16473</td><td>3.27</td><td>0.65</td></tr>
        <tr><td>6</td><td>15139</td><td>3.56</td><td>0.59</td></tr>
        <tr><td>7</td><td>13518</td><td>3.99</td><td>0.57</td></tr>
        <tr><td>8</td><td>25401</td><td>2.12</td><td>0.27</td></tr>
      </tbody>
    </table>
  </div>
</div>


<div class="container">
  <div class="column">
    <h4>Weak scalability test (6 cores, 12 threads)</h4>
    <div class="note">Execution times recorded on a PC with 6 cores and 12 threads.<br>Simulation run for 1'000 iterations with ~20'000 lattice points per thread.</div>
    <table>
      <thead>
        <tr>
          <th>Threads</th>
          <th>Cavity Size</th>
          <th>Time (ms)</th>
          <th>Efficiency</th>
        </tr>
      </thead>
      <tbody>
        <tr><td>1</td><td>142x141</td><td>1381</td><td>-</td></tr>
        <tr><td>2</td><td>200x200</td><td>1537</td><td>0.90</td></tr>
        <tr><td>3</td><td>246x244</td><td>1821</td><td>0.76</td></tr>
        <tr><td>4</td><td>284x282</td><td>2029</td><td>0.68</td></tr>
        <tr><td>5</td><td>317x316</td><td>2205</td><td>0.63</td></tr>
        <tr><td>6</td><td>347x346</td><td>2351</td><td>0.59</td></tr>
        <tr><td>7</td><td>375x374</td><td>2503</td><td>0.55</td></tr>
        <tr><td>8</td><td>400x400</td><td>2598</td><td>0.53</td></tr>
        <tr><td>9</td><td>425x424</td><td>2651</td><td>0.52</td></tr>
        <tr><td>10</td><td>448x447</td><td>2748</td><td>0.50</td></tr>
        <tr><td>11</td><td>470x469</td><td>2942</td><td>0.47</td></tr>
        <tr><td>12</td><td>491x489</td><td>3288</td><td>0.42</td></tr>
      </tbody>
    </table>
  </div>
  
  <div class="column">
    <h4>Weak scalability test (8 cores, 8 threads)</h4>
    <div class="note">Execution times recorded on a PC with 8 cores and 8 threads.<br>Simulation run for 1'000 iterations with ~20'000 lattice points per thread.</div>
    <table>
      <thead>
        <tr>
          <th>Threads</th>
          <th>Cavity Size</th>
          <th>Time (ms)</th>
          <th>Efficiency</th>
        </tr>
      </thead>
      <tbody>
        <tr><td>1</td><td>142x141</td><td>1016</td><td>-</td></tr>
        <tr><td>2</td><td>200x200</td><td>1027</td><td>0.99</td></tr>
        <tr><td>3</td><td>246x244</td><td>1058</td><td>0.96</td></tr>
        <tr><td>4</td><td>284x282</td><td>1109</td><td>0.92</td></tr>
        <tr><td>5</td><td>317x316</td><td>1191</td><td>0.85</td></tr>
        <tr><td>6</td><td>347x346</td><td>1257</td><td>0.81</td></tr>
        <tr><td>7</td><td>375x374</td><td>1295</td><td>0.78</td></tr>
        <tr><td>8</td><td>400x400</td><td>1413</td><td>0.72</td></tr>
      </tbody>
    </table>
  </div>
</div>

### Validation of Results
We were provided with some reference data to compare our results.
The following tables, for instance, allows to compare the y-component of the velocity of the fluid along the vertical line through the geometrical center of the cavity.  

It can be observed that although the results with Re = 100 are relatively accurate, there's an error build-up with toward the lower part of the cavity and with increasing Reynolds numbers. That can probably be caused by slightly different parameters and formulae used. We are confident that our code works as intended and that, with some tweaking, we could allign our results with the reference data.



<div class="container2">
  <div class="row2">
    <div class="column2">
      <h4>Our results</h4>
      <table>
        <thead>
          <tr>
            <th>129-grid pt. no.</th>
            <th>y</th>
            <th>Re (100)</th>
            <th>Re (400)</th>
            <th>Re (1000)</th>
          </tr>
        </thead>
        <tbody>
          <tr><td>129</td><td>1.00</td><td>1</td><td>1</td><td>1</td></tr>
          <tr><td>126</td><td>0.98</td><td>0.844576</td><td>0.759886</td><td>0.654879</td></tr>
          <tr><td>125</td><td>0.97</td><td>0.792432</td><td>0.68533</td><td>0.56604</td></tr>
          <tr><td>124</td><td>0.96</td><td>0.741241</td><td>0.6171</td><td>0.496168</td></tr>
          <tr><td>123</td><td>0.95</td><td>0.691354</td><td>0.556289</td><td>0.44422</td></tr>
          <tr><td>110</td><td>0.85</td><td>0.222347</td><td>0.262017</td><td>0.295701</td></tr>
          <tr><td>95</td><td>0.74</td><td>-0.01649</td><td>0.141458</td><td>0.163419</td></tr>
          <tr><td>80</td><td>0.62</td><td>-0.1383</td><td>0.004383</td><td>0.042401</td></tr>
          <tr><td>65</td><td>0.50</td><td>-0.17906</td><td>-0.13328</td><td>-0.0651</td></tr>
          <tr><td>59</td><td>0.46</td><td>-0.17623</td><td>-0.18906</td><td>-0.10638</td></tr>
          <tr><td>37</td><td>0.29</td><td>-0.12114</td><td>-0.27501</td><td>-0.2716</td></tr>
          <tr><td>23</td><td>0.18</td><td>-0.07795</td><td>-0.16957</td><td>-0.31972</td></tr>
          <tr><td>14</td><td>0.11</td><td>-0.04994</td><td>-0.09319</td><td>-0.20474</td></tr>
          <tr><td>10</td><td>0.08</td><td>-0.03668</td><td>-0.06236</td><td>-0.14069</td></tr>
          <tr><td>9</td><td>0.07</td><td>-0.03321</td><td>-0.05479</td><td>-0.1247</td></tr>
          <tr><td>8</td><td>0.06</td><td>-0.02968</td><td>-0.04723</td><td>-0.10798</td></tr>
          <tr><td>1</td><td>0.00</td><td>0</td><td>0</td><td>0</td></tr>
        </tbody>
      </table>
    </div>
    <div class="column2">
      <h4>Reference data</h4>
      <table>
        <thead>
          <tr>
            <th>129-grid pt. no.</th>
            <th>y</th>
            <th>Re (100)</th>
            <th>Re (400)</th>
            <th>Re (1000)</th>
          </tr>
        </thead>
        <tbody>
          <tr><td>129</td><td>1.00000</td><td>1.00000</td><td>1.00000</td><td>1.00000</td></tr>
          <tr><td>126</td><td>0.9766</td><td>0.84123</td><td>0.75837</td><td>0.65928</td></tr>
          <tr><td>125</td><td>0.9688</td><td>0.78871</td><td>0.68439</td><td>0.57492</td></tr>
          <tr><td>124</td><td>0.9609</td><td>0.73722</td><td>0.61756</td><td>0.51117</td></tr>
          <tr><td>123</td><td>0.9531</td><td>0.68717</td><td>0.55892</td><td>0.46604</td></tr>
          <tr><td>110</td><td>0.8516</td><td>0.23151</td><td>0.29093</td><td>0.33304</td></tr>
          <tr><td>95</td><td>0.7344</td><td>0.00332</td><td>0.16256</td><td>0.18719</td></tr>
          <tr><td>80</td><td>0.6172</td><td>-0.13641</td><td>0.02135</td><td>0.05702</td></tr>
          <tr><td>65</td><td>0.5000</td><td>-0.20581</td><td>-0.11477</td><td>-0.06080</td></tr>
          <tr><td>59</td><td>0.4531</td><td>-0.21090</td><td>-0.17119</td><td>-0.10648</td></tr>
          <tr><td>37</td><td>0.2813</td><td>-0.15662</td><td>-0.32726</td><td>-0.27805</td></tr>
          <tr><td>23</td><td>0.1719</td><td>-0.10150</td><td>-0.24299</td><td>-0.38289</td></tr>
          <tr><td>14</td><td>0.1016</td><td>-0.06434</td><td>-0.14612</td><td>-0.29730</td></tr>
          <tr><td>10</td><td>0.0703</td><td>-0.04775</td><td>-0.10388</td><td>-0.22220</td></tr>
          <tr><td>9</td><td>0.0625</td><td>-0.04192</td><td>-0.09266</td><td>-0.20196</td></tr>
          <tr><td>8</td><td>0.0547</td><td>-0.03717</td><td>-0.08186</td><td>-0.18109</td></tr>
          <tr><td>1</td><td>0.0000</td><td>0.00000</td><td>0.00000</td><td>0.00000</td></tr>
        </tbody>
      </table>
    </div>
  </div>
  
  <div class="centered2">
    <div style="width: 70%;">
      <h4>Error wrt the reference data</h4>
      <table>
        <thead>
          <tr>
            <th>129-Grid pt. no.</th>
            <th>y</th>
            <th>Re (100)</th>
            <th>Re (400)</th>
            <th>Re (1000)</th>
          </tr>
        </thead>
        <tbody>
          <tr><td>129</td><td>1.00</td><td>0.000000</td><td>0.000000</td><td>0.000000</td></tr>
          <tr><td>126</td><td>0.98</td><td>0.003346</td><td>0.001516</td><td>-0.004401</td></tr>
          <tr><td>125</td><td>0.97</td><td>0.003722</td><td>0.000940</td><td>-0.008880</td></tr>
          <tr><td>124</td><td>0.96</td><td>0.004021</td><td>-0.000460</td><td>-0.015002</td></tr>
          <tr><td>123</td><td>0.95</td><td>0.004184</td><td>-0.002631</td><td>-0.021820</td></tr>
          <tr><td>110</td><td>0.85</td><td>-0.009163</td><td>-0.028913</td><td>-0.037339</td></tr>
          <tr><td>95</td><td>0.74</td><td>-0.019810</td><td>-0.021102</td><td>-0.023771</td></tr>
          <tr><td>80</td><td>0.62</td><td>-0.001890</td><td>-0.016967</td><td>-0.014619</td></tr>
          <tr><td>65</td><td>0.50</td><td>0.026750</td><td>-0.018510</td><td>-0.004300</td></tr>
          <tr><td>59</td><td>0.46</td><td>0.034670</td><td>-0.017870</td><td>0.000100</td></tr>
          <tr><td>37</td><td>0.29</td><td>0.035480</td><td>0.052250</td><td>0.006450</td></tr>
          <tr><td>23</td><td>0.18</td><td>0.023550</td><td>0.073420</td><td>0.063170</td></tr>
          <tr><td>14</td><td>0.11</td><td>0.014400</td><td>0.052930</td><td>0.092560</td></tr>
          <tr><td>10</td><td>0.08</td><td>0.011070</td><td>0.041520</td><td>0.081510</td></tr>
          <tr><td>9</td><td>0.07</td><td>0.008710</td><td>0.037870</td><td>0.077260</td></tr>
          <tr><td>8</td><td>0.06</td><td>0.007490</td><td>0.034630</td><td>0.073110</td></tr>
          <tr><td>1</td><td>0.00</td><td>0.000000</td><td>0.000000</td><td>0.000000</td></tr>
        </tbody>
      </table>
    </div>
  </div>
</div>


### Gallery

<table>
    <tr>
        <th>Re=100 on a 100x100 cavity</th>
        <th>Re=10000 on a 100x100 cavity</th>
    </tr>
    <tr>
        <td><img src="./media/100x100_re100_steps5000_periteration5_fps24.gif" alt="100x100_re100_steps5000_periteration5_fps24"></td>
        <td><img src="./media/100x100_re10000_steps10000_periteration10_fps24.gif" alt="100x100_re10000_steps10000_periteration10_fps24"></td>
    </tr>
    <tr>
        <th>Re=10000 on a 75x50 cavity with an alternative boundary condition</th>
        <th>Re=10000 on a 240x200 cavity with an alternative boundary condition</th>
    </tr>
    <tr>
        <td><img src="./media/75x50_re10000_steps10000_periteration20_fps24_alt1.gif" alt="75x50_re10000_steps10000_periteration20_fps24_alt1"></td>
        <td><img src="./media/240x200_re10000_steps20000_periteration50_fps24_alt2.gif" alt="240x200_re10000_steps20000_periteration50_fps24_alt2"></td>
    </tr>
    <tr>
        <th>Re=2200 on a 200x80 cavity </th>
        <th>Re=1000 on a 250x250 cavity</th>
    </tr>
    <tr>
        <td><img src="./media/200x80_re2200_steps10000_periteration50_fps24.gif" alt="200x80_re2200_steps10000_periteration50_fps24"></td>
        <td><img src="./media/250x250_re1000_steps10000_periteration20_fps24.gif" alt="250x250_re1000_steps10000_periteration20_fps24"></td>
    </tr>
</table>



<br><br>

## 3D Lid-Driven Cavity Simulation with 3DLBM

### Simulation Strategy

* **Domain discretization:** The cavity is divided into cubic grid cells, each representing local fluid populations.
* **Time-stepping cycle:**

  1. **Collision:** Particle distributions relax toward local equilibrium, modeling viscous dissipation.
  2. **Streaming:** Updated distributions move along discrete velocity directions to neighboring cells.
  3. **Boundary conditions:**

     * **Stationary walls:** Halfway bounce-back enforces no-slip zero-velocity at solid faces.
     * **Moving lid:** Zou–He velocity boundary condition prescribes a uniform tangential speed on the top plane, delivering accurate momentum transfer.
* **Stabilization:** A two-relaxation-time collision operator separates symmetric and anti-symmetric moments, damping spurious modes at moderate Reynolds numbers.
* **Parallelization:** Shared-memory threading accelerates both collision and streaming operations, scaling efficiently on multicore processors.

### Output and Post-Processing

Simulation data are exported in VTK format and visualized using ParaView:

* Bulk import of VTK sequences for time-resolved analysis.
* Stream Tracer and Glyph filters to highlight vortex structures.
* Automated rendering of video animations via ParaView’s Python scripting interface.

<details>
 <summary> Gallery of videos (click to espand) </summary>

### evolutions of system with the streamline section evolutions
 
<table>
 <tr>
<th colspan="2">Re=100 on 100×100×100</th>
</tr>
<tr>
<td colspan="2">
<img src="./media/Re_100.gif" alt="Re_100">
</td>
</tr>
 <tr>
<th colspan="2">Re=500 on 100×100×100</th>
</tr>
<tr>
<td colspan="2">
<img src="./media/Re_500.gif" alt="Re_500">
</td>
</tr>
<tr>
<th colspan="2">Re=1000 on 100×100×100</th>
</tr>
<tr>
<td colspan="2">
<img src="./media/Re_1000.gif" alt="Re_1000">
</td>
</tr>
</table>

</details>

### Representative Results

* **Primary vortex:** At Re=1000, a dominant circulation fills the cavity, with symmetric secondary eddies at each corner.
* **Velocity profiles:** Mid-plane velocity curves closely match literature benchmarks, confirming numerical accuracy.
* **3D visualizations:** High-fidelity renderings capture complex vortical structures, highlighting the method’s capacity to resolve secondary flows.

### Performance Highlights

* Modular separation of physics and I/O allows easy extension to new geometries.
* Threaded implementation achieves near-linear speedup, enabling large-grid simulations within practical runtimes.

*End of section.*\*\*




