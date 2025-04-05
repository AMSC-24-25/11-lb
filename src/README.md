# Lattice Boltzmann Method (D2Q9) Implementation

This repository contains a C++ implementation of the Lattice Boltzmann Method (LBM) using the D2Q9 model for fluid dynamics simulations. The implementation focuses on simulating lid-driven cavity flow.

## Mathematical Foundation

### D2Q9 Model Description

The D2Q9 model uses 9 discrete velocities in a 2D lattice. The velocities are defined as:

```
e₀ = (1,1)    e₁ = (1,0)    e₂ = (1,-1)
e₃ = (0,1)    e₄ = (0,0)    e₅ = (0,-1)
e₆ = (-1,1)   e₇ = (-1,0)   e₈ = (-1,-1)
```

With corresponding weights:

```
w₀,w₂,w₆,w₈ = 1/36
w₁,w₃,w₅,w₇ = 1/9
w₄ = 4/9
```

### Core Equations

1. **Equilibrium Distribution Function**

The equilibrium distribution function is calculated as:

$f_i^{eq} = w_i \rho \left[1 + 3(\mathbf{e}_i \cdot \mathbf{u}) + \frac{9}{2}(\mathbf{e}_i \cdot \mathbf{u})^2 - \frac{3}{2}\mathbf{u}^2\right]$

where:
- $w_i$ are the weights
- $\rho$ is the density
- $\mathbf{e}_i$ are the discrete velocities
- $\mathbf{u}$ is the macroscopic velocity

2. **Collision Step**

The BGK collision operator is implemented as:

$f_i(\mathbf{x} + \mathbf{e}_i\Delta t, t + \Delta t) = f_i(\mathbf{x}, t) + \frac{1}{\tau_f}[f_i^{eq}(\mathbf{x}, t) - f_i(\mathbf{x}, t)]$

where:
- $\tau_f$ is the relaxation time
- $\Delta t$ is the time step

3. **Macroscopic Quantities**

The macroscopic density and velocity are computed as:

$\rho = \sum_i f_i$

$\mathbf{u} = \frac{1}{\rho}\sum_i \mathbf{e}_i f_i$

### Key Parameters

- Reynolds number: $Re = \frac{u_{lid}L_x}{\nu}$
- Kinematic viscosity: $\nu = \frac{u_{lid}L_x}{Re}$
- Relaxation time: $\tau_f = 3\nu + 0.5$

## Implementation Details

### Class Structure

The `LBM` class implements the D2Q9 model with the following key components:

1. **State Variables**
- `rho`, `rho2`: Current and updated density fields
- `u`, `u2`: Current and updated velocity fields
- `f`, `F`: Current and updated distribution functions

2. **Physical Parameters**
- `dx`, `dy`: Spatial steps
- `dt`: Time step
- `c`: Speed of sound
- `nu`: Kinematic viscosity
- `tau_f`: Relaxation time
- `Re`: Reynolds number
- `u_lid`: Lid velocity

### Key Methods

1. **Evolution**
Allows the user to perform one or multiple iterations
```cpp
// single iteration
lbm.evolution();

// multiple iterations
lbm.evolution(iterations)
```

2. **Collision and Streaming**
The `compute()` method implements:
- Distribution function updates using BGK collision operator
- Macroscopic quantity calculations
- Streaming step through index manipulation
And all these operations are performed within a single (nested) parallelized for-loop

3. **Boundary Conditions**
The code implements:
- No-slip conditions on walls
- Moving lid at the top boundary
- Density extrapolation at boundaries

### Parallelization

The implementation uses OpenMP for parallel computation:
- Collision and streaming steps are parallelized using `#pragma omp parallel for`
- Boundary conditions are handled with parallel regions
- Memory access is optimized to reduce cache misses

## Usage

Initialize the LBM solver with:

```cpp
LBM solver(nx, ny, u_lid, Re);
```

Run the simulation:

```cpp
// Single step
solver.evolution();

// Multiple steps
solver.evolution(num_iterations);
```

## Notes

- The code uses lattice units with $\Delta x = \Delta y = \Delta t = 1$
- Memory management is handled manually with dynamic allocation
- The implementation is optimized for cache efficiency through careful array ordering
- There are two set of state variables arrays, so that the whole computation can be performed within a single parallelizable loop 