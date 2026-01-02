# 2D Inviscid Burgers Equation CFD Project

This project implements a numerical solver for the **2D Inviscid Burgers Equation** using the **Finite Volume Method (FVM)** with a first-order **Upwind Scheme**. It also includes an exact analytical solution for validation purposes.

The governing equation solved is:
```math
\frac{\partial u}{\partial t} + a u \frac{\partial u}{\partial x} + b \frac{\partial u}{\partial y} = 0
```

Where:
- $a$: Coefficient for nonlinear convection term in x-direction (Case 1: $a=1.0$)
- $b$: Coefficient for linear advection term in y-direction (Case 1: $b=1.0$)

## Features

- **Finite Volume Method (FVM)** discretization.
- **First-Order Upwind Scheme** for convective fluxes to handle shock waves.
- **Explicit Euler** time integration with CFL-based adaptive time stepping.
- **Grid Generation Visualization** for N=21 and N=41 meshes.
- **Verification** against an exact analytical solution (Case 1).
- **Grid Convergence Analysis** to demonstrate solution accuracy.

## Project Structure

```
├── src/
│   ├── analytical/          # Analytical solution implementation
│   └── numerical/           # Numerical FVM solver implementation
├── demo_numerical_solver.py # Main script to run numerical analysis
├── demo_analytical_case1.py # Script to visualize analytical solution
├── plots/                   # Generated output plots
└── requirements.txt         # Python dependencies
```

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/FlameEnjoyer/Major-Task-CFD-1-Burger-Equation.git
   cd Major-Task-CFD-1-Burger-Equation
   ```

2. **Install dependencies:**
   It is recommended to use a virtual environment.
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### 1. Run Numerical Analysis (Main Task)

To run the complete numerical simulation, generate grids, and perform error analysis:

```bash
python demo_numerical_solver.py
```

**What this does:**
- Generates and plots the computational grid for **N=21** and **N=41**.
- Solves the inviscid Burgers equation until **steady state** is reached.
- Computes error metrics (L1, L2, L∞) against the analytical solution.
- Generates the following visualizations in `plots/numerical/`:
  - `grid_N21.png`, `grid_N41.png`: Physical/Computational domain grids.
  - `solution_N*.png`: Filled contour and 3D surface plots of the solution.
  - `convergence_N*.png`: Residual history showing convergence to steady state.
  - `comparison_N*.png`: Side-by-side comparison with the analytical solution.
  - `grid_comparison.png`: Visual comparison of results across grid sizes.

### 2. Run Analytical Solution

To visualize the exact analytical solution and understand the shock wave structures:

```bash
python demo_analytical_case1.py
```

**What this does:**
- Verifies boundary conditions.
- Generates contour plots showing the shock discontinuities.
- Plots solution slices at constant x and y values.
- Saves figures to `plots/analytical/`.

## Dependencies

- **NumPy**: Matrix operations and grid handling.
- **Matplotlib**: Visualization and plotting.
- **SciPy**: Scientific computing utilities.
