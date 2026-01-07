# 2D Inviscid Burgers Equation CFD Project

This project implements a numerical solver for the **2D Inviscid Burgers Equation** using the **Finite Volume Method (FVM)** with **TVD (Total Variation Diminishing)** schemes. It includes multiple flux limiters and an exact analytical solution for validation.

The governing equation solved is:
```math
\frac{\partial u}{\partial t} + a u \frac{\partial u}{\partial x} + b \frac{\partial u}{\partial y} = 0
```

Where:
- $a$: Coefficient for nonlinear convection term in x-direction (Case 1: $a=1.0$)
- $b$: Coefficient for linear advection term in y-direction (Case 1: $b=1.0$)

## Features

- **Finite Volume Method (FVM)** discretization
- **TVD Schemes** with multiple flux limiters:
  - Minmod, Van Albada, UMIST, Superbee
  - Van Leer, Koren, QUICK, SMART
- **First-Order Upwind Scheme** for baseline comparison
- **Explicit Euler** time integration with CFL-based adaptive time stepping
- **Verification** against an exact analytical solution (Case 1)

## Project Structure

```
├── demos/
│   ├── demo_analytical_solution.py              # Visualize analytical solution
│   ├── demo_tvd_single_limiter.py               # Solve with user-selected limiter
│   ├── demo_tvd_all_limiters.py                 # Compare all TVD limiters
│   ├── demo_quick_upwind_comparison.py          # Run QUICK vs Upwind comparison
│   ├── demo_grid_convergence.py                 # Run grid convergence study with Upwind
│   ├── demo_grid_visualization.py               # Visualize computational grids
│   ├── demo_flowcharts.py                       # Generate academic flowcharts
│   ├── demo_tvd_limiter_comparisons.py          # Detailed limiter comparisons
│   ├── demo_analytical_numerical_comparisons.py # Numerical vs analytical comparisons
│   ├── demo_all_limiters_comprehensive_plot.py  # Comprehensive limiter analysis
│   └── run_cases.py                             # Run predefined test cases
├── src/
│   ├── analytical/                              # Analytical solution implementation
│   │   └── case1_solution.py
│   └── numerical/                               # Numerical solver implementations
│       ├── solver.py                            # Original FVM solver
│       └── solver_fvm_tvd.py                    # TVD solver with flux limiters
├── scripts/                                     # Utility scripts
│   ├── analysis/                                # Analysis tools
│   ├── demos/                                   # Additional demo utilities
│   └── tests/                                   # Test scripts
├── docs/                                        # Documentation
│   ├── FINAL_REPORT.md                          # Comprehensive project report
│   ├── DERIVATIONS.md                           # Mathematical derivations
│   ├── HOW_TO_RUN.md                            # Running instructions
│   ├── PLOTS_GUIDE.md                           # Plotting guidelines
│   ├── PROJECT_SUMMARY.md                       # Project summary
│   └── GIT_WORKFLOW.md                          # Git workflow guide
├── References/                                  # Reference materials
├── plots/                                       # Generated output plots
├── cases.yaml                                   # Test case configurations
├── regenerate_upwind_corrected.py               # Upwind regeneration script
├── verify_analytical_grid_independence.py       # Grid independence verification
├── test_tdma_solver.py                          # TDMA solver tests
├── TDMA_IMPLEMENTATION_GUIDE.md                 # TDMA implementation guide
├── .gitignore
├── README.md
└── requirements.txt
```

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/FlameEnjoyer/Major-Task-CFD-1-Burger-Equation.git
   cd Major-Task-CFD-1-Burger-Equation
   ```

2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Demo 1: Analytical Solution

Visualize the exact analytical solution to understand the shock wave structure:

```bash
python demos/demo_analytical_solution.py
```

**Output:**
- Contour plot of the analytical solution
- 3D surface plot
- Solution profiles at constant x and y values

### Demo 2: Single TVD Limiter (Interactive)

Solve using a user-selected TVD flux limiter and compare with the analytical solution:

```bash
python demos/demo_tvd_single_limiter.py
```

**What this does:**
- Presents a menu of available flux limiters
- Solves the equation using your chosen limiter
- Displays side-by-side comparison (Numerical vs Analytical)
- Reports error metrics (L₁, L₂, L∞)

### Demo 3: All TVD Limiters Comparison

Comprehensive comparison of all available TVD flux limiters:

```bash
python demos/demo_tvd_all_limiters.py
```

**What this does:**
- Solves with all 8 TVD limiters
- Generates a grid of contour plots for visual comparison
- Creates error bar chart ranking all limiters
- Prints a summary table with error metrics

### Demo 4: Upwind vs QUICK Comparison

Compare first-order Upwind scheme with QUICK method:

```bash
python demos/demo_quick_upwind_comparison.py
```

**What this does:**
- Solves using both Upwind and QUICK methods
- Side-by-side contour plots
- Error metrics comparison
- Demonstrates accuracy differences between 1st and higher-order schemes

### Demo 5: Grid Convergence Study

Perform grid convergence analysis using Upwind scheme:

```bash
python demos/demo_grid_convergence.py
```

**What this does:**
- Tests multiple grid resolutions (N = 21 to 201)
- Computes L2 error norms for each grid
- Generates convergence plot with reference slopes
- Calculates numerical order of accuracy

## Available Flux Limiters

| Limiter | Description |
|---------|-------------|
| Minmod | Most diffusive, very stable |
| Van Albada | Smooth, balanced accuracy/stability |
| UMIST | Sharp shock capturing |
| Superbee | Most compressive |
| Van Leer | Good balance |
| Koren | Third-order accurate |
| QUICK | Quadratic upstream interpolation |
| SMART | High resolution |

## Dependencies

- **NumPy**: Matrix operations and grid handling
- **Matplotlib**: Visualization and plotting
- **SciPy**: Scientific computing utilities
