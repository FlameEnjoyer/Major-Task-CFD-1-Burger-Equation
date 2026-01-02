# Quick Start Guide

## 2D Inviscid Burgers Equation Solver

This repository contains a Finite Volume Method (FVM) solver for the 2D Inviscid Burgers Equation.

### Prerequisites

- Python 3.8 or higher
- Git

### Installation

1. **Clone the repository:**
   ```bash
   git clone <repository_url>
   cd burgers-equation-cfd
   ```

2. **Install dependencies:**
   It is recommended to use a virtual environment.
   ```bash
   pip install -r requirements.txt
   ```

### Running the Code

#### 1. Numerical Solution (Main Task)
To run the numerical solver, generate grids, and perform convergence analysis:
```bash
python demo_numerical_solver.py
```
This will:
- Generate grid visualizations for N=21 and N=41.
- Solve the equation until steady state.
- Generate solution plots and convergence graphs.
- Save all figures to `plots/numerical/`.

#### 2. Analytical Solution
To generate the exact analytical solution and visualize shock waves:
```bash
python demo_analytical_case1.py
```
This will:
- Verify boundary conditions.
- Generate contour plots and solution slices.
- Save all figures to `plots/analytical/`.

### Output

All generated plots are saved in the `plots/` directory (created automatically).
