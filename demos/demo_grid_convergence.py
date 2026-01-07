"""
=================================================================================
Demo 6: Grid Convergence Study (Order of Accuracy)
=================================================================================
This script performs a grid dependence test to verify the order of accuracy
of the numerical scheme.

It runs the solver on three successively refined grids:
1. Coarse (N=21)
2. Medium (N=41)
3. Fine   (N=81)

It computes the L2 Error norm against the Analytical Solution for each.
The Order of Accuracy (p) is estimated using the ratio of errors.

Usage:
    python demos/demo_grid_convergence.py
=================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import math

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.numerical.solver_fvm_tvd import BurgersSolverFVM_TVD, FluxLimiter, SolverMethod
from src.analytical.case1_solution import analytical_solution_case1

def run_grid_convergence():
    # 1. Define Grids to test
    # grids = [21, 26, 31, 36, 41, 44, 46, 48, 50, 51]
    # make a grid list from 21 to 41 with step 2
    grids = np.arange(21, 201, 5)
    results = []

    print("\n" + "="*80)
    print("GRID CONVERGENCE STUDY")
    print("="*80)
    print(f"{'Grid (NxN)':<15} {'Dx':<12} {'Iterations':<12} {'L2 Error':<12} {'Order (p)':<10}")
    print("-" * 80)

    # 2. Loop through grids
    prev_error = None
    
    for N in grids:
        # We use VAN_LEER as the representative scheme (stable and 2nd order)
        # You can change this to FluxLimiter.UPWIND to see 1st order behavior
        solver = BurgersSolverFVM_TVD(
            N=N, 
            a=1.0, b=1.0, c=1.5, d=-0.5, 
            cfl=0.5,                    
            limiter=FluxLimiter.UPWIND   
        )
        
        # Solve
        # FIX: Removed 'verbose=False' and added 'print_interval'
        stats = solver.solve_steady_state(
            max_iterations=10000, 
            tolerance=1e-6, 
            method=SolverMethod.EXPLICIT_EULER,
            print_interval=10001 # Set higher than max_iter to silence output
        )
        
        # Compute Error
        errors = solver.compute_error(analytical_solution_case1)
        l2_error = errors['L2']
        
        # Calculate Order of Accuracy (p)
        order_str = "N/A"
        if prev_error is not None:
            # Formula: p = log(Error_coarse / Error_fine) / log(Refinement_ratio)
            # Refinement ratio is roughly 2 (since dx halves)
            # Accurate ratio: ratio = dx_old / dx_new
            
            # Note: We use the previous loop's dx, but since we didn't save it, 
            # we assume ~2.0 refinement.
            p = math.log(prev_error / l2_error) / math.log(2.0)
            order_str = f"{p:.2f}"
            
        print(f"{f'{N}x{N}':<15} {solver.dx:<12.4f} {stats['iterations']:<12} {l2_error:<12.6f} {order_str:<10}")
        
        results.append({
            'N': N,
            'dx': solver.dx,
            'error': l2_error
        })
        
        prev_error = l2_error

    # 3. Plot Convergence (Log-Log Plot)
    dx_values = [r['dx'] for r in results]
    N_grid = [r['N'] for r in results]
    error_values = [r['error'] for r in results]

    plt.figure(figsize=(10, 6))
    # plt.loglog(dx_values, error_values, 'bo-', linewidth=2, markersize=8, label='Numerical Error')
    plt.plot(N_grid, error_values, 'bo-', linewidth=2, markersize=8, label='Numerical Error')
    
    # Add reference slope (Order = 1)
    # y = mx + c  => log(err) = 1 * log(dx) + offset
    if len(dx_values) > 0:
        ref_x = np.array(dx_values)
        ref_N = np.array(N_grid)
        
        # 1st Order Reference Line
        ref_y1 = error_values[0] * (ref_x / ref_x[0])**1
        # plt.loglog(ref_x, ref_y1, 'k--', label='1st Order Ref.')
        plt.plot(ref_N, ref_y1, 'k--', label='1st Order Ref.')
        
        # 2nd Order Reference Line
        ref_y2 = error_values[0] * (ref_x / ref_x[0])**2
        # plt.loglog(ref_x, ref_y2, 'r--', label='2nd Order Ref.')
        plt.plot(ref_N, ref_y2, 'r--', label='2nd Order Ref.')

    # plt.xlabel('Grid Spacing (dx)')
    plt.xlabel('Grid Number (N)')
    plt.ylabel('L2 Error Norm')
    plt.title('Grid Convergence Study (UPWIND)')
    plt.grid(True, which="both", ls="-", alpha=0.4)
    plt.legend()
    
    # Save plot
    save_path = 'plots/convergence/grid_convergence.png'
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=150)
    print(f"\nPlot saved to {save_path}")
    
    # Show plot
    plt.show()

if __name__ == "__main__":
    run_grid_convergence()