"""
Multi-Case Study Runner for 2D Inviscid Burgers Equation
(Prompt 4: Convergence and Multi-Case Study)

This script executes the solver for 4 different parameter cases:
1. Case 1 (Reference): a=1, b=1, c=1.5, d=-0.5
2. Case 2: a=1.5, b=0.5, c=1, d=-0.5
3. Case 3: a=1.5, b=2, c=1, d=-0.5
4. Case 4: a=1, b=1, c=1, d=-0.6

It generates:
- Convergence plots (Residual vs Iterations)
- Side-by-side plots of Numerical Result vs Grid Density
- Error Map for Case 1 (Numerical - Analytical)
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Add src to path to import solver
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from numerical.solver import BurgersSolver2D
from analytical.case1_solution import analytical_solution_case1

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def run_multicase_study():
    # Output directory
    plot_dir = 'plots/multicase'
    ensure_dir(plot_dir)

    print("\n" + "="*80)
    print(" MULTI-CASE STUDY RUNNER")
    print("="*80)

    # Define cases
    cases = [
        {
            'name': 'Case 1',
            'params': {'a': 1.0, 'b': 1.0, 'c': 1.5, 'd': -0.5},
            'desc': 'Reference (a=1, b=1, c=1.5, d=-0.5)',
            'analytical': analytical_solution_case1
        },
        {
            'name': 'Case 2',
            'params': {'a': 1.5, 'b': 0.5, 'c': 1.0, 'd': -0.5},
            'desc': 'a=1.5, b=0.5, c=1, d=-0.5',
            'analytical': None
        },
        {
            'name': 'Case 3',
            'params': {'a': 1.5, 'b': 2.0, 'c': 1.0, 'd': -0.5},
            'desc': 'a=1.5, b=2.0, c=1, d=-0.5',
            'analytical': None
        },
        {
            'name': 'Case 4',
            'params': {'a': 1.0, 'b': 1.0, 'c': 1.0, 'd': -0.6},
            'desc': 'a=1, b=1, c=1, d=-0.6',
            'analytical': None
        }
    ]

    # Solver settings (Standardized)
    N = 41 # Use fine grid for better visualization
    max_iter = 10000
    tol = 1e-6
    cfl = 0.5

    for case in cases:
        name = case['name']
        print(f"\nRunning {name}: {case['desc']}")
        print("-" * 40)
        
        # Initialize solver
        p = case['params']
        solver = BurgersSolver2D(N=N, a=p['a'], b=p['b'], c=p['c'], d=p['d'], cfl=cfl)
        
        # Solve
        result = solver.solve_steady_state(
            max_iterations=max_iter, 
            tolerance=tol, 
            print_interval=2000
        )
        
        # 1. Plot Convergence (Residual vs Iterations)
        print(f"  Generating convergence plot...")
        solver.plot_convergence(
            save_path=f'{plot_dir}/{name.replace(" ", "_")}_convergence.png',
            show_plot=False
        )

        # 2. Side-by-Side: Numerical Result + Grid Density
        print(f"  Generating Result + Grid density comparison...")
        fig = plt.figure(figsize=(12, 5))
        
        # Left: Numerical Result
        ax1 = fig.add_subplot(121)
        cf = ax1.contourf(solver.X, solver.Y, solver.u, levels=20, cmap='RdYlBu_r')
        plt.colorbar(cf, ax=ax1, label='u')
        ax1.set_title(f'{name}: Numerical Solution')
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_aspect('equal')

        # Right: Grid Density (Visualized as mesh lines)
        ax2 = fig.add_subplot(122)
        # Plot mesh lines
        for i in range(solver.N):
            ax2.plot(solver.X[:, i], solver.Y[:, i], 'k-', linewidth=0.5, alpha=0.5)
        for j in range(solver.N):
            ax2.plot(solver.X[j, :], solver.Y[j, :], 'k-', linewidth=0.5, alpha=0.5)
            
        ax2.set_title(f'{name}: Grid Density (N={N})')
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        ax2.set_xlim(0, 1)
        ax2.set_ylim(0, 1)
        ax2.set_aspect('equal')
        
        plt.tight_layout()
        plt.savefig(f'{plot_dir}/{name.replace(" ", "_")}_result_grid.png', dpi=300)
        plt.close()

        # 3. Special Task for Case 1: Error Map
        if case['analytical']:
            print(f"  Generating Error Map for {name}...")
            # Compute analytical on numerical grid
            u_exact = case['analytical'](solver.X, solver.Y)
            error = np.abs(solver.u - u_exact)
            
            fig_err = plt.figure(figsize=(6, 5))
            ax_err = fig_err.add_subplot(111)
            cf_err = ax_err.contourf(solver.X, solver.Y, error, levels=20, cmap='inferno')
            plt.colorbar(cf_err, ax=ax_err, label='Absolute Error')
            ax_err.set_title(f'{name}: Error Map |u_num - u_exact|')
            ax_err.set_xlabel('x')
            ax_err.set_ylabel('y')
            ax_err.set_aspect('equal')
            
            plt.tight_layout()
            plt.savefig(f'{plot_dir}/{name.replace(" ", "_")}_error_map.png', dpi=300)
            plt.close()

    print("\n" + "="*80)
    print(f" ALL CASES COMPLETED. Results saved in '{plot_dir}/'")
    print("="*80 + "\n")

if __name__ == "__main__":
    run_multicase_study()
