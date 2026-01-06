"""
=================================================================================
Demo 3: Comprehensive Comparison of All TVD Flux Limiters
=================================================================================

This script solves the 2D Inviscid Burgers Equation using ALL available TVD flux
limiters and presents a comprehensive comparison with the analytical solution.

Limiters Compared:
    - Minmod, Van Albada, UMIST, Superbee
    - Van Leer, Koren, QUICK, SMART

Output:
    - Grid of contour plots showing all limiters + analytical
    - Error comparison table
    - Summary statistics

Usage:
    python demos/demo_tvd_all_limiters.py

=================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import os
import time

# Add project root to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.numerical.solver_fvm_tvd import BurgersSolverFVM_TVD, FluxLimiter, SolverMethod
from src.analytical.case1_solution import analytical_solution_case1


# All TVD limiters to compare
ALL_LIMITERS = [
    (FluxLimiter.MINMOD, 'Minmod'),
    (FluxLimiter.VAN_ALBADA, 'Van Albada'),
    (FluxLimiter.UMIST, 'UMIST'),
    (FluxLimiter.SUPERBEE, 'Superbee'),
    (FluxLimiter.VAN_LEER, 'Van Leer'),
    (FluxLimiter.KOREN, 'Koren'),
    (FluxLimiter.QUICK, 'QUICK'),
    (FluxLimiter.SMART, 'SMART'),
]


def solve_all_limiters(N: int = 41, max_iterations: int = 10000, 
                       tolerance: float = 1e-6):
    """
    Solve the equation using all available limiters.
    
    Returns
    -------
    results : dict
        Dictionary containing solver and metrics for each limiter
    """
    results = {}
    
    print(f"\nSolving with all {len(ALL_LIMITERS)} limiters on {N}×{N} grid...")
    print("-" * 70)
    
    for limiter_enum, limiter_name in ALL_LIMITERS:
        print(f"  Solving with {limiter_name:12}...", end=" ", flush=True)
        start_time = time.time()
        
        # Create solver
        solver = BurgersSolverFVM_TVD(
            N=N, a=1.0, b=1.0, c=1.5, d=-0.5,
            cfl=0.5, limiter=limiter_enum
        )
        
        # Solve to steady state
        converged, iterations, final_residual = solver.solve_steady_state(
            max_iterations=max_iterations,
            tolerance=tolerance,
            print_interval=0,  # Silent
            method=SolverMethod.EXPLICIT_EULER
        )
        
        elapsed = time.time() - start_time
        
        # Compute errors
        X, Y = np.meshgrid(solver.x, solver.y)
        U_numerical = solver.u[1:-1, 1:-1]
        U_analytical = analytical_solution_case1(X, Y)
        error = np.abs(U_numerical - U_analytical)
        
        L1 = np.mean(error)
        L2 = np.sqrt(np.mean(error**2))
        Linf = np.max(error)
        
        results[limiter_name] = {
            'solver': solver,
            'converged': converged,
            'iterations': iterations,
            'time': elapsed,
            'L1': L1,
            'L2': L2,
            'Linf': Linf,
            'U': U_numerical,
            'X': X,
            'Y': Y,
        }
        
        status = "✓" if converged else "⚠"
        print(f"{status} {iterations:5d} iter, L₂={L2:.4f}, {elapsed:.2f}s")
    
    print("-" * 70)
    
    return results


def print_comparison_table(results: dict):
    """Print a formatted comparison table."""
    print("\n" + "=" * 80)
    print(" " * 20 + "ERROR COMPARISON TABLE")
    print("=" * 80)
    print(f"{'Limiter':<12} {'Status':<10} {'Iterations':>10} {'L₁ Error':>12} "
          f"{'L₂ Error':>12} {'L∞ Error':>12}")
    print("-" * 80)
    
    # Sort by L2 error
    sorted_results = sorted(results.items(), key=lambda x: x[1]['L2'])
    
    for name, data in sorted_results:
        status = "Converged" if data['converged'] else "Max iter"
        print(f"{name:<12} {status:<10} {data['iterations']:>10} {data['L1']:>12.6f} "
              f"{data['L2']:>12.6f} {data['Linf']:>12.6f}")
    
    print("-" * 80)
    
    # Find best limiter
    best = sorted_results[0]
    print(f"\nBest performer (lowest L₂ error): {best[0]} with L₂ = {best[1]['L2']:.6f}")
    print("=" * 80)


def plot_all_comparisons(results: dict, save_path: str = None):
    """
    Create a grid of contour plots comparing all limiters.
    
    Layout: 3x3 grid (8 limiters + 1 analytical)
    """
    fig, axes = plt.subplots(3, 3, figsize=(15, 14))
    axes = axes.flatten()
    
    levels = np.linspace(-0.5, 1.5, 21)
    
    # Plot analytical solution first
    first_result = list(results.values())[0]
    X, Y = first_result['X'], first_result['Y']
    U_analytical = analytical_solution_case1(X, Y)
    
    cf = axes[0].contourf(X, Y, U_analytical, levels=levels, cmap='RdYlBu_r')
    axes[0].contour(X, Y, U_analytical, levels=levels, colors='black', linewidths=0.3, alpha=0.5)
    axes[0].set_title('Analytical Solution', fontweight='bold', fontsize=11)
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    axes[0].set_aspect('equal')
    
    # Plot each limiter
    sorted_results = sorted(results.items(), key=lambda x: x[1]['L2'])
    
    for idx, (name, data) in enumerate(sorted_results):
        ax = axes[idx + 1]
        cf = ax.contourf(X, Y, data['U'], levels=levels, cmap='RdYlBu_r')
        ax.contour(X, Y, data['U'], levels=levels, colors='black', linewidths=0.3, alpha=0.5)
        ax.set_title(f"{name}\nL₂ = {data['L2']:.4f}", fontsize=10)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal')
    
    # Add colorbar
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    fig.colorbar(cf, cax=cbar_ax, label='u(x, y)')
    
    plt.suptitle('2D Inviscid Burgers Equation - All TVD Limiters Comparison\n'
                 '(Sorted by L₂ error, best first)', fontsize=14, fontweight='bold')
    
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"\nComparison plot saved to: {save_path}")
    
    return fig


def plot_error_bars(results: dict, save_path: str = None):
    """Create a bar chart comparing errors across limiters."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    
    # Sort by L2 error
    sorted_items = sorted(results.items(), key=lambda x: x[1]['L2'])
    names = [item[0] for item in sorted_items]
    
    L1_errors = [item[1]['L1'] for item in sorted_items]
    L2_errors = [item[1]['L2'] for item in sorted_items]
    Linf_errors = [item[1]['Linf'] for item in sorted_items]
    
    colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(names)))
    
    # L1 Error
    axes[0].barh(names, L1_errors, color=colors)
    axes[0].set_xlabel('L₁ Error')
    axes[0].set_title('L₁ Error Comparison')
    axes[0].invert_yaxis()
    
    # L2 Error
    axes[1].barh(names, L2_errors, color=colors)
    axes[1].set_xlabel('L₂ Error')
    axes[1].set_title('L₂ Error Comparison')
    axes[1].invert_yaxis()
    
    # L∞ Error
    axes[2].barh(names, Linf_errors, color=colors)
    axes[2].set_xlabel('L∞ Error')
    axes[2].set_title('L∞ Error Comparison')
    axes[2].invert_yaxis()
    
    plt.suptitle('Error Comparison Across All TVD Limiters\n(Sorted from best to worst)',
                 fontsize=12, fontweight='bold')
    plt.tight_layout()
    
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Error bar chart saved to: {save_path}")
    
    return fig


def main():
    """Main function to run comprehensive limiter comparison."""
    
    print("\n" + "=" * 70)
    print(" " * 5 + "2D INVISCID BURGERS EQUATION - ALL TVD LIMITERS COMPARISON")
    print("=" * 70)
    print("\nThis demo solves the equation using ALL available TVD flux limiters")
    print("and creates a comprehensive comparison with the analytical solution.")
    print("\nCase 1 Parameters: a=1, b=1, c=1.5, d=-0.5")
    print("=" * 70)
    
    # Solve with all limiters
    results = solve_all_limiters(N=41, max_iterations=10000)
    
    # Print comparison table
    print_comparison_table(results)
    
    # Generate comparison plot
    print("\nGenerating comparison plots...")
    fig1 = plot_all_comparisons(
        results,
        save_path='plots/tvd_comparison/all_limiters_comparison.png'
    )
    
    # Generate error bar chart
    fig2 = plot_error_bars(
        results,
        save_path='plots/tvd_comparison/all_limiters_errors.png'
    )
    
    print("\n" + "=" * 70)
    print("Comprehensive comparison complete!")
    print("=" * 70 + "\n")
    
    plt.show()


if __name__ == "__main__":
    main()
