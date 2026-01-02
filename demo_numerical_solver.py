"""
Demonstration Script for Numerical Solver
2D Inviscid Burgers Equation using Finite Volume Method

This script demonstrates the numerical solution implementation for Case 1
using the Finite Volume Method with First-Order Upwind Scheme.

Tasks:
- Task 2: Derive discrete equations (FVM)
- Task 3: Implement upwind scheme
- Task 4: Solve for N=21 and N=41 grids
- Task 5: Compute convergence solution
- Task 6: Compare numerical to analytical solution
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from numerical.solver import BurgersSolver2D
from analytical.case1_solution import analytical_solution_case1


def run_case1_analysis():
    """
    Complete analysis for Case 1 with N=21 and N=41 grids.
    """
    print("\n" + "="*70)
    print(" "*15 + "2D INVISCID BURGERS EQUATION")
    print(" "*15 + "Finite Volume Method - Upwind Scheme")
    print("="*70)
    print("\nEquation: du/dt + a*u*du/dx + b*du/dy = 0")
    print("\nCase 1 Parameters:")
    print("  a = 1.0 (nonlinear convection coefficient)")
    print("  b = 1.0 (linear advection coefficient)")
    print("  c = 1.5 (left boundary)")
    print("  d = -0.5 (right boundary)")
    print("\nBoundary Conditions:")
    print("  u(0, y) = 1.5")
    print("  u(1, y) = -0.5")
    print("  u(x, 0) = 1.5 - 2.0*x")
    print("\nDomain: [0, 1] x [0, 1]")
    print("="*70 + "\n")

    # Parameters
    max_iterations = 50000
    tolerance = 1e-6
    cfl = 0.5

    results = {}

    # ===== Solve with N = 21 =====
    print("\n" + "#"*70)
    print(" GRID: N = 21 (21 x 21 = 441 cells)")
    print("#"*70)

    solver_21 = BurgersSolver2D(N=21, a=1.0, b=1.0, c=1.5, d=-0.5, cfl=cfl)
    result_21 = solver_21.solve_steady_state(
        max_iterations=max_iterations,
        tolerance=tolerance,
        print_interval=2000
    )

    # Compute error against analytical solution
    errors_21 = solver_21.compute_error(analytical_solution_case1)
    print("\nError Analysis (N=21):")
    print(f"  L1 error:    {errors_21['L1']:.6e}")
    print(f"  L2 error:    {errors_21['L2']:.6e}")
    print(f"  Linf error:  {errors_21['Linf']:.6e}")
    print(f"  Relative L2: {errors_21['relative_L2']:.6e}")

    results['N21'] = {
        'solver': solver_21,
        'result': result_21,
        'errors': errors_21
    }

    # ===== Solve with N = 41 =====
    print("\n" + "#"*70)
    print(" GRID: N = 41 (41 x 41 = 1681 cells)")
    print("#"*70)

    solver_41 = BurgersSolver2D(N=41, a=1.0, b=1.0, c=1.5, d=-0.5, cfl=cfl)
    result_41 = solver_41.solve_steady_state(
        max_iterations=max_iterations,
        tolerance=tolerance,
        print_interval=2000
    )

    # Compute error against analytical solution
    errors_41 = solver_41.compute_error(analytical_solution_case1)
    print("\nError Analysis (N=41):")
    print(f"  L1 error:    {errors_41['L1']:.6e}")
    print(f"  L2 error:    {errors_41['L2']:.6e}")
    print(f"  Linf error:  {errors_41['Linf']:.6e}")
    print(f"  Relative L2: {errors_41['relative_L2']:.6e}")

    results['N41'] = {
        'solver': solver_41,
        'result': result_41,
        'errors': errors_41
    }

    # ===== Grid Convergence Analysis =====
    print("\n" + "="*70)
    print(" GRID CONVERGENCE ANALYSIS")
    print("="*70)

    print("\nComparison of errors for different grid sizes:")
    print("-"*60)
    print(f"{'Grid':<10} {'L1 Error':<15} {'L2 Error':<15} {'Linf Error':<15}")
    print("-"*60)
    print(f"{'N=21':<10} {errors_21['L1']:<15.6e} {errors_21['L2']:<15.6e} {errors_21['Linf']:<15.6e}")
    print(f"{'N=41':<10} {errors_41['L1']:<15.6e} {errors_41['L2']:<15.6e} {errors_41['Linf']:<15.6e}")
    print("-"*60)

    # Compute convergence order
    h1 = 1.0 / 20  # dx for N=21
    h2 = 1.0 / 40  # dx for N=41

    if errors_21['L2'] > 0 and errors_41['L2'] > 0:
        order_L2 = np.log(errors_21['L2'] / errors_41['L2']) / np.log(h1 / h2)
        print(f"\nObserved convergence order (L2): {order_L2:.2f}")
        print("(Expected: ~1 for first-order upwind scheme)")

    return results


def generate_plots(results: dict, show_plots: bool = True):
    """
    Generate all visualization plots.
    """
    print("\n" + "="*70)
    print(" GENERATING PLOTS")
    print("="*70)

    # Create output directory
    os.makedirs('plots/numerical', exist_ok=True)

    solver_21 = results['N21']['solver']
    solver_41 = results['N41']['solver']

    # 1. Solution plots for N=21
    print("\n1. Generating solution plot for N=21...")
    solver_21.plot_solution(
        save_path='plots/numerical/solution_N21.png',
        show_plot=False,
        title_suffix='Case 1: a=1, b=1, c=1.5, d=-0.5'
    )
    print("   Saved: plots/numerical/solution_N21.png")

    # 2. Solution plots for N=41
    print("\n2. Generating solution plot for N=41...")
    solver_41.plot_solution(
        save_path='plots/numerical/solution_N41.png',
        show_plot=False,
        title_suffix='Case 1: a=1, b=1, c=1.5, d=-0.5'
    )
    print("   Saved: plots/numerical/solution_N41.png")

    # 3. Convergence history
    print("\n3. Generating convergence plots...")
    solver_21.plot_convergence(
        save_path='plots/numerical/convergence_N21.png',
        show_plot=False
    )
    print("   Saved: plots/numerical/convergence_N21.png")

    solver_41.plot_convergence(
        save_path='plots/numerical/convergence_N41.png',
        show_plot=False
    )
    print("   Saved: plots/numerical/convergence_N41.png")

    # 4. Comparison with analytical solution
    print("\n4. Generating comparison plots...")
    solver_21.plot_comparison(
        analytical_solution_case1,
        save_path='plots/numerical/comparison_N21.png',
        show_plot=False
    )
    print("   Saved: plots/numerical/comparison_N21.png")

    solver_41.plot_comparison(
        analytical_solution_case1,
        save_path='plots/numerical/comparison_N41.png',
        show_plot=False
    )
    print("   Saved: plots/numerical/comparison_N41.png")

    # 5. Combined convergence comparison
    print("\n5. Generating combined convergence comparison...")
    fig, ax = plt.subplots(figsize=(10, 6))

    iterations_21 = range(1, len(results['N21']['result']['residual_history']) + 1)
    iterations_41 = range(1, len(results['N41']['result']['residual_history']) + 1)

    ax.semilogy(iterations_21, results['N21']['result']['residual_history'],
                'b-', linewidth=1.5, label='N=21', alpha=0.8)
    ax.semilogy(iterations_41, results['N41']['result']['residual_history'],
                'r-', linewidth=1.5, label='N=41', alpha=0.8)

    ax.set_xlabel('Iteration', fontsize=12)
    ax.set_ylabel('Residual (L2 norm)', fontsize=12)
    ax.set_title('Convergence History Comparison\nCase 1: a=1, b=1, c=1.5, d=-0.5', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11)

    plt.tight_layout()
    plt.savefig('plots/numerical/convergence_comparison.png', dpi=300, bbox_inches='tight')
    print("   Saved: plots/numerical/convergence_comparison.png")
    plt.close()

    # 6. Grid comparison (side by side)
    print("\n6. Generating grid comparison plot...")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # N=21
    cf1 = axes[0].contourf(solver_21.X, solver_21.Y, solver_21.u, levels=20, cmap='RdYlBu_r')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    axes[0].set_title(f'Numerical (N=21)\nL2 Error: {results["N21"]["errors"]["L2"]:.4e}')
    axes[0].set_aspect('equal')
    plt.colorbar(cf1, ax=axes[0])

    # N=41
    cf2 = axes[1].contourf(solver_41.X, solver_41.Y, solver_41.u, levels=20, cmap='RdYlBu_r')
    axes[1].set_xlabel('x')
    axes[1].set_ylabel('y')
    axes[1].set_title(f'Numerical (N=41)\nL2 Error: {results["N41"]["errors"]["L2"]:.4e}')
    axes[1].set_aspect('equal')
    plt.colorbar(cf2, ax=axes[1])

    # Analytical (on fine grid)
    x_fine = np.linspace(0, 1, 100)
    y_fine = np.linspace(0, 1, 100)
    X_fine, Y_fine = np.meshgrid(x_fine, y_fine)
    U_analytical = analytical_solution_case1(X_fine, Y_fine)

    cf3 = axes[2].contourf(X_fine, Y_fine, U_analytical, levels=20, cmap='RdYlBu_r')
    axes[2].set_xlabel('x')
    axes[2].set_ylabel('y')
    axes[2].set_title('Analytical Solution')
    axes[2].set_aspect('equal')
    plt.colorbar(cf3, ax=axes[2])

    plt.tight_layout()
    plt.savefig('plots/numerical/grid_comparison.png', dpi=300, bbox_inches='tight')
    print("   Saved: plots/numerical/grid_comparison.png")
    plt.close()

    print("\n" + "="*70)
    print(" ALL PLOTS GENERATED SUCCESSFULLY")
    print("="*70)

    if show_plots:
        print("\nDisplaying plots...")
        plt.show()


def print_summary_table(results: dict):
    """
    Print a comprehensive summary table.
    """
    print("\n" + "="*70)
    print(" SUMMARY TABLE - CASE 1")
    print("="*70)

    print("\n" + "-"*70)
    print(f"{'Metric':<30} {'N=21':<20} {'N=41':<20}")
    print("-"*70)

    # Grid info
    print(f"{'Grid cells':<30} {'21 x 21 = 441':<20} {'41 x 41 = 1681':<20}")
    print(f"{'Grid spacing (dx)':<30} {1/20:<20.6f} {1/40:<20.6f}")

    # Convergence info
    r21 = results['N21']['result']
    r41 = results['N41']['result']
    print(f"{'Iterations to converge':<30} {r21['iterations']:<20} {r41['iterations']:<20}")
    print(f"{'Final residual':<30} {r21['final_residual']:<20.6e} {r41['final_residual']:<20.6e}")
    print(f"{'Converged':<30} {str(r21['converged']):<20} {str(r41['converged']):<20}")
    print(f"{'Time elapsed (s)':<30} {r21['time_elapsed']:<20.2f} {r41['time_elapsed']:<20.2f}")

    # Error metrics
    e21 = results['N21']['errors']
    e41 = results['N41']['errors']
    print("-"*70)
    print(f"{'L1 Error':<30} {e21['L1']:<20.6e} {e41['L1']:<20.6e}")
    print(f"{'L2 Error':<30} {e21['L2']:<20.6e} {e41['L2']:<20.6e}")
    print(f"{'Linf Error':<30} {e21['Linf']:<20.6e} {e41['Linf']:<20.6e}")
    print(f"{'Relative L2 Error':<30} {e21['relative_L2']:<20.6e} {e41['relative_L2']:<20.6e}")

    print("-"*70)

    # Solution statistics
    s21 = results['N21']['solver']
    s41 = results['N41']['solver']
    print(f"{'Solution min':<30} {s21.u.min():<20.6f} {s41.u.min():<20.6f}")
    print(f"{'Solution max':<30} {s21.u.max():<20.6f} {s41.u.max():<20.6f}")
    print(f"{'Solution mean':<30} {s21.u.mean():<20.6f} {s41.u.mean():<20.6f}")

    print("-"*70 + "\n")


def main():
    """Main function to run the complete numerical analysis."""

    # Create output directories
    os.makedirs('plots/numerical', exist_ok=True)

    # Run analysis
    results = run_case1_analysis()

    # Print summary
    print_summary_table(results)

    # Generate plots
    generate_plots(results, show_plots=False)

    print("\n" + "="*70)
    print(" ANALYSIS COMPLETED SUCCESSFULLY!")
    print("="*70)
    print("\nOutput files saved in plots/numerical/:")
    print("  - solution_N21.png")
    print("  - solution_N41.png")
    print("  - convergence_N21.png")
    print("  - convergence_N41.png")
    print("  - comparison_N21.png")
    print("  - comparison_N41.png")
    print("  - convergence_comparison.png")
    print("  - grid_comparison.png")
    print("\nRun with show_plots=True in generate_plots() to display interactively.")
    print("="*70 + "\n")

    return results


if __name__ == "__main__":
    results = main()
