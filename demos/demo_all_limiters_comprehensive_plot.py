"""
Generate comprehensive comparison plot showing all TVD limiters + analytical solution.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.numerical.solver_fvm_tvd import BurgersSolverFVM_TVD, FluxLimiter
from src.analytical.case1_solution import analytical_solution_case1


def create_comprehensive_comparison(N=41, save_path=None):
    """Create comprehensive 2D contour comparison plot for all TVD limiters."""
    print("Generating comprehensive TVD limiter comparison...")

    # Define all limiters to compare
    limiters = [
        (None, "Analytical"),  # None means analytical
        (FluxLimiter.UPWIND, "Upwind"),
        (FluxLimiter.MINMOD, "Min-Mod"),
        (FluxLimiter.VAN_LEER, "Van Leer"),
        (FluxLimiter.VAN_ALBADA, "Albada"),
        (FluxLimiter.SUPERBEE, "Superbee"),
        (FluxLimiter.QUICK, "QUICK"),
        (FluxLimiter.UMIST, "UMIST"),
    ]

    # Create figure with 2x4 subplots
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()

    # Solve all limiters and store solutions
    solutions = []
    X, Y = None, None

    for i, (limiter, name) in enumerate(limiters):
        print(f"\n{i+1}. Processing {name}...")

        if limiter is None:
            # Analytical solution
            if X is None:
                # Create grid for analytical solution
                x = np.linspace(0, 1, N)
                y = np.linspace(0, 1, N)
                X, Y = np.meshgrid(x, y)

            U = analytical_solution_case1(X, Y)
            solutions.append((X, Y, U, name, None))
        else:
            # Numerical solution
            solver = BurgersSolverFVM_TVD(
                N=N, a=1.0, b=1.0, c=1.5, d=-0.5, cfl=0.5, limiter=limiter
            )

            result = solver.solve_steady_state(
                max_iterations=20000,
                tolerance=1e-6,
                print_interval=10000
            )

            errors = solver.compute_error(analytical_solution_case1)

            X, Y = solver.X, solver.Y
            U = solver.u

            print(f"   Converged: {result['converged']}, Iterations: {result['iterations']}")
            print(f"   L2 error: {errors['L2']:.6e}")

            solutions.append((X, Y, U, name, errors['L2']))

    # Find global min/max for consistent colorbar
    all_values = [sol[2] for sol in solutions]
    vmin = min([v.min() for v in all_values])
    vmax = max([v.max() for v in all_values])
    levels = np.linspace(vmin, vmax, 25)

    # Plot all solutions
    for i, (X_sol, Y_sol, U_sol, name, error) in enumerate(solutions):
        ax = axes[i]

        contourf = ax.contourf(X_sol, Y_sol, U_sol, levels=levels,
                               cmap='RdYlBu_r', vmin=vmin, vmax=vmax)

        ax.set_xlabel('$x$', fontsize=11)
        ax.set_ylabel('$y$', fontsize=11)

        if error is not None:
            title = f'{name}\n(L2 error: {error:.4f})'
        else:
            title = f'{name}\n(Exact Solution)'

        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_aspect('equal')
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])

    # Add shared colorbar
    fig.subplots_adjust(right=0.92, wspace=0.3, hspace=0.35)
    cbar_ax = fig.add_axes([0.94, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(contourf, cax=cbar_ax)
    cbar.set_label('$u$', fontsize=13, rotation=0, labelpad=15)

    # Overall title
    fig.suptitle(f'Comprehensive Comparison: Analytical and TVD FVM Solutions ({N}x{N} Grid)',
                 fontsize=16, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 0.92, 0.96])

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"\n[SAVED] {save_path}")

    plt.close()

    # Print summary
    print("\n" + "="*80)
    print(" SUMMARY OF ALL LIMITERS (N=" + str(N) + ")")
    print("="*80)
    print(f"{'Limiter':<20} {'L2 Error':<15}")
    print("-"*80)
    for X_sol, Y_sol, U_sol, name, error in solutions:
        if error is not None:
            print(f"{name:<20} {error:<15.6e}")
        else:
            print(f"{name:<20} {'(Analytical)':<15}")
    print("="*80 + "\n")


def main():
    """Main function to generate comprehensive comparison plot."""
    output_dir = 'plots/tvd_comparisons'
    os.makedirs(output_dir, exist_ok=True)

    save_path = os.path.join(output_dir, 'all_limiters_comprehensive_N41.png')

    create_comprehensive_comparison(N=41, save_path=save_path)

    print("\n" + "="*80)
    print(" COMPREHENSIVE COMPARISON PLOT GENERATED SUCCESSFULLY!")
    print("="*80)
    print(f"\nOutput file: {save_path}")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()
