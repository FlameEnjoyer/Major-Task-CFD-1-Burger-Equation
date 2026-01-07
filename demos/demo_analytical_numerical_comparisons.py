"""
Generate Side-by-Side Comparisons of Analytical and Numerical Solutions
For 2D Inviscid Burgers Equation

This script creates publication-quality comparison plots showing:
- 2D contour plots (Analytical vs Numerical) side-by-side
- 3D surface plots (Analytical vs Numerical) side-by-side
- For both N=21 and N=41 grid sizes

Following the format from the reference figures.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.numerical.solver import BurgersSolver2D
from src.analytical.case1_solution import analytical_solution_case1


def plot_2d_comparison(solver, N, save_path=None):
    """
    Create side-by-side 2D contour comparison plot.

    Left: Analytical Solution
    Right: Numerical Solution (Upwind FVM)

    Parameters
    ----------
    solver : BurgersSolver2D
        Solved numerical solver instance
    N : int
        Grid size
    save_path : str, optional
        Path to save the figure
    """
    # Get numerical solution
    X_num = solver.X
    Y_num = solver.Y
    U_num = solver.u

    # Compute analytical solution on same grid
    U_analytical = analytical_solution_case1(X_num, Y_num)

    # Create figure with 2 subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # Common colorbar limits
    vmin = min(U_analytical.min(), U_num.min())
    vmax = max(U_analytical.max(), U_num.max())
    levels = np.linspace(vmin, vmax, 25)

    # Left subplot: Analytical Solution
    contourf1 = ax1.contourf(X_num, Y_num, U_analytical, levels=levels,
                             cmap='RdYlBu_r', vmin=vmin, vmax=vmax)
    ax1.set_xlabel('$x$', fontsize=13)
    ax1.set_ylabel('$y$', fontsize=13)
    ax1.set_title('Analytical Solution', fontsize=13, fontweight='bold')
    ax1.set_aspect('equal')
    ax1.set_xlim([0, 1])
    ax1.set_ylim([0, 1])

    # Right subplot: Numerical Solution
    contourf2 = ax2.contourf(X_num, Y_num, U_num, levels=levels,
                             cmap='RdYlBu_r', vmin=vmin, vmax=vmax)
    ax2.set_xlabel('$x$', fontsize=13)
    ax2.set_ylabel('$y$', fontsize=13)
    ax2.set_title('Numerical Solution', fontsize=13, fontweight='bold')
    ax2.set_aspect('equal')
    ax2.set_xlim([0, 1])
    ax2.set_ylim([0, 1])

    # Add shared colorbar
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(contourf2, cax=cbar_ax)
    cbar.set_label('$u$', fontsize=13, rotation=0, labelpad=15)

    # Overall title
    fig.suptitle(f'Comparison of Analytical and Upwind FVM Solutions in 2D Contour with {N} Grid',
                 fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 0.88, 0.96])

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"[SAVED] {save_path}")

    plt.close()

    return fig


def plot_3d_comparison(solver, N, save_path=None):
    """
    Create side-by-side 3D surface comparison plot.

    Left: Analytical Solution
    Right: Numerical Solution (Upwind FVM)

    Parameters
    ----------
    solver : BurgersSolver2D
        Solved numerical solver instance
    N : int
        Grid size
    save_path : str, optional
        Path to save the figure
    """
    # Get numerical solution
    X_num = solver.X
    Y_num = solver.Y
    U_num = solver.u

    # Compute analytical solution on same grid
    U_analytical = analytical_solution_case1(X_num, Y_num)

    # Create figure with 2 subplots side by side
    fig = plt.figure(figsize=(16, 7))

    # Common colorbar limits
    vmin = min(U_analytical.min(), U_num.min())
    vmax = max(U_analytical.max(), U_num.max())

    # Left subplot: Analytical Solution (3D)
    ax1 = fig.add_subplot(121, projection='3d')
    surf1 = ax1.plot_surface(X_num, Y_num, U_analytical, cmap='RdYlBu_r',
                             edgecolor='none', alpha=0.95, vmin=vmin, vmax=vmax,
                             linewidth=0, antialiased=True, shade=True)

    ax1.set_xlabel('$x$', fontsize=12, labelpad=8)
    ax1.set_ylabel('$y$', fontsize=12, labelpad=8)
    ax1.set_zlabel('$u$', fontsize=12, labelpad=8)
    ax1.set_title('Analytical Solution', fontsize=13, fontweight='bold', pad=15)
    ax1.view_init(elev=25, azim=45)
    ax1.set_xlim([0, 1])
    ax1.set_ylim([0, 1])
    ax1.set_zlim([vmin - 0.1, vmax + 0.1])

    # Right subplot: Numerical Solution (3D)
    ax2 = fig.add_subplot(122, projection='3d')
    surf2 = ax2.plot_surface(X_num, Y_num, U_num, cmap='RdYlBu_r',
                             edgecolor='none', alpha=0.95, vmin=vmin, vmax=vmax,
                             linewidth=0, antialiased=True, shade=True)

    ax2.set_xlabel('$x$', fontsize=12, labelpad=8)
    ax2.set_ylabel('$y$', fontsize=12, labelpad=8)
    ax2.set_zlabel('$u$', fontsize=12, labelpad=8)
    ax2.set_title('Numerical Solution', fontsize=13, fontweight='bold', pad=15)
    ax2.view_init(elev=25, azim=45)
    ax2.set_xlim([0, 1])
    ax2.set_ylim([0, 1])
    ax2.set_zlim([vmin - 0.1, vmax + 0.1])

    # Add shared colorbar
    fig.subplots_adjust(right=0.88, wspace=0.15)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(surf2, cax=cbar_ax, shrink=0.8, aspect=20)
    cbar.set_label('$u$', fontsize=13, rotation=0, labelpad=15)

    # Overall title
    fig.suptitle(f'Comparison of Analytical and Upwind FVM Solutions in 3D Surface with {N} Grid',
                 fontsize=14, fontweight='bold', y=0.95)

    plt.tight_layout(rect=[0, 0, 0.88, 0.93])

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"[SAVED] {save_path}")

    plt.close()

    return fig


def run_comparison_for_grid_size(N, output_dir='plots/comparisons'):
    """
    Run complete comparison for a given grid size.

    Parameters
    ----------
    N : int
        Grid size (N x N)
    output_dir : str
        Directory to save output plots
    """
    print("\n" + "="*70)
    print(f" RUNNING COMPARISON FOR N={N}")
    print("="*70)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Initialize and solve numerically using explicit FVM method
    print(f"\n1. Solving numerically with Explicit Upwind FVM (N={N})...")
    solver = BurgersSolver2D(N=N, a=1.0, b=1.0, c=1.5, d=-0.5, cfl=0.5)

    result = solver.solve_steady_state(
        max_iterations=20000,
        tolerance=1e-6,
        print_interval=5000,
        use_implicit=False,  # Use explicit method
        use_fvm=True,        # Use proper FVM
        use_vectorized=True  # Use vectorized computation
    )

    # Compute error
    errors = solver.compute_error(analytical_solution_case1)

    print(f"\nNumerical Solution Statistics:")
    print(f"  Converged: {result['converged']}")
    print(f"  Iterations: {result['iterations']}")
    print(f"  Final residual: {result['final_residual']:.6e}")
    print(f"  Time elapsed: {result['time_elapsed']:.2f} s")

    print(f"\nError vs Analytical Solution:")
    print(f"  L1 error:  {errors['L1']:.6e}")
    print(f"  L2 error:  {errors['L2']:.6e}")
    print(f"  Linf error: {errors['Linf']:.6e}")

    # Generate 2D comparison plot
    print(f"\n2. Generating 2D contour comparison plot...")
    save_path_2d = f'{output_dir}/comparison_2D_N{N}.png'
    plot_2d_comparison(solver, N, save_path=save_path_2d)

    # Generate 3D comparison plot
    print(f"\n3. Generating 3D surface comparison plot...")
    save_path_3d = f'{output_dir}/comparison_3D_N{N}.png'
    plot_3d_comparison(solver, N, save_path=save_path_3d)

    print(f"\n" + "="*70)
    print(f" COMPLETED COMPARISON FOR N={N}")
    print("="*70)

    return solver, errors, result


def main():
    """
    Main function to generate all comparisons.
    """
    print("\n" + "#"*70)
    print("#" + " "*68 + "#")
    print("#" + " "*10 + "ANALYTICAL vs NUMERICAL SOLUTION COMPARISON" + " "*15 + "#")
    print("#" + " "*15 + "2D Inviscid Burgers Equation" + " "*26 + "#")
    print("#" + " "*68 + "#")
    print("#"*70)

    output_dir = 'plots/comparisons'

    # Run for N=21
    print("\n\n")
    solver_21, errors_21, result_21 = run_comparison_for_grid_size(N=21, output_dir=output_dir)

    # Run for N=41
    print("\n\n")
    solver_41, errors_41, result_41 = run_comparison_for_grid_size(N=41, output_dir=output_dir)

    # Summary
    print("\n\n" + "="*70)
    print(" FINAL SUMMARY")
    print("="*70)

    print("\nComparison for N=21:")
    print(f"  Converged: {result_21['converged']}, Iterations: {result_21['iterations']}")
    print(f"  L2 Error: {errors_21['L2']:.6e}")
    print(f"  Time: {result_21['time_elapsed']:.2f} s")

    print("\nComparison for N=41:")
    print(f"  Converged: {result_41['converged']}, Iterations: {result_41['iterations']}")
    print(f"  L2 Error: {errors_41['L2']:.6e}")
    print(f"  Time: {result_41['time_elapsed']:.2f} s")

    print("\n" + "="*70)
    print(" ALL PLOTS GENERATED SUCCESSFULLY!")
    print("="*70)
    print(f"\nOutput files saved in {output_dir}/:")
    print("  - comparison_2D_N21.png  (Figure 4.2 style)")
    print("  - comparison_3D_N21.png  (Figure 4.4 style)")
    print("  - comparison_2D_N41.png")
    print("  - comparison_3D_N41.png")
    print("="*70 + "\n")

    return {
        'N21': {'solver': solver_21, 'errors': errors_21, 'result': result_21},
        'N41': {'solver': solver_41, 'errors': errors_41, 'result': result_41}
    }


if __name__ == "__main__":
    results = main()
