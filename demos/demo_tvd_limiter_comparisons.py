"""
Generate Comprehensive TVD Limiter Comparisons.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.numerical.solver_fvm_tvd import BurgersSolverFVM_TVD, FluxLimiter
from src.analytical.case1_solution import analytical_solution_case1


def plot_analytical_vs_numerical_2d(solver, limiter_name, N, save_path=None):
    """Create side-by-side 2D contour comparison plot."""
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
    fig.suptitle(f'Comparison of Analytical and TVD {limiter_name} FVM Solutions in 2D Contour with {N} Grid',
                 fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 0.88, 0.96])

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"[SAVED] {save_path}")

    plt.close()

    return fig


def plot_analytical_vs_numerical_3d(solver, limiter_name, N, save_path=None):
    """Create side-by-side 3D surface comparison plot."""
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
    fig.suptitle(f'Comparison of Analytical and TVD {limiter_name} FVM Solutions in 3D Surface with {N} Grid',
                 fontsize=14, fontweight='bold', y=0.95)

    plt.tight_layout(rect=[0, 0, 0.88, 0.93])

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"[SAVED] {save_path}")

    plt.close()

    return fig


def plot_error_comparison(solver_21, solver_41, limiter_name, save_path=None):
    """Create side-by-side error plots for N=21 and N=41."""
    # Compute errors for N=21
    X_21 = solver_21.X
    Y_21 = solver_21.Y
    U_num_21 = solver_21.u
    U_analytical_21 = analytical_solution_case1(X_21, Y_21)
    error_21 = U_num_21 - U_analytical_21

    # Compute errors for N=41
    X_41 = solver_41.X
    Y_41 = solver_41.Y
    U_num_41 = solver_41.u
    U_analytical_41 = analytical_solution_case1(X_41, Y_41)
    error_41 = U_num_41 - U_analytical_41

    # Create figure with 2 subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # Common colorbar limits (use symmetric range)
    vmin = min(error_21.min(), error_41.min())
    vmax = max(error_21.max(), error_41.max())
    abs_max = max(abs(vmin), abs(vmax))
    levels = np.linspace(-abs_max, abs_max, 25)

    # Left subplot: Error for N=21
    contourf1 = ax1.contourf(X_21, Y_21, error_21, levels=levels,
                             cmap='RdYlBu_r', extend='both')
    ax1.set_xlabel('$x$', fontsize=13)
    ax1.set_ylabel('$y$', fontsize=13)
    ax1.set_title(f'Velocity (u) Differences in TVD {limiter_name} FVM\nto Analytical with 21 grid',
                  fontsize=12, fontweight='bold')
    ax1.set_aspect('equal')
    ax1.set_xlim([0, 1])
    ax1.set_ylim([0, 1])

    # Right subplot: Error for N=41
    contourf2 = ax2.contourf(X_41, Y_41, error_41, levels=levels,
                             cmap='RdYlBu_r', extend='both')
    ax2.set_xlabel('$x$', fontsize=13)
    ax2.set_ylabel('$y$', fontsize=13)
    ax2.set_title(f'Velocity (u) Differences in TVD {limiter_name} FVM\nto Analytical with 41 grid',
                  fontsize=12, fontweight='bold')
    ax2.set_aspect('equal')
    ax2.set_xlim([0, 1])
    ax2.set_ylim([0, 1])

    # Add shared colorbar
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(contourf2, cax=cbar_ax)
    cbar.set_label('Error', fontsize=13, rotation=90, labelpad=15)

    # Overall title
    fig.suptitle(f'Contour Plot of error',
                 fontsize=16, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 0.88, 0.96])

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"[SAVED] {save_path}")

    plt.close()

    return fig


def run_limiter_comparison(limiter, limiter_name, output_dir='plots/tvd_comparisons'):
    """
    Run complete comparison for a given flux limiter.

    Parameters
    ----------
    limiter : FluxLimiter
        Flux limiter enum value
    limiter_name : str
        Display name for the limiter
    output_dir : str
        Directory to save output plots
    """
    print("\n" + "="*80)
    print(f" RUNNING COMPARISON FOR TVD {limiter_name}")
    print("="*80)

    # Create output directory
    limiter_dir = os.path.join(output_dir, limiter_name.lower().replace(' ', '_').replace('-', '_'))
    os.makedirs(limiter_dir, exist_ok=True)

    # -------------------------------------------------------------------------
    # 1. Run solver for N=41 (for analytical comparison)
    # -------------------------------------------------------------------------
    print(f"\n1. Solving with TVD {limiter_name} for N=41...")
    solver_41 = BurgersSolverFVM_TVD(
        N=41,
        a=1.0,
        b=1.0,
        c=1.5,
        d=-0.5,
        cfl=0.5,
        limiter=limiter
    )

    result_41 = solver_41.solve_steady_state(
        max_iterations=20000,
        tolerance=1e-6,
        print_interval=5000
    )

    errors_41 = solver_41.compute_error(analytical_solution_case1)

    print(f"\nN=41 Results:")
    print(f"  Converged: {result_41['converged']}")
    print(f"  Iterations: {result_41['iterations']}")
    print(f"  Time: {result_41['time_elapsed']:.2f} s")
    print(f"  L1 error:  {errors_41['L1']:.6e}")
    print(f"  L2 error:  {errors_41['L2']:.6e}")
    print(f"  Linf error: {errors_41['Linf']:.6e}")

    # -------------------------------------------------------------------------
    # 2. Run solver for N=21 (for error comparison)
    # -------------------------------------------------------------------------
    print(f"\n2. Solving with TVD {limiter_name} for N=21...")
    solver_21 = BurgersSolverFVM_TVD(
        N=21,
        a=1.0,
        b=1.0,
        c=1.5,
        d=-0.5,
        cfl=0.5,
        limiter=limiter
    )

    result_21 = solver_21.solve_steady_state(
        max_iterations=20000,
        tolerance=1e-6,
        print_interval=5000
    )

    errors_21 = solver_21.compute_error(analytical_solution_case1)

    print(f"\nN=21 Results:")
    print(f"  Converged: {result_21['converged']}")
    print(f"  Iterations: {result_21['iterations']}")
    print(f"  Time: {result_21['time_elapsed']:.2f} s")
    print(f"  L1 error:  {errors_21['L1']:.6e}")
    print(f"  L2 error:  {errors_21['L2']:.6e}")
    print(f"  Linf error: {errors_21['Linf']:.6e}")

    # -------------------------------------------------------------------------
    # 3. Generate 2D comparison plot (N=41)
    # -------------------------------------------------------------------------
    print(f"\n3. Generating 2D contour comparison plot (N=41)...")
    save_path_2d = os.path.join(limiter_dir, 'comparison_2D_N41.png')
    plot_analytical_vs_numerical_2d(solver_41, limiter_name, 41, save_path=save_path_2d)

    # -------------------------------------------------------------------------
    # 4. Generate 3D comparison plot (N=41)
    # -------------------------------------------------------------------------
    print(f"\n4. Generating 3D surface comparison plot (N=41)...")
    save_path_3d = os.path.join(limiter_dir, 'comparison_3D_N41.png')
    plot_analytical_vs_numerical_3d(solver_41, limiter_name, 41, save_path=save_path_3d)

    # -------------------------------------------------------------------------
    # 5. Generate error comparison plot (N=21 and N=41)
    # -------------------------------------------------------------------------
    print(f"\n5. Generating error comparison plot (N=21 vs N=41)...")
    save_path_error = os.path.join(limiter_dir, 'error_comparison_N21_N41.png')
    plot_error_comparison(solver_21, solver_41, limiter_name, save_path=save_path_error)

    print(f"\n" + "="*80)
    print(f" COMPLETED COMPARISON FOR TVD {limiter_name}")
    print("="*80)

    return {
        'solver_21': solver_21,
        'solver_41': solver_41,
        'errors_21': errors_21,
        'errors_41': errors_41,
        'result_21': result_21,
        'result_41': result_41
    }


def main():
    """Main function to generate all TVD limiter comparisons."""
    print("Starting TVD limiter comparison...")

    output_dir = 'plots/tvd_comparisons'

    # Define limiters to compare
    limiters = [
        (FluxLimiter.UPWIND, "Upwind"),
        (FluxLimiter.MINMOD, "Min-Mod"),
        (FluxLimiter.VAN_LEER, "Van Leer"),
        (FluxLimiter.VAN_ALBADA, "Albada"),
        (FluxLimiter.SUPERBEE, "Superbee"),
        (FluxLimiter.QUICK, "QUICK"),
        (FluxLimiter.UMIST, "UMIST"),
    ]

    # Note: Sweby is not implemented in the enum, using available limiters

    results = {}

    # Run comparison for each limiter
    for limiter, limiter_name in limiters:
        results[limiter_name] = run_limiter_comparison(limiter, limiter_name, output_dir)

    # -------------------------------------------------------------------------
    # Final Summary
    # -------------------------------------------------------------------------
    print("\n\n" + "="*80)
    print(" FINAL SUMMARY - ALL TVD LIMITERS")
    print("="*80)

    print(f"\n{'Limiter':<15} {'N=21 L2 Error':<15} {'N=41 L2 Error':<15} {'Iterations (N=41)':<20}")
    print("-"*80)

    for limiter, limiter_name in limiters:
        res = results[limiter_name]
        print(f"{limiter_name:<15} {res['errors_21']['L2']:<15.6e} {res['errors_41']['L2']:<15.6e} {res['result_41']['iterations']:<20}")

    print("\n" + "="*80)
    print(" ALL PLOTS GENERATED SUCCESSFULLY!")
    print("="*80)
    print(f"\nOutput directory: {output_dir}/")
    print("Each limiter has its own subdirectory with:")
    print("  - comparison_2D_N41.png  (Analytical vs Numerical 2D)")
    print("  - comparison_3D_N41.png  (Analytical vs Numerical 3D)")
    print("  - error_comparison_N21_N41.png  (Error plots side-by-side)")
    print("="*80 + "\n")

    return results


if __name__ == "__main__":
    results = main()
