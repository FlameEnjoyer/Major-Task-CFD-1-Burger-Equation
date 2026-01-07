"""
Generate Side-by-Side Comparisons of Analytical and Numerical Solutions.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
import yaml
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.numerical.solver import BurgersSolver2D
from src.analytical.case1_solution import analytical_solution_case1


def plot_2d_comparison(solver, N, save_path=None):
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
    fig.suptitle(f'Comparison of Analytical and Upwind FVM Solutions in 2D Contour with {N} Grid',
                 fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 0.88, 0.96])

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"[SAVED] {save_path}")

    plt.close()

    return fig

def plot_error_field(solver, N, save_path=None):
    """Figure B: Standalone Error Field Plot (Log-Scaled)."""
    # Get solutions
    X_num = solver.X
    Y_num = solver.Y
    U_num = solver.u
    U_analytical = analytical_solution_case1(X_num, Y_num)

    # Compute Squared Error
    Squared_Error = (U_num - U_analytical)**2
    
    # Avoid log(0)
    epsilon = 1e-16
    Squared_Error_Log = Squared_Error + epsilon

    # Setup Figure
    fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)

    # Log Normalization for visibility
    # We clip the bottom range at 1e-10 to hide machine noise
    err_vmin = max(Squared_Error.min(), 1e-10)
    err_vmax = Squared_Error.max()

    contourf = ax.contourf(X_num, Y_num, Squared_Error,
                           levels=25,
                           cmap='jet')

    ax.set_xlabel('$x$', fontsize=14)
    ax.set_ylabel('$y$', fontsize=14)
    ax.set_title(f'Squared Error Field $(u_{{num}} - u_{{ana}})^2$\n(Normal Scale, N={N})', 
                 fontsize=16, fontweight='bold')
    ax.set_aspect('equal')

    # Colorbar
    cbar = fig.colorbar(contourf, ax=ax, aspect=20, pad=0.05)
    cbar.set_label('Squared Error', fontsize=14)

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"[SAVED] {save_path}")

    plt.close()
    return fig

def plot_3d_comparison(solver, N, save_path=None):
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
    fig.suptitle(f'Comparison of Analytical and Upwind FVM Solutions in 3D Surface with {N} Grid',
                 fontsize=14, fontweight='bold', y=0.95)

    plt.tight_layout(rect=[0, 0, 0.88, 0.93])

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"[SAVED] {save_path}")

    plt.close()

    return fig



def plot_convergence_history(result, N, save_path=None):
    """
    Plot the residual error convergence history.
    """
    residual_history = result['residual_history']
    iterations = range(1, len(residual_history) + 1)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogy(iterations, residual_history, 'b-', linewidth=1.5)
    ax.set_xlabel('Iteration', fontsize=12)
    ax.set_ylabel('Residual Error (L2 Norm)', fontsize=12)
    ax.set_title(f'Residual Convergence History (N={N})', fontsize=14, fontweight='bold')
    ax.grid(True, which='both', linestyle='--', alpha=0.7)
    
    # Mark final value
    if residual_history:
        ax.axhline(y=residual_history[-1], color='r', linestyle='--', alpha=0.5, 
                  label=f'Final Residual: {residual_history[-1]:.2e}')
        ax.legend(fontsize=11)

    plt.tight_layout()

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

    # Generate Error Field plot
    print(f"\n2. Generating Error Field plot...")
    save_path_error = f'{output_dir}/error_field_N{N}.png'
    plot_error_field(solver, N, save_path=save_path_error)

    # Generate 3D comparison plot
    print(f"\n3. Generating 3D surface comparison plot...")
    save_path_3d = f'{output_dir}/comparison_3D_N{N}.png'
    plot_3d_comparison(solver, N, save_path=save_path_3d)

    # Generate Convergence plot
    print(f"\n4. Generating Residual Convergence plot...")
    save_path_conv = f'{output_dir}/convergence_N{N}.png'
    plot_convergence_history(result, N, save_path=save_path_conv)

    print(f"\n" + "="*70)
    print(f" COMPLETED COMPARISON FOR N={N}")
    print("="*70)

    return solver, errors, result


def main():
    """Main function to generate all comparisons."""
    print("Starting analytical vs numerical comparison...")

    output_dir = 'plots/comparisons'
    grid_sizes = [21, 41, 201]
    
    comparisons_history = {}

    for N in grid_sizes:
        print("\n\n")
        solver, errors, result = run_comparison_for_grid_size(N=N, output_dir=output_dir)
        
        # Structure data for YAML
        comparisons_history[f'N_{N}'] = {
            'simulation_results': {
                'converged': bool(result['converged']),
                'iterations': int(result['iterations']),
                'final_residual': float(result['final_residual']),
                'time_elapsed': float(result['time_elapsed'])
            },
            'error_metrics': {
                'L1_error': float(errors['L1']),
                'L2_error': float(errors['L2']),
                'Linf_error': float(errors['Linf']),
                'relative_L2': float(errors['relative_L2'])
            }
        }

    # Save history to YAML
    yaml_path = f'{output_dir}/case1_comparison_history.yaml'
    with open(yaml_path, 'w') as f:
        yaml.dump({'CASE1_COMPARISONS': comparisons_history}, f, sort_keys=False)

    # Summary
    print("\n\n" + "="*70)
    print(" FINAL SUMMARY")
    print("="*70)
    
    for N in grid_sizes:
        res = comparisons_history[f'N_{N}']['simulation_results']
        err = comparisons_history[f'N_{N}']['error_metrics']
        print(f"\nComparison for N={N}:")
        print(f"  Converged: {res['converged']}, Iterations: {res['iterations']}")
        print(f"  L2 Error: {err['L2_error']:.6e}")
        print(f"  Time: {res['time_elapsed']:.2f} s")

    print("\n" + "="*70)
    print(" ALL PLOTS GENERATED SUCCESSFULLY!")
    print("="*70)
    print(f"\nOutput files saved in {output_dir}/:")
    for N in grid_sizes:
        print(f"  - comparison_2D_N{N}.png")
        print(f"  - comparison_3D_N{N}.png")
        print(f"  - error_field_N{N}.png")
        print(f"  - convergence_N{N}.png")
    print(f"  - case1_comparison_history.yaml")
    print("="*70 + "\n")

    return comparisons_history

if __name__ == "__main__":
    main()
