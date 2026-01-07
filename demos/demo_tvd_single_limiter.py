"""
Demo 2: Single TVD Limiter Comparison with Analytical Solution.

Allows user to select ONE TVD flux limiter and compare with analytical solution.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import os

# Add project root to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.numerical.solver_fvm_tvd import BurgersSolverFVM_TVD, FluxLimiter, SolverMethod
from src.analytical.case1_solution import analytical_solution_case1


# Available limiters with descriptions
LIMITERS = {
    '1': (FluxLimiter.MINMOD, 'Minmod', 'Most diffusive, very stable'),
    '2': (FluxLimiter.VAN_ALBADA, 'Van Albada', 'Smooth, balanced accuracy/stability'),
    '3': (FluxLimiter.UMIST, 'UMIST', 'Sharp shock capturing'),
    '4': (FluxLimiter.SUPERBEE, 'Superbee', 'Most compressive'),
    '5': (FluxLimiter.VAN_LEER, 'Van Leer', 'Good balance'),
    '6': (FluxLimiter.KOREN, 'Koren', 'Third-order accurate'),
    '7': (FluxLimiter.QUICK, 'QUICK', 'Quadratic upstream interpolation'),
    '8': (FluxLimiter.SMART, 'SMART', 'High resolution'),
}


def display_menu():
    """Display the limiter selection menu."""
    print("Select a TVD Flux Limiter:")
    for key, (_, name, desc) in LIMITERS.items():
        print(f"  [{key}] {name:12} - {desc}")


def get_user_choice():
    """Get user's limiter choice."""
    while True:
        choice = input("\nEnter your choice (1-8): ").strip()
        if choice in LIMITERS:
            return LIMITERS[choice]
        print("Invalid choice. Please enter a number between 1 and 8.")


def solve_with_limiter(limiter: FluxLimiter, N: int = 41, 
                       max_iterations: int = 10000, tolerance: float = 1e-6):
    """
    Solve the 2D Inviscid Burgers Equation using the specified limiter.
    
    Parameters
    ----------
    limiter : FluxLimiter
        The flux limiter to use
    N : int
        Grid size (N x N)
    max_iterations : int
        Maximum iterations for convergence
    tolerance : float
        Convergence tolerance
        
    Returns
    -------
    solver : BurgersSolverFVM_TVD
        The solver object with solution
    converged : bool
        Whether the solution converged
    iterations : int
        Number of iterations used
    """
    print(f"Solving with {limiter.value} limiter on {N}×{N} grid...")
    
    # Create solver
    solver = BurgersSolverFVM_TVD(
        N=N, a=1.0, b=1.0, c=1.5, d=-0.5,
        cfl=0.5, limiter=limiter
    )
    
    # Solve to steady state
    converged, iterations, final_residual = solver.solve_steady_state(
        max_iterations=max_iterations,
        tolerance=tolerance,
        print_interval=1000,
        method=SolverMethod.EXPLICIT_EULER
    )
    
    return solver, converged, iterations


def compute_errors(solver, X, Y):
    """Compute error metrics between numerical and analytical solutions."""
    U_numerical = solver.u[1:-1, 1:-1]  # Interior points
    U_analytical = analytical_solution_case1(X, Y)
    
    error = np.abs(U_numerical - U_analytical)
    
    L1_error = np.mean(error)
    L2_error = np.sqrt(np.mean(error**2))
    Linf_error = np.max(error)
    
    return L1_error, L2_error, Linf_error, error


def plot_comparison(solver, limiter_name: str, save_path: str = None):
    """
    Plot side-by-side comparison of numerical and analytical solutions.
    
    Parameters
    ----------
    solver : BurgersSolverFVM_TVD
        The solver object with solution
    limiter_name : str
        Name of the limiter used
    save_path : str, optional
        Path to save the figure
    """
    # Get grid
    x = solver.x
    y = solver.y
    X, Y = np.meshgrid(x, y)
    
    # Get solutions
    U_numerical = solver.u[1:-1, 1:-1]
    U_analytical = analytical_solution_case1(X, Y)
    
    # Compute error
    L1, L2, Linf, error = compute_errors(solver, X, Y)
    
    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    levels = np.linspace(-0.5, 1.5, 21)
    
    # 1. Numerical Solution
    cf1 = axes[0].contourf(X, Y, U_numerical, levels=levels, cmap='RdYlBu_r')
    axes[0].contour(X, Y, U_numerical, levels=levels, colors='black', linewidths=0.3, alpha=0.5)
    plt.colorbar(cf1, ax=axes[0], label='u(x, y)')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    axes[0].set_title(f'Numerical Solution\n({limiter_name} Limiter)')
    axes[0].set_aspect('equal')
    
    # 2. Analytical Solution
    cf2 = axes[1].contourf(X, Y, U_analytical, levels=levels, cmap='RdYlBu_r')
    axes[1].contour(X, Y, U_analytical, levels=levels, colors='black', linewidths=0.3, alpha=0.5)
    plt.colorbar(cf2, ax=axes[1], label='u(x, y)')
    axes[1].set_xlabel('x')
    axes[1].set_ylabel('y')
    axes[1].set_title('Analytical Solution')
    axes[1].set_aspect('equal')
    
    # 3. Error Distribution
    cf3 = axes[2].contourf(X, Y, error, levels=20, cmap='hot_r')
    plt.colorbar(cf3, ax=axes[2], label='|Error|')
    axes[2].set_xlabel('x')
    axes[2].set_ylabel('y')
    axes[2].set_title(f'Error Distribution\nL₁={L1:.4f}, L₂={L2:.4f}, L∞={Linf:.4f}')
    axes[2].set_aspect('equal')
    
    plt.suptitle(f'2D Inviscid Burgers Equation - {limiter_name} Limiter vs Analytical',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"\nFigure saved to: {save_path}")
    
    return fig


def main():
    """Main function for single limiter demonstration."""
    print("TVD Solver Demo")
    
    # Display menu and get choice
    display_menu()
    limiter_enum, limiter_name, limiter_desc = get_user_choice()
    
    print(f"Selected: {limiter_name}")
    
    # Solve
    solver, converged, iterations = solve_with_limiter(limiter_enum, N=41)
    
    if converged:
        print(f"Converged in {iterations} iterations")
    else:
        print(f"Did not converge after {iterations} iterations")
    
    # Compute and display errors
    X, Y = np.meshgrid(solver.x, solver.y)
    L1, L2, Linf, _ = compute_errors(solver, X, Y)
    
    print(f"Errors: L1={L1:.6f}, L2={L2:.6f}, Linf={Linf:.6f}")
    
    # Plot comparison
    save_name = limiter_name.lower().replace(' ', '_')
    fig = plot_comparison(
        solver, limiter_name,
        save_path=f'plots/tvd_comparison/demo_{save_name}_comparison.png'
    )
    
    print("Done.")
    
    plt.show()


if __name__ == "__main__":
    main()
