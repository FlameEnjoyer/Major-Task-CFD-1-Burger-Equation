"""
Analytical Solution for 2D Inviscid Burgers Equation - Case 1

This module implements the piecewise analytical solution for:
    ∂u/∂t + a*u*∂u/∂x + b*∂u/∂y = 0

Case 1 Parameters:
    a = 1.0, b = 1.0, c = 1.5, d = -0.5

Boundary Conditions:
    u(0, y) = c = 1.5
    u(1, y) = d = -0.5
    u(x, 0) = c - (c - d)*x = 1.5 - 2.0*x

Domain: 0 ≤ x ≤ 1, 0 ≤ y ≤ 1
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from typing import Tuple, Union


def analytical_solution_case1(x: Union[float, np.ndarray],
                               y: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Compute the analytical solution for Case 1 of the 2D Inviscid Burgers Equation.

    Parameters
    ----------
    a = 1.0, b = 1.0, c = 1.5, d = -0.5

    Piecewise solution:

    For y < 0.5:
        u(x,y) = 1.5                    if x ≤ 1.5y
        u(x,y) = (1.5 - 2.0x)/(1 - 2y)  if 1.5y ≤ x ≤ (1 - 0.5y)
        u(x,y) = -0.5                   if x ≥ (1 - 0.5y)

    For y ≥ 0.5:
        u(x,y) = 1.5                    if x ≤ (0.5 + 0.5y)
        u(x,y) = -0.5                   if x > (0.5 + 0.5y)

    Parameters
    ----------
    x : float or np.ndarray
        x-coordinate(s) in the domain [0, 1]
    y : float or np.ndarray
        y-coordinate(s) in the domain [0, 1]

    Returns
    -------
    u : float or np.ndarray
        Solution value(s) at the given point(s)

    Examples
    --------
    >>> # Single point evaluation
    >>> u = analytical_solution_case1(0.5, 0.3)

    >>> # Mesh evaluation
    >>> x = np.linspace(0, 1, 100)
    >>> y = np.linspace(0, 1, 100)
    >>> X, Y = np.meshgrid(x, y)
    >>> U = analytical_solution_case1(X, Y)
    """
    # Parameters
    c = 1.5
    d = -0.5

    # Convert to numpy arrays for vectorized operations
    x = np.asarray(x)
    y = np.asarray(y)

    # Ensure arrays are broadcastable and get the output shape
    x_arr, y_arr = np.broadcast_arrays(x, y)

    # Initialize solution array with the correct shape
    u = np.zeros_like(x_arr, dtype=float)

    # Case 1: y < 0.5
    mask_y_low = y_arr < 0.5

    if np.any(mask_y_low):
        # Region 1: x ≤ 1.5y
        mask_region1 = mask_y_low & (x_arr <= 1.5 * y_arr)
        u[mask_region1] = c  # u = 1.5

        # Region 2: 1.5y ≤ x ≤ (1 - 0.5y)
        boundary1 = 1.5 * y_arr
        boundary2 = 1.0 - 0.5 * y_arr
        mask_region2 = mask_y_low & (x_arr >= boundary1) & (x_arr <= boundary2)

        # Avoid division by zero
        denominator = 1.0 - 2.0 * y_arr[mask_region2]
        # For y very close to 0.5, handle carefully
        safe_denom = np.where(np.abs(denominator) < 1e-10, 1e-10, denominator)
        u[mask_region2] = (1.5 - 2.0 * x_arr[mask_region2]) / safe_denom

        # Region 3: x ≥ (1 - 0.5y)
        mask_region3 = mask_y_low & (x_arr >= boundary2)
        u[mask_region3] = d  # u = -0.5

    # Case 2: y ≥ 0.5
    mask_y_high = y_arr >= 0.5

    if np.any(mask_y_high):
        # Region 1: x < (0.5 + 0.5y)
        boundary = 0.5 + 0.5 * y_arr
        mask_region1 = mask_y_high & (x_arr < boundary)
        u[mask_region1] = c  # u = 1.5

        # Region 2: x >= (0.5 + 0.5y)
        # Note: Using >= to ensure boundary condition u(1,y) = -0.5 is satisfied
        mask_region2 = mask_y_high & (x_arr >= boundary)
        u[mask_region2] = d  # u = -0.5

    # Return scalar if inputs were scalars, otherwise return array
    return u if u.shape else u.item()


def verify_boundary_conditions(N: int = 100) -> bool:
    """
    Verify that the analytical solution satisfies the boundary conditions.

    Boundary conditions:
        u(0, y) = 1.5
        u(1, y) = -0.5
        u(x, 0) = 1.5 - 2.0*x

    Parameters
    ----------
    N : int, optional
        Number of points to test along each boundary (default: 100)

    Returns
    -------
    bool
        True if all boundary conditions are satisfied within tolerance
    """
    tolerance = 1e-10

    # Test points
    y_test = np.linspace(0, 1, N)
    x_test = np.linspace(0, 1, N)

    # BC 1: u(0, y) = 1.5
    u_left = analytical_solution_case1(0.0, y_test)
    bc1_satisfied = np.allclose(u_left, 1.5, atol=tolerance)

    # BC 2: u(1, y) = -0.5
    u_right = analytical_solution_case1(1.0, y_test)
    bc2_satisfied = np.allclose(u_right, -0.5, atol=tolerance)

    # BC 3: u(x, 0) = 1.5 - 2.0*x
    u_bottom = analytical_solution_case1(x_test, 0.0)
    u_bottom_expected = 1.5 - 2.0 * x_test
    bc3_satisfied = np.allclose(u_bottom, u_bottom_expected, atol=tolerance)

    print("Boundary Condition Verification:")
    print(f"  BC1: u(0, y) = 1.5        -> {'PASS' if bc1_satisfied else 'FAIL'}")
    print(f"  BC2: u(1, y) = -0.5       -> {'PASS' if bc2_satisfied else 'FAIL'}")
    print(f"  BC3: u(x, 0) = 1.5 - 2x   -> {'PASS' if bc3_satisfied else 'FAIL'}")

    return bc1_satisfied and bc2_satisfied and bc3_satisfied


def plot_analytical_solution(N: int = 100,
                             save_path: str = None,
                             show_plot: bool = True) -> Tuple[plt.Figure, np.ndarray]:
    """
    Generate a 2D contour plot of the analytical solution on an N×N mesh.

    Parameters
    ----------
    N : int, optional
        Number of grid points in each direction (default: 100)
    save_path : str, optional
        Path to save the figure. If None, figure is not saved.
    show_plot : bool, optional
        Whether to display the plot (default: True)

    Returns
    -------
    fig : matplotlib.figure.Figure
        The generated figure
    U : np.ndarray
        The solution values on the mesh (shape: N×N)

    Examples
    --------
    >>> fig, U = plot_analytical_solution(N=100, save_path='plots/analytical/case1.png')
    """
    # Create mesh
    x = np.linspace(0, 1, N)
    y = np.linspace(0, 1, N)
    X, Y = np.meshgrid(x, y)

    # Compute solution
    U = analytical_solution_case1(X, Y)

    # Create figure with multiple subplots
    fig = plt.figure(figsize=(16, 5))

    # Subplot 1: Filled contour plot
    ax1 = fig.add_subplot(131)
    contourf = ax1.contourf(X, Y, U, levels=20, cmap='RdYlBu_r')
    ax1.set_xlabel('x', fontsize=12)
    ax1.set_ylabel('y', fontsize=12)
    ax1.set_title('Analytical Solution - Filled Contour\nCase 1: a=1, b=1, c=1.5, d=-0.5', fontsize=11)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    cbar1 = plt.colorbar(contourf, ax=ax1)
    cbar1.set_label('u(x, y)', fontsize=11)

    # Subplot 2: Contour lines
    ax2 = fig.add_subplot(132)
    contour = ax2.contour(X, Y, U, levels=15, colors='black', linewidths=0.5)
    ax2.clabel(contour, inline=True, fontsize=8, fmt='%.2f')
    contourf2 = ax2.contourf(X, Y, U, levels=20, cmap='RdYlBu_r', alpha=0.6)
    ax2.set_xlabel('x', fontsize=12)
    ax2.set_ylabel('y', fontsize=12)
    ax2.set_title('Analytical Solution - Contour Lines', fontsize=11)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)

    # Mark shock boundaries
    y_plot = np.linspace(0, 0.5, 100)
    x_shock1 = 1.5 * y_plot  # Left shock for y < 0.5
    x_shock2 = 1.0 - 0.5 * y_plot  # Right shock for y < 0.5
    ax2.plot(x_shock1, y_plot, 'r--', linewidth=2, label='Shock boundaries (y<0.5)')
    ax2.plot(x_shock2, y_plot, 'r--', linewidth=2)

    y_plot2 = np.linspace(0.5, 1, 100)
    x_shock3 = 0.5 + 0.5 * y_plot2  # Shock for y ≥ 0.5
    ax2.plot(x_shock3, y_plot2, 'b--', linewidth=2, label='Shock boundary (y≥0.5)')
    ax2.legend(fontsize=9)

    # Subplot 3: 3D surface plot
    ax3 = fig.add_subplot(133, projection='3d')
    surf = ax3.plot_surface(X, Y, U, cmap='RdYlBu_r', edgecolor='none', alpha=0.9)
    ax3.set_xlabel('x', fontsize=10)
    ax3.set_ylabel('y', fontsize=10)
    ax3.set_zlabel('u(x, y)', fontsize=10)
    ax3.set_title('Analytical Solution - 3D Surface', fontsize=11)
    ax3.view_init(elev=25, azim=45)
    cbar3 = plt.colorbar(surf, ax=ax3, shrink=0.5, aspect=5)
    cbar3.set_label('u(x, y)', fontsize=9)

    plt.tight_layout()

    # Save figure if path provided
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")

    # Show plot if requested
    if show_plot:
        plt.show()

    return fig, U


def plot_solution_slices(N: int = 100,
                         save_path: str = None,
                         show_plot: bool = True) -> plt.Figure:
    """
    Plot solution profiles along constant x and y lines.

    Parameters
    ----------
    N : int, optional
        Number of grid points (default: 100)
    save_path : str, optional
        Path to save the figure
    show_plot : bool, optional
        Whether to display the plot (default: True)

    Returns
    -------
    fig : matplotlib.figure.Figure
        The generated figure
    """
    x = np.linspace(0, 1, N)
    y = np.linspace(0, 1, N)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Plot 1: u vs x at different y values
    y_values = [0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9]
    for y_val in y_values:
        u = analytical_solution_case1(x, y_val)
        ax1.plot(x, u, label=f'y = {y_val}', linewidth=2)

    ax1.set_xlabel('x', fontsize=12)
    ax1.set_ylabel('u(x, y)', fontsize=12)
    ax1.set_title('Solution Profiles at Different y Values', fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9)
    ax1.set_xlim([0, 1])

    # Plot 2: u vs y at different x values
    x_values = [0.1, 0.3, 0.5, 0.7, 0.9]
    for x_val in x_values:
        u = analytical_solution_case1(x_val, y)
        ax2.plot(y, u, label=f'x = {x_val}', linewidth=2)

    ax2.set_xlabel('y', fontsize=12)
    ax2.set_ylabel('u(x, y)', fontsize=12)
    ax2.set_title('Solution Profiles at Different x Values', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9)
    ax2.set_xlim([0, 1])

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")

    if show_plot:
        plt.show()

    return fig


def print_solution_summary(N: int = 100):
    """
    Print summary statistics of the analytical solution.

    Parameters
    ----------
    N : int, optional
        Number of grid points for evaluation (default: 100)
    """
    # Create mesh
    x = np.linspace(0, 1, N)
    y = np.linspace(0, 1, N)
    X, Y = np.meshgrid(x, y)

    # Compute solution
    U = analytical_solution_case1(X, Y)

    print("\n" + "="*60)
    print("Analytical Solution Summary - Case 1")
    print("="*60)
    print(f"Parameters: a = 1.0, b = 1.0, c = 1.5, d = -0.5")
    print(f"Domain: [0, 1] × [0, 1]")
    print(f"Grid size: {N} × {N}")
    print("-"*60)
    print(f"Solution statistics:")
    print(f"  Min value:  {U.min():.6f}")
    print(f"  Max value:  {U.max():.6f}")
    print(f"  Mean value: {U.mean():.6f}")
    print(f"  Std dev:    {U.std():.6f}")
    print("-"*60)

    # Verify boundary conditions
    verify_boundary_conditions(N)
    print("="*60 + "\n")


if __name__ == "__main__":
    # Example usage and testing
    print("Testing Analytical Solution - Case 1\n")

    # Test single point
    print("Single point evaluation:")
    u_test = analytical_solution_case1(0.5, 0.3)
    print(f"u(0.5, 0.3) = {u_test:.6f}\n")

    # Print summary
    print_solution_summary(N=100)

    # Generate plots
    print("Generating plots...")
    fig, U = plot_analytical_solution(
        N=100,
        save_path='../../plots/analytical/case1_contour.png',
        show_plot=False
    )

    fig_slices = plot_solution_slices(
        N=100,
        save_path='../../plots/analytical/case1_slices.png',
        show_plot=False
    )

    print("\nPlots generated successfully!")
    print("  - case1_contour.png")
    print("  - case1_slices.png")
