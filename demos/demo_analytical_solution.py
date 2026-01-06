"""
=================================================================================
Demo 1: Analytical Solution for 2D Inviscid Burgers Equation
=================================================================================

This script visualizes the exact analytical solution for Case 1 of the
2D Inviscid Burgers Equation.

Equation: ∂u/∂t + a*u*∂u/∂x + b*∂u/∂y = 0

Case 1 Parameters:
    a = 1.0, b = 1.0, c = 1.5, d = -0.5

Boundary Conditions:
    u(0, y) = c = 1.5   (left boundary)
    u(1, y) = d = -0.5  (right boundary)
    u(x, 0) = 1.5 - 2x  (bottom boundary)

Usage:
    python demos/demo_analytical_solution.py

=================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import os

# Add project root to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.analytical.case1_solution import analytical_solution_case1


def plot_analytical_solution(N: int = 100, save_path: str = None):
    """
    Generate visualizations of the analytical solution.
    
    Parameters
    ----------
    N : int
        Grid resolution (N x N points)
    save_path : str, optional
        Path to save the figure
    """
    # Create mesh
    x = np.linspace(0, 1, N)
    y = np.linspace(0, 1, N)
    X, Y = np.meshgrid(x, y)
    
    # Compute analytical solution
    U = analytical_solution_case1(X, Y)
    
    # Create figure with 2x2 subplots
    fig = plt.figure(figsize=(14, 12))
    
    # 1. Filled Contour Plot
    ax1 = fig.add_subplot(2, 2, 1)
    levels = np.linspace(-0.5, 1.5, 21)
    cf = ax1.contourf(X, Y, U, levels=levels, cmap='RdYlBu_r')
    ax1.contour(X, Y, U, levels=levels, colors='black', linewidths=0.3, alpha=0.5)
    plt.colorbar(cf, ax=ax1, label='u(x, y)')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_title('Analytical Solution - Filled Contour')
    ax1.set_aspect('equal')
    
    # 2. 3D Surface Plot
    ax2 = fig.add_subplot(2, 2, 2, projection='3d')
    surf = ax2.plot_surface(X, Y, U, cmap='RdYlBu_r', edgecolor='none', alpha=0.9)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_zlabel('u(x, y)')
    ax2.set_title('Analytical Solution - 3D Surface')
    ax2.view_init(elev=25, azim=45)
    
    # 3. Solution Slices at constant y
    ax3 = fig.add_subplot(2, 2, 3)
    y_slices = [0.0, 0.25, 0.5, 0.75, 1.0]
    colors = plt.cm.viridis(np.linspace(0, 1, len(y_slices)))
    for i, y_val in enumerate(y_slices):
        u_slice = analytical_solution_case1(x, np.full_like(x, y_val))
        ax3.plot(x, u_slice, color=colors[i], linewidth=2, label=f'y = {y_val:.2f}')
    ax3.set_xlabel('x')
    ax3.set_ylabel('u(x, y)')
    ax3.set_title('Solution Profiles at Constant y')
    ax3.legend(loc='best')
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0, 1)
    ax3.set_ylim(-0.7, 1.7)
    
    # 4. Solution Slices at constant x
    ax4 = fig.add_subplot(2, 2, 4)
    x_slices = [0.0, 0.25, 0.5, 0.75, 1.0]
    colors = plt.cm.plasma(np.linspace(0, 1, len(x_slices)))
    for i, x_val in enumerate(x_slices):
        u_slice = analytical_solution_case1(np.full_like(y, x_val), y)
        ax4.plot(y, u_slice, color=colors[i], linewidth=2, label=f'x = {x_val:.2f}')
    ax4.set_xlabel('y')
    ax4.set_ylabel('u(x, y)')
    ax4.set_title('Solution Profiles at Constant x')
    ax4.legend(loc='best')
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim(0, 1)
    ax4.set_ylim(-0.7, 1.7)
    
    plt.suptitle('2D Inviscid Burgers Equation - Analytical Solution (Case 1)\n'
                 'a=1, b=1, c=1.5, d=-0.5', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")
    
    return fig, U


def main():
    """Main function to demonstrate the analytical solution."""
    
    print("\n" + "=" * 70)
    print(" " * 10 + "2D INVISCID BURGERS EQUATION - ANALYTICAL SOLUTION")
    print("=" * 70)
    print("\nEquation: ∂u/∂t + a*u*∂u/∂x + b*∂u/∂y = 0")
    print("\nCase 1 Parameters:")
    print("  a = 1.0  (nonlinear convection coefficient)")
    print("  b = 1.0  (linear advection coefficient)")
    print("  c = 1.5  (left boundary value)")
    print("  d = -0.5 (right boundary value)")
    print("\nDomain: [0, 1] × [0, 1]")
    print("=" * 70)
    
    # Shock wave information
    print("\nShock Wave Structure:")
    print("-" * 70)
    print("For y < 0.5:")
    print("  - Left shock boundary:  x = 1.5y")
    print("  - Right shock boundary: x = 1 - 0.5y")
    print("  - Solution: u = 1.5 (left), varies (middle), u = -0.5 (right)")
    print("\nFor y ≥ 0.5:")
    print("  - Single shock boundary: x = 0.5 + 0.5y")
    print("  - Solution: u = 1.5 (left), u = -0.5 (right)")
    print("-" * 70)
    
    # Generate plots
    print("\nGenerating visualization...")
    fig, U = plot_analytical_solution(
        N=100, 
        save_path='plots/analytical/demo_analytical_solution.png'
    )
    
    # Print solution statistics
    print("\nSolution Statistics:")
    print(f"  Minimum value: {U.min():.4f}")
    print(f"  Maximum value: {U.max():.4f}")
    print(f"  Mean value:    {U.mean():.4f}")
    
    print("\n" + "=" * 70)
    print("Demonstration complete!")
    print("=" * 70 + "\n")
    
    plt.show()


if __name__ == "__main__":
    main()
