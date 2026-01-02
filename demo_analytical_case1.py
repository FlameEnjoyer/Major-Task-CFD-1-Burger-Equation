"""
Demonstration Script for Analytical Solution - Case 1
2D Inviscid Burgers Equation

This script demonstrates the analytical solution implementation for Case 1
with parameters: a=1, b=1, c=1.5, d=-0.5
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from analytical.case1_solution import (
    analytical_solution_case1,
    plot_analytical_solution,
    plot_solution_slices,
    print_solution_summary,
    verify_boundary_conditions
)


def main():
    """Main demonstration function."""

    print("\n" + "="*70)
    print(" " * 15 + "2D INVISCID BURGERS EQUATION")
    print(" " * 15 + "Analytical Solution - Case 1")
    print("="*70)
    print("\nEquation: du/dt + a*u*du/dx + b*du/dy = 0")
    print("\nParameters:")
    print("  a = 1.0")
    print("  b = 1.0")
    print("  c = 1.5 (left boundary)")
    print("  d = -0.5 (right boundary)")
    print("\nDomain: [0, 1] × [0, 1]")
    print("="*70 + "\n")

    # Test individual points
    print("Testing individual points:")
    print("-" * 70)

    test_points = [
        (0.2, 0.2),
        (0.5, 0.3),
        (0.7, 0.4),
        (0.3, 0.6),
        (0.8, 0.8)
    ]

    for x, y in test_points:
        u = analytical_solution_case1(x, y)
        print(f"  u({x:.1f}, {y:.1f}) = {u:8.4f}")

    print("-" * 70 + "\n")

    # Verify boundary conditions
    print("Verifying boundary conditions...")
    print("-" * 70)
    verify_boundary_conditions(N=100)
    print("-" * 70 + "\n")

    # Print solution summary
    print_solution_summary(N=100)

    # Generate visualizations
    print("\nGenerating visualizations...")
    print("-" * 70)

    # Main contour plot with 100x100 mesh
    print("\n1. Creating contour plots (100×100 mesh)...")
    fig1, U = plot_analytical_solution(
        N=100,
        save_path='plots/analytical/case1_contour_100x100.png',
        show_plot=False
    )
    print(f"   Saved: plots/analytical/case1_contour_100x100.png")

    # Solution slices
    print("\n2. Creating solution profile plots...")
    fig2 = plot_solution_slices(
        N=100,
        save_path='plots/analytical/case1_slices.png',
        show_plot=False
    )
    print(f"   Saved: plots/analytical/case1_slices.png")

    # High resolution contour
    print("\n3. Creating high-resolution contour (200×200 mesh)...")
    fig3, U_hires = plot_analytical_solution(
        N=200,
        save_path='plots/analytical/case1_contour_200x200.png',
        show_plot=False
    )
    print(f"   Saved: plots/analytical/case1_contour_200x200.png")

    print("-" * 70)

    # Additional analysis
    print("\nAdditional Analysis:")
    print("-" * 70)

    # Analyze shock regions
    print("\nShock Wave Boundaries:")
    print("  For y < 0.5:")
    print("    Left shock:  x = 1.5y")
    print("    Right shock: x = 1 - 0.5y")
    print("  For y >= 0.5:")
    print("    Single shock: x = 0.5 + 0.5y")

    # Check solution values in different regions
    print("\nSolution values in different regions:")
    print(f"  Region 1 (left):   u = {1.5:.1f}")
    print(f"  Region 2 (middle): u varies from {1.5:.1f} to {-0.5:.1f}")
    print(f"  Region 3 (right):  u = {-0.5:.1f}")

    print("\n" + "="*70)
    print("Demonstration completed successfully!")
    print("="*70 + "\n")

    # Show all plots
    print("Displaying plots...")
    plt.show()


if __name__ == "__main__":
    # Create plots directory if it doesn't exist
    os.makedirs('plots/analytical', exist_ok=True)

    # Run demonstration
    main()
