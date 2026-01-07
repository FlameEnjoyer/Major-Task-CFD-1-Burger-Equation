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
from typing import Union


def analytical_solution_case1(x: Union[float, np.ndarray],
                               y: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    Compute the analytical solution for Case 1 of the 2D Inviscid Burgers Equation.

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
    >>> u = analytical_solution_case1(0.5, 0.3)
    >>> X, Y = np.meshgrid(np.linspace(0, 1, 100), np.linspace(0, 1, 100))
    >>> U = analytical_solution_case1(X, Y)
    """
    c = 1.5
    d = -0.5

    x = np.asarray(x)
    y = np.asarray(y)

    x_arr, y_arr = np.broadcast_arrays(x, y)
    u = np.zeros_like(x_arr, dtype=float)

    # For y < 0.5
    mask_y_low = y_arr < 0.5

    if np.any(mask_y_low):
        # Region 1: x <= 1.5y
        mask_region1 = mask_y_low & (x_arr <= 1.5 * y_arr)
        u[mask_region1] = c

        # Region 2: 1.5y <= x <= (1 - 0.5y)
        boundary1 = 1.5 * y_arr
        boundary2 = 1.0 - 0.5 * y_arr
        mask_region2 = mask_y_low & (x_arr >= boundary1) & (x_arr <= boundary2)

        denominator = 1.0 - 2.0 * y_arr[mask_region2]
        safe_denom = np.where(np.abs(denominator) < 1e-10, 1e-10, denominator)
        u[mask_region2] = (1.5 - 2.0 * x_arr[mask_region2]) / safe_denom

        # Region 3: x >= (1 - 0.5y)
        mask_region3 = mask_y_low & (x_arr >= boundary2)
        u[mask_region3] = d

    # For y >= 0.5
    mask_y_high = y_arr >= 0.5

    if np.any(mask_y_high):
        boundary = 0.5 + 0.5 * y_arr

        # Region 1: x < (0.5 + 0.5y)
        mask_region1 = mask_y_high & (x_arr < boundary)
        u[mask_region1] = c

        # Region 2: x >= (0.5 + 0.5y)
        mask_region2 = mask_y_high & (x_arr >= boundary)
        u[mask_region2] = d

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

    y_test = np.linspace(0, 1, N)
    x_test = np.linspace(0, 1, N)

    # BC 1
    u_left = analytical_solution_case1(0.0, y_test)
    bc1_satisfied = np.allclose(u_left, 1.5, atol=tolerance)

    # BC 2
    u_right = analytical_solution_case1(1.0, y_test)
    bc2_satisfied = np.allclose(u_right, -0.5, atol=tolerance)

    # BC 3
    u_bottom = analytical_solution_case1(x_test, 0.0)
    u_bottom_expected = 1.5 - 2.0 * x_test
    bc3_satisfied = np.allclose(u_bottom, u_bottom_expected, atol=tolerance)

    print("Boundary Condition Verification:")
    print(f"  BC1: {'PASS' if bc1_satisfied else 'FAIL'}")
    print(f"  BC2: {'PASS' if bc2_satisfied else 'FAIL'}")
    print(f"  BC3: {'PASS' if bc3_satisfied else 'FAIL'}")

    return bc1_satisfied and bc2_satisfied and bc3_satisfied
