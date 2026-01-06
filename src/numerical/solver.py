"""
Finite Volume Method Solver for 2D Inviscid Burgers Equation

This module implements a numerical solver using:
    - Finite Volume Method (FVM) discretization
    - First-order Upwind Scheme for convective fluxes
    - Explicit Euler time integration
    - CFL-based adaptive time stepping

Equation:
    du/dt + a*u*du/dx + b*du/dy = 0

In conservative form:
    du/dt + d(a*u^2/2)/dx + d(b*u)/dy = 0

Domain: [0, 1] x [0, 1]

Boundary Conditions:
    u(0, y) = c       (left boundary)
    u(1, y) = d       (right boundary)
    u(x, 0) = c - (c-d)*x   (bottom boundary)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from typing import Tuple, Optional, Dict, List
import time


class BurgersSolver2D:
    """
    2D Inviscid Burgers Equation Solver using Finite Volume Method.

    Uses first-order upwind scheme for spatial discretization and
    explicit Euler method for time integration.

    Parameters
    ----------
    N : int
        Number of grid points in each direction (N x N grid)
    a : float
        Coefficient for nonlinear convection term (u * du/dx)
    b : float
        Coefficient for linear advection term (du/dy)
    c : float
        Left boundary condition value u(0, y) = c
    d : float
        Right boundary condition value u(1, y) = d
    cfl : float, optional
        CFL number for stability (default: 0.5)

    Attributes
    ----------
    x : np.ndarray
        x-coordinates of cell centers
    y : np.ndarray
        y-coordinates of cell centers
    X : np.ndarray
        Meshgrid of x-coordinates
    Y : np.ndarray
        Meshgrid of y-coordinates
    u : np.ndarray
        Solution array (N x N)
    dx : float
        Grid spacing in x-direction
    dy : float
        Grid spacing in y-direction
    dt : float
        Current time step

    Examples
    --------
    >>> # Case 1: a=1, b=1, c=1.5, d=-0.5
    >>> solver = BurgersSolver2D(N=41, a=1.0, b=1.0, c=1.5, d=-0.5)
    >>> solver.solve_steady_state(max_iterations=10000, tolerance=1e-6)
    >>> solver.plot_solution()
    """

    def __init__(self, N: int, a: float, b: float, c: float, d: float,
                 cfl: float = 0.5):
        """Initialize the Burgers equation solver."""

        # Store parameters
        self.N = N
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.cfl = cfl

        # Domain bounds
        self.x_min, self.x_max = 0.0, 1.0
        self.y_min, self.y_max = 0.0, 1.0

        # Grid spacing
        self.dx = (self.x_max - self.x_min) / (N - 1)
        self.dy = (self.y_max - self.y_min) / (N - 1)

        # Create grid (cell-centered values at nodes)
        self.x = np.linspace(self.x_min, self.x_max, N)
        self.y = np.linspace(self.y_min, self.y_max, N)
        self.X, self.Y = np.meshgrid(self.x, self.y)

        # Initialize solution array
        self.u = np.zeros((N, N))

        # Time step (will be computed based on CFL)
        self.dt = 0.0

        # Convergence history
        self.residual_history = []
        self.iteration_count = 0

        # Initialize with boundary and initial conditions
        self._apply_initial_conditions()
        self._apply_boundary_conditions()

        # Compute initial time step
        self._compute_time_step()

        print(f"BurgersSolver2D initialized:")
        print(f"  Grid: {N} x {N}")
        print(f"  Parameters: a={a}, b={b}, c={c}, d={d}")
        print(f"  Grid spacing: dx={self.dx:.6f}, dy={self.dy:.6f}")
        print(f"  CFL number: {cfl}")
        print(f"  Initial dt: {self.dt:.6e}")

    def _apply_initial_conditions(self):
        """
        Apply initial conditions to the solution array.

        Initial condition: u(x, 0) = c - (c - d) * x
        This is extended throughout the domain as an initial guess.
        """
        # Use bottom boundary condition as initial guess for entire domain
        for i in range(self.N):
            for j in range(self.N):
                x = self.x[i]
                self.u[j, i] = self.c - (self.c - self.d) * x

    def _apply_boundary_conditions(self):
        """
        Apply boundary conditions to the solution array.

        BCs (as specified in the PDF, Figure 4.1):
            u(0, y) = c                   (left boundary, Dirichlet)
            u(1, y) = d                   (right boundary, Dirichlet)
            u(x, 0) = c - (c-d)*x         (bottom boundary, Dirichlet)
            du/dy = 0 at y = 1            (top boundary, Neumann - zero gradient)

        Since v = b > 0, flow goes upward (from bottom to top).
        Top boundary is an outflow boundary, so we use zero gradient (Neumann).
        """
        # Left boundary: u(0, y) = c (Dirichlet)
        self.u[:, 0] = self.c

        # Right boundary: u(1, y) = d (Dirichlet)
        self.u[:, -1] = self.d

        # Bottom boundary: u(x, 0) = c - (c - d) * x (Dirichlet)
        for i in range(self.N):
            self.u[0, i] = self.c - (self.c - self.d) * self.x[i]

        # Top boundary: du/dy = 0 (Neumann - zero gradient, outflow BC)
        # Copy interior values to top row (excluding corners which are set by left/right BC)
        self.u[-1, 1:-1] = self.u[-2, 1:-1]

    def _compute_time_step(self):
        """
        Compute the time step based on the CFL condition.

        CFL condition for stability:
            dt <= CFL * min(dx / |a*u|_max, dy / |b|)

        Or equivalently:
            CFL = dt * (|a*u|_max / dx + |b| / dy) <= 1
        """
        # Maximum velocity for CFL calculation
        u_max = np.max(np.abs(self.u))
        if u_max < 1e-10:
            u_max = max(abs(self.c), abs(self.d))

        # Wave speeds
        speed_x = np.abs(self.a) * u_max
        speed_y = np.abs(self.b)

        # Avoid division by zero
        if speed_x < 1e-10:
            speed_x = 1e-10
        if speed_y < 1e-10:
            speed_y = 1e-10

        # CFL-based time step
        dt_x = self.dx / speed_x
        dt_y = self.dy / speed_y

        self.dt = self.cfl * min(dt_x, dt_y)

    def _compute_fvm_flux_divergence(self, u: np.ndarray) -> np.ndarray:
        """
        Compute the flux divergence using proper Finite Volume Method with Upwind scheme.

        FVM formulation for 2D Inviscid Burgers Equation:
            ∂u/∂t + ∂(a*u²/2)/∂x + ∂(b*u)/∂y = 0

        In integral form over control volume V_c:
            ∫∫ ∂u/∂t dV + ∮ H⃗·dS⃗ = 0

        Where H⃗ = E*î + F*ĵ with:
            E = a*u²/2  (x-direction flux)
            F = b*u     (y-direction flux)

        The flux at each face using first-order upwind:
            - Face velocity: u_face = (u_C + u_neighbor)/2
            - If face_velocity > 0: use upwind (upstream) cell value
            - If face_velocity < 0: use downwind (downstream) cell value

        Parameters
        ----------
        u : np.ndarray
            Current solution array

        Returns
        -------
        flux_div : np.ndarray
            Total flux divergence (F_e + F_w + F_n + F_s) / cell_area
        """
        N = self.N
        flux_div = np.zeros_like(u)

        # Cell area (for uniform grid)
        cell_area = self.dx * self.dy

        # Surface areas
        # k_e = +dy (east face, normal points +x)
        # k_w = -dy (west face, normal points -x)
        # l_n = +dx (north face, normal points +y)
        # l_s = -dx (south face, normal points -y)
        k_e = self.dy
        k_w = -self.dy
        l_n = self.dx
        l_s = -self.dx

        # Interior points (excluding boundaries which are fixed)
        for j in range(1, N - 1):  # y direction (row)
            for i in range(1, N - 1):  # x direction (column)
                u_C = u[j, i]      # Center
                u_E = u[j, i + 1]  # East
                u_W = u[j, i - 1]  # West
                u_N = u[j + 1, i]  # North
                u_S = u[j - 1, i]  # South

                # ============================================
                # East face flux: F_e = k_e * u_e * phi_e
                # ============================================
                # Face velocity (average)
                u_e = 0.5 * (u_C + u_E)

                # Upwind selection for phi_e based on k_e * u_e
                # k_e > 0, so sign depends on u_e
                if self.a * u_e >= 0:
                    # Flow goes from C to E, use upwind value = u_C
                    phi_e = u_C
                else:
                    # Flow goes from E to C, use upwind value = u_E
                    phi_e = u_E

                # Flux at east face: E_e * area = (a/2)*u_e*phi_e * dy
                # Using the formulation from PDF: k_e * u_e * phi_e where the flux is a*u*phi
                # But for Burgers: E = a*u²/2, so at face: E_e = (a/2)*u_e*u_e = (a/2)*phi_e*phi_e
                # Actually, the proper upwind flux for Burgers is:
                # F_{i+1/2} = (a/2) * phi_e² where phi_e is the upwind value
                flux_e = k_e * (self.a / 2.0) * phi_e * phi_e

                # ============================================
                # West face flux: F_w = k_w * u_w * phi_w
                # ============================================
                # Face velocity (average)
                u_w = 0.5 * (u_C + u_W)

                # Upwind selection for phi_w
                # k_w < 0, so we need to consider the flux direction
                # At west face, positive flow (u_w > 0) means flow from W to C
                if self.a * u_w >= 0:
                    # Flow goes from W to C, use upwind value = u_W
                    phi_w = u_W
                else:
                    # Flow goes from C to W, use upwind value = u_C
                    phi_w = u_C

                # Flux at west face
                flux_w = k_w * (self.a / 2.0) * phi_w * phi_w

                # ============================================
                # North face flux: F_n = l_n * phi_n
                # ============================================
                # For y-direction: F = b*u, so flux = b*phi_n * dx
                # v = b is constant (given as 1 in the problem)

                # Upwind selection for phi_n
                if self.b >= 0:
                    # Flow goes from C to N (upward), use upwind value = u_C
                    phi_n = u_C
                else:
                    # Flow goes from N to C (downward), use upwind value = u_N
                    phi_n = u_N

                # Flux at north face
                flux_n = l_n * self.b * phi_n

                # ============================================
                # South face flux: F_s = l_s * phi_s
                # ============================================
                # Upwind selection for phi_s
                if self.b >= 0:
                    # Flow goes from S to C (upward), use upwind value = u_S
                    phi_s = u_S
                else:
                    # Flow goes from C to S (downward), use upwind value = u_C
                    phi_s = u_C

                # Flux at south face
                flux_s = l_s * self.b * phi_s

                # ============================================
                # Total flux divergence = (sum of fluxes) / cell_area
                # ============================================
                flux_div[j, i] = (flux_e + flux_w + flux_n + flux_s) / cell_area

        return flux_div

    def _compute_fvm_flux_divergence_vectorized(self, u: np.ndarray) -> np.ndarray:
        """
        Vectorized computation of FVM flux divergence using upwind scheme.

        More efficient implementation using numpy operations.
        """
        N = self.N
        flux_div = np.zeros_like(u)

        # Cell area
        cell_area = self.dx * self.dy

        # Surface areas
        k_e = self.dy
        k_w = -self.dy
        l_n = self.dx
        l_s = -self.dx

        # Interior slices
        int_j = slice(1, N - 1)  # rows 1 to N-2
        int_i = slice(1, N - 1)  # cols 1 to N-2

        # Get neighbor values
        u_C = u[int_j, int_i]
        u_E = u[int_j, 2:]       # i+1
        u_W = u[int_j, :-2]      # i-1
        u_N = u[2:, int_i]       # j+1
        u_S = u[:-2, int_i]      # j-1

        # ============================================
        # East face
        # ============================================
        u_e = 0.5 * (u_C + u_E)
        phi_e = np.where(self.a * u_e >= 0, u_C, u_E)
        flux_e = k_e * (self.a / 2.0) * phi_e * phi_e

        # ============================================
        # West face
        # ============================================
        u_w = 0.5 * (u_C + u_W)
        phi_w = np.where(self.a * u_w >= 0, u_W, u_C)
        flux_w = k_w * (self.a / 2.0) * phi_w * phi_w

        # ============================================
        # North face
        # ============================================
        if self.b >= 0:
            phi_n = u_C
        else:
            phi_n = u_N
        flux_n = l_n * self.b * phi_n

        # ============================================
        # South face
        # ============================================
        if self.b >= 0:
            phi_s = u_S
        else:
            phi_s = u_C
        flux_s = l_s * self.b * phi_s

        # Total flux divergence
        flux_div[int_j, int_i] = (flux_e + flux_w + flux_n + flux_s) / cell_area

        return flux_div

    def _compute_upwind_flux_x(self, u: np.ndarray) -> np.ndarray:
        """
        Compute the convective flux in x-direction using first-order upwind.

        DEPRECATED: Use _compute_fvm_flux_divergence instead for proper FVM.
        Kept for backward compatibility.

        For the term a*u*du/dx = d(a*u^2/2)/dx:
            - If a*u > 0: use left (backward) value
            - If a*u < 0: use right (forward) value

        Parameters
        ----------
        u : np.ndarray
            Current solution array

        Returns
        -------
        flux_x : np.ndarray
            Flux contribution in x-direction

        FVM Derivation Note:
        --------------------
        The discrete equation for cell (i,j) is:
        u_new = u_old - (dt/dx)*(F_{i+1/2} - F_{i-1/2}) - ...

        Where the numerical flux F_{i+1/2} at the interface is approximated using
        the First-Order Upwind Scheme based on the wave speed 'a*u':
        - If a*u > 0: F_{i+1/2} ≈ a * u_{i}   (Information flows from left)
        - If a*u < 0: F_{i+1/2} ≈ a * u_{i+1} (Information flows from right)

        This code implements this by computing the net flux divergence directly.
        """
        N = self.N
        flux = np.zeros_like(u)

        # Interior points (excluding boundaries)
        for j in range(1, N - 1):  # y direction
            for i in range(1, N - 1):  # x direction
                # Local velocity for upwind determination
                u_local = u[j, i]
                wave_speed = self.a * u_local

                if wave_speed >= 0:
                    # Use backward difference (upwind from left)
                    flux[j, i] = self.a * u_local * (u[j, i] - u[j, i-1]) / self.dx
                else:
                    # Use forward difference (upwind from right)
                    flux[j, i] = self.a * u_local * (u[j, i+1] - u[j, i]) / self.dx

        return flux

    def _compute_upwind_flux_y(self, u: np.ndarray) -> np.ndarray:
        """
        Compute the advective flux in y-direction using first-order upwind.

        DEPRECATED: Use _compute_fvm_flux_divergence instead for proper FVM.
        Kept for backward compatibility.

        For the term b*du/dy:
            - If b > 0: use backward difference
            - If b < 0: use forward difference

        Parameters
        ----------
        u : np.ndarray
            Current solution array

        Returns
        -------
        flux_y : np.ndarray
            Flux contribution in y-direction
        """
        N = self.N
        flux = np.zeros_like(u)

        # Interior points (excluding boundaries)
        for j in range(1, N - 1):  # y direction
            for i in range(1, N - 1):  # x direction
                if self.b >= 0:
                    # Use backward difference (upwind from bottom)
                    flux[j, i] = self.b * (u[j, i] - u[j-1, i]) / self.dy
                else:
                    # Use forward difference (upwind from top)
                    flux[j, i] = self.b * (u[j+1, i] - u[j, i]) / self.dy

        return flux

    def _compute_upwind_flux_x_vectorized(self, u: np.ndarray) -> np.ndarray:
        """
        Vectorized computation of x-direction flux using upwind scheme.

        DEPRECATED: Use _compute_fvm_flux_divergence_vectorized instead for proper FVM.
        Kept for backward compatibility.

        More efficient implementation using numpy operations.
        """
        flux = np.zeros_like(u)

        # Interior points
        interior_j = slice(1, self.N - 1)
        interior_i = slice(1, self.N - 1)

        # Wave speed based on local velocity
        wave_speed = self.a * u[interior_j, interior_i]

        # Backward difference (for wave_speed >= 0)
        du_backward = (u[interior_j, interior_i] - u[interior_j, :-2]) / self.dx

        # Forward difference (for wave_speed < 0)
        du_forward = (u[interior_j, 2:] - u[interior_j, interior_i]) / self.dx

        # Upwind selection
        flux[interior_j, interior_i] = np.where(
            wave_speed >= 0,
            self.a * u[interior_j, interior_i] * du_backward,
            self.a * u[interior_j, interior_i] * du_forward
        )

        return flux

    def _compute_upwind_flux_y_vectorized(self, u: np.ndarray) -> np.ndarray:
        """
        Vectorized computation of y-direction flux using upwind scheme.

        DEPRECATED: Use _compute_fvm_flux_divergence_vectorized instead for proper FVM.
        Kept for backward compatibility.
        """
        flux = np.zeros_like(u)

        # Interior points
        interior_j = slice(1, self.N - 1)
        interior_i = slice(1, self.N - 1)

        if self.b >= 0:
            # Backward difference (upwind from bottom)
            flux[interior_j, interior_i] = self.b * (
                u[interior_j, interior_i] - u[:-2, interior_i]
            ) / self.dy
        else:
            # Forward difference (upwind from top)
            flux[interior_j, interior_i] = self.b * (
                u[2:, interior_i] - u[interior_j, interior_i]
            ) / self.dy

        return flux

    def _solve_tdma(self, a: np.ndarray, b: np.ndarray, c: np.ndarray,
                    d: np.ndarray) -> np.ndarray:
        """
        Solve a tri-diagonal system using the Thomas Algorithm (TDMA).

        The system is:
            a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]

        Algorithm:
            Forward sweep: Eliminate lower diagonal
            Backward substitution: Solve for x

        Parameters
        ----------
        a : np.ndarray
            Lower diagonal (a[0] = 0)
        b : np.ndarray
            Main diagonal
        c : np.ndarray
            Upper diagonal (c[n-1] = 0)
        d : np.ndarray
            Right-hand side

        Returns
        -------
        x : np.ndarray
            Solution vector
        """
        n = len(d)
        c_prime = np.zeros(n)
        d_prime = np.zeros(n)
        x = np.zeros(n)

        # Forward sweep
        c_prime[0] = c[0] / b[0] if abs(b[0]) > 1e-15 else 0
        d_prime[0] = d[0] / b[0] if abs(b[0]) > 1e-15 else 0

        for i in range(1, n):
            denom = b[i] - a[i] * c_prime[i-1]
            if abs(denom) < 1e-15:
                denom = 1e-15
            c_prime[i] = c[i] / denom
            d_prime[i] = (d[i] - a[i] * d_prime[i-1]) / denom

        # Backward substitution
        x[n-1] = d_prime[n-1]
        for i in range(n-2, -1, -1):
            x[i] = d_prime[i] - c_prime[i] * x[i+1]

        return x

    def time_step_implicit(self) -> float:
        """
        Perform one time step using fully implicit method with direct TDMA.

        Equation: ∂u/∂t + ∂(a*u²/2)/∂x + ∂(b*u)/∂y = 0

        For implicit method: (u^(n+1) - u^n)/dt + flux_div(u^(n+1)) = 0
        Linearize nonlinear x-flux: d(a*u²/2)/dx ≈ a*u^n * d(u^(n+1))/dx

        Solve each row and column iteratively using TDMA.

        Returns
        -------
        residual : float
            L2 norm of the change in solution
        """
        u_old = self.u.copy()
        N = self.N

        # Compute time step
        self._compute_time_step()
        dt = self.dt

        # For better accuracy, do a few sub-iterations (Gauss-Seidel style)
        n_sweeps = 5  # Number of ADI sweeps per time step

        u_new = u_old.copy()

        for sweep in range(n_sweeps):
            # Step 1: Implicit solve in x-direction (line-by-line for each row)
            for j in range(1, N-1):  # For each interior row
                # Build tri-diagonal system for this row
                aa = np.zeros(N)
                bb = np.zeros(N)
                cc = np.zeros(N)
                dd = np.zeros(N)

                for i in range(N):
                    if i == 0:  # Left boundary
                        bb[i] = 1.0
                        dd[i] = self.c
                    elif i == N-1:  # Right boundary
                        bb[i] = 1.0
                        dd[i] = self.d
                    else:  # Interior points
                        # Conservative form: d(a*u²/2)/dx
                        # Linearize: d(a*u²/2)/dx ≈ a*u^n * du^(n+1)/dx
                        u_ref = u_new[j, i]  # Use latest available value
                        coeff_x = self.a * u_ref

                        # Upwind based on wave speed
                        if coeff_x >= 0:
                            # Backward difference
                            alpha_x = dt * coeff_x / self.dx
                            aa[i] = -alpha_x
                            bb[i] = 1.0 + alpha_x
                            cc[i] = 0.0
                        else:
                            # Forward difference
                            alpha_x = -dt * coeff_x / self.dx
                            aa[i] = 0.0
                            bb[i] = 1.0 + alpha_x
                            cc[i] = -alpha_x

                        # RHS includes old value and y-direction flux (explicit for now)
                        flux_y = 0.0
                        if self.b >= 0 and j > 0:
                            flux_y = dt * self.b * (u_new[j, i] - u_new[j-1, i]) / self.dy
                        elif self.b < 0 and j < N-1:
                            flux_y = dt * self.b * (u_new[j+1, i] - u_new[j, i]) / self.dy

                        dd[i] = u_old[j, i] - flux_y

                # Solve TDMA for this row
                u_new[j, :] = self._solve_tdma(aa, bb, cc, dd)

            # Step 2: Implicit solve in y-direction (line-by-line for each column)
            for i in range(1, N-1):  # For each interior column
                # Build tri-diagonal system for this column
                aa = np.zeros(N)
                bb = np.zeros(N)
                cc = np.zeros(N)
                dd = np.zeros(N)

                for j in range(N):
                    if j == 0:  # Bottom boundary
                        bb[j] = 1.0
                        dd[j] = self.c - (self.c - self.d) * self.x[i]
                    elif j == N-1:  # Top boundary (Neumann)
                        bb[j] = 1.0
                        cc[j] = -1.0
                        dd[j] = 0.0
                    else:  # Interior points
                        # Y-direction: LINEAR advection b*du/dy
                        beta_y = dt * self.b / self.dy if self.b >= 0 else -dt * self.b / self.dy

                        if self.b >= 0:
                            aa[j] = -beta_y
                            bb[j] = 1.0 + beta_y
                            cc[j] = 0.0
                        else:
                            aa[j] = 0.0
                            bb[j] = 1.0 + beta_y
                            cc[j] = -beta_y

                        # RHS includes x-direction update
                        flux_x = 0.0
                        u_ref = u_new[j, i]
                        coeff_x = self.a * u_ref

                        if coeff_x >= 0 and i > 0:
                            flux_x = dt * coeff_x * (u_new[j, i] - u_new[j, i-1]) / self.dx
                        elif coeff_x < 0 and i < N-1:
                            flux_x = dt * coeff_x * (u_new[j, i+1] - u_new[j, i]) / self.dx

                        dd[j] = u_old[j, i] - flux_x

                # Solve TDMA for this column
                u_new[:, i] = self._solve_tdma(aa, bb, cc, dd)

        self.u = u_new

        # Apply boundary conditions
        self._apply_boundary_conditions()

        # Compute residual
        diff = self.u - u_old
        residual = np.sqrt(np.sum(diff**2)) / self.N

        return residual

    def time_step(self, use_vectorized: bool = True, use_fvm: bool = True,
                  use_implicit: bool = False) -> float:
        """
        Perform one time step using explicit Euler or implicit TDMA method.

        Explicit: u^{n+1} = u^n - dt * (flux_divergence)
        Implicit: Uses ADI with TDMA solver

        Parameters
        ----------
        use_vectorized : bool
            Use vectorized flux computation (faster) - for explicit only
        use_fvm : bool
            Use proper FVM flux divergence (True) or old method (False) - for explicit only
        use_implicit : bool
            Use implicit TDMA method (True) or explicit Euler (False)

        Returns
        -------
        residual : float
            L2 norm of the change in solution
        """
        if use_implicit:
            return self.time_step_implicit()

        # Original explicit method
        # Store old solution
        u_old = self.u.copy()

        # Compute time step based on current CFL
        self._compute_time_step()

        # Compute flux divergence
        if use_fvm:
            # Proper FVM with upwind scheme
            if use_vectorized:
                flux_div = self._compute_fvm_flux_divergence_vectorized(self.u)
            else:
                flux_div = self._compute_fvm_flux_divergence(self.u)
        else:
            # Old method (kept for comparison)
            if use_vectorized:
                flux_x = self._compute_upwind_flux_x_vectorized(self.u)
                flux_y = self._compute_upwind_flux_y_vectorized(self.u)
            else:
                flux_x = self._compute_upwind_flux_x(self.u)
                flux_y = self._compute_upwind_flux_y(self.u)
            flux_div = flux_x + flux_y

        # Update solution (explicit Euler)
        self.u = u_old - self.dt * flux_div

        # Apply boundary conditions
        self._apply_boundary_conditions()

        # Compute residual (L2 norm of change)
        diff = self.u - u_old
        residual = np.sqrt(np.sum(diff**2)) / self.N

        return residual

    def solve_steady_state(self, max_iterations: int = 10000,
                           tolerance: float = 1e-6,
                           print_interval: int = 500,
                           use_vectorized: bool = True,
                           use_fvm: bool = True,
                           use_implicit: bool = False) -> Dict:
        """
        Solve until steady state is reached.

        Iterates until the residual falls below the tolerance or
        maximum iterations is reached.

        Parameters
        ----------
        max_iterations : int
            Maximum number of time steps
        tolerance : float
            Convergence tolerance for residual
        print_interval : int
            Print progress every this many iterations
        use_vectorized : bool
            Use vectorized flux computation (explicit only)
        use_fvm : bool
            Use proper FVM flux divergence (explicit only)
        use_implicit : bool
            Use implicit TDMA method (ADI) instead of explicit Euler

        Returns
        -------
        result : dict
            Dictionary containing:
            - 'converged': bool
            - 'iterations': int
            - 'final_residual': float
            - 'residual_history': list
            - 'time_elapsed': float
            - 'method': str
        """
        method_name = "Implicit TDMA (ADI)" if use_implicit else "Explicit Euler"
        print("\n" + "="*60)
        print("Starting steady-state solver...")
        print(f"  Method: {method_name}")
        print(f"  Max iterations: {max_iterations}")
        print(f"  Tolerance: {tolerance:.2e}")
        if not use_implicit:
            print(f"  Using FVM: {use_fvm}")
        print("="*60 + "\n")

        start_time = time.time()
        self.residual_history = []
        converged = False

        for iteration in range(1, max_iterations + 1):
            residual = self.time_step(
                use_vectorized=use_vectorized,
                use_fvm=use_fvm,
                use_implicit=use_implicit
            )
            self.residual_history.append(residual)
            self.iteration_count = iteration

            # Print progress
            if iteration % print_interval == 0 or iteration == 1:
                print(f"  Iteration {iteration:6d}: residual = {residual:.6e}, "
                      f"dt = {self.dt:.6e}")

            # Check convergence
            if residual < tolerance:
                converged = True
                print(f"\n  CONVERGED at iteration {iteration}!")
                print(f"  Final residual: {residual:.6e}")
                break

        elapsed_time = time.time() - start_time

        if not converged:
            print(f"\n  WARNING: Did not converge within {max_iterations} iterations")
            print(f"  Final residual: {self.residual_history[-1]:.6e}")

        print(f"\n  Time elapsed: {elapsed_time:.2f} seconds")
        print("="*60 + "\n")

        return {
            'converged': converged,
            'iterations': self.iteration_count,
            'final_residual': self.residual_history[-1] if self.residual_history else 0,
            'residual_history': self.residual_history.copy(),
            'time_elapsed': elapsed_time,
            'method': method_name
        }

    def get_solution(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Get the current solution.

        Returns
        -------
        X : np.ndarray
            Meshgrid of x-coordinates
        Y : np.ndarray
            Meshgrid of y-coordinates
        u : np.ndarray
            Solution array
        """
        return self.X, self.Y, self.u.copy()

    def compute_error(self, analytical_solution: callable) -> Dict:
        """
        Compute error metrics against an analytical solution.

        Parameters
        ----------
        analytical_solution : callable
            Function that takes (x, y) and returns analytical u

        Returns
        -------
        errors : dict
            Dictionary containing:
            - 'L1': L1 norm error
            - 'L2': L2 norm error
            - 'Linf': Maximum absolute error
            - 'relative_L2': Relative L2 error
        """
        u_analytical = analytical_solution(self.X, self.Y)
        diff = self.u - u_analytical

        L1 = np.mean(np.abs(diff))
        L2 = np.sqrt(np.mean(diff**2))
        Linf = np.max(np.abs(diff))

        # Relative L2 error
        u_analytical_norm = np.sqrt(np.mean(u_analytical**2))
        if u_analytical_norm > 1e-10:
            relative_L2 = L2 / u_analytical_norm
        else:
            relative_L2 = L2

        return {
            'L1': L1,
            'L2': L2,
            'Linf': Linf,
            'relative_L2': relative_L2
        }

    def plot_solution(self, save_path: str = None,
                      show_plot: bool = True,
                      title_suffix: str = "") -> plt.Figure:
        """
        Plot the numerical solution.

        Parameters
        ----------
        save_path : str, optional
            Path to save the figure
        show_plot : bool
            Whether to display the plot
        title_suffix : str
            Additional text to add to plot titles

        Returns
        -------
        fig : matplotlib.figure.Figure
            The generated figure
        """
        fig = plt.figure(figsize=(16, 5))

        # Subplot 1: Filled contour
        ax1 = fig.add_subplot(131)
        contourf = ax1.contourf(self.X, self.Y, self.u, levels=20, cmap='RdYlBu_r')
        ax1.set_xlabel('x', fontsize=12)
        ax1.set_ylabel('y', fontsize=12)
        ax1.set_title(f'Numerical Solution (N={self.N})\n{title_suffix}', fontsize=11)
        ax1.set_aspect('equal')
        ax1.grid(True, alpha=0.3)
        plt.colorbar(contourf, ax=ax1, label='u(x, y)')

        # Subplot 2: Contour lines
        ax2 = fig.add_subplot(132)
        contour = ax2.contour(self.X, self.Y, self.u, levels=15, colors='black', linewidths=0.5)
        ax2.clabel(contour, inline=True, fontsize=8, fmt='%.2f')
        ax2.contourf(self.X, self.Y, self.u, levels=20, cmap='RdYlBu_r', alpha=0.6)
        ax2.set_xlabel('x', fontsize=12)
        ax2.set_ylabel('y', fontsize=12)
        ax2.set_title('Contour Lines', fontsize=11)
        ax2.set_aspect('equal')
        ax2.grid(True, alpha=0.3)

        # Subplot 3: 3D surface
        ax3 = fig.add_subplot(133, projection='3d')
        surf = ax3.plot_surface(self.X, self.Y, self.u, cmap='RdYlBu_r',
                                 edgecolor='none', alpha=0.9)
        ax3.set_xlabel('x', fontsize=10)
        ax3.set_ylabel('y', fontsize=10)
        ax3.set_zlabel('u(x, y)', fontsize=10)
        ax3.set_title('3D Surface', fontsize=11)
        ax3.view_init(elev=25, azim=45)
        plt.colorbar(surf, ax=ax3, shrink=0.5, aspect=5)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved to: {save_path}")

        if show_plot:
            plt.show()

        return fig

    def plot_convergence(self, save_path: str = None,
                         show_plot: bool = True) -> plt.Figure:
        """
        Plot the convergence history.

        Parameters
        ----------
        save_path : str, optional
            Path to save the figure
        show_plot : bool
            Whether to display the plot

        Returns
        -------
        fig : matplotlib.figure.Figure
            The generated figure
        """
        if not self.residual_history:
            print("No convergence history available. Run solve_steady_state first.")
            return None

        fig, ax = plt.subplots(figsize=(10, 6))

        iterations = range(1, len(self.residual_history) + 1)
        ax.semilogy(iterations, self.residual_history, 'b-', linewidth=1.5)
        ax.set_xlabel('Iteration', fontsize=12)
        ax.set_ylabel('Residual (L2 norm)', fontsize=12)
        ax.set_title(f'Convergence History (N={self.N})\n'
                     f'a={self.a}, b={self.b}, c={self.c}, d={self.d}', fontsize=12)
        ax.grid(True, alpha=0.3)

        # Mark final value
        final_iter = len(self.residual_history)
        final_res = self.residual_history[-1]
        ax.axhline(y=final_res, color='r', linestyle='--', alpha=0.5,
                   label=f'Final: {final_res:.2e}')
        ax.legend()

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved to: {save_path}")

        if show_plot:
            plt.show()

        return fig

    def plot_comparison(self, analytical_solution: callable,
                        save_path: str = None,
                        show_plot: bool = True) -> plt.Figure:
        """
        Plot comparison between numerical and analytical solutions.

        Parameters
        ----------
        analytical_solution : callable
            Function that takes (x, y) and returns analytical u
        save_path : str, optional
            Path to save the figure
        show_plot : bool
            Whether to display the plot

        Returns
        -------
        fig : matplotlib.figure.Figure
            The generated figure
        """
        u_analytical = analytical_solution(self.X, self.Y)
        error = self.u - u_analytical

        fig = plt.figure(figsize=(16, 10))

        # Row 1: Solutions
        # Numerical
        ax1 = fig.add_subplot(231)
        cf1 = ax1.contourf(self.X, self.Y, self.u, levels=20, cmap='RdYlBu_r')
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_title(f'Numerical Solution (N={self.N})')
        ax1.set_aspect('equal')
        plt.colorbar(cf1, ax=ax1)

        # Analytical
        ax2 = fig.add_subplot(232)
        cf2 = ax2.contourf(self.X, self.Y, u_analytical, levels=20, cmap='RdYlBu_r')
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        ax2.set_title('Analytical Solution')
        ax2.set_aspect('equal')
        plt.colorbar(cf2, ax=ax2)

        # Error
        ax3 = fig.add_subplot(233)
        cf3 = ax3.contourf(self.X, self.Y, error, levels=20, cmap='coolwarm')
        ax3.set_xlabel('x')
        ax3.set_ylabel('y')
        ax3.set_title(f'Error (Numerical - Analytical)')
        ax3.set_aspect('equal')
        plt.colorbar(cf3, ax=ax3)

        # Row 2: Line profiles
        # Profile at y = 0.25
        y_idx_025 = np.argmin(np.abs(self.y - 0.25))
        ax4 = fig.add_subplot(234)
        ax4.plot(self.x, self.u[y_idx_025, :], 'b-', linewidth=2, label='Numerical')
        ax4.plot(self.x, u_analytical[y_idx_025, :], 'r--', linewidth=2, label='Analytical')
        ax4.set_xlabel('x')
        ax4.set_ylabel('u')
        ax4.set_title(f'Profile at y = {self.y[y_idx_025]:.2f}')
        ax4.legend()
        ax4.grid(True, alpha=0.3)

        # Profile at y = 0.5
        y_idx_05 = np.argmin(np.abs(self.y - 0.5))
        ax5 = fig.add_subplot(235)
        ax5.plot(self.x, self.u[y_idx_05, :], 'b-', linewidth=2, label='Numerical')
        ax5.plot(self.x, u_analytical[y_idx_05, :], 'r--', linewidth=2, label='Analytical')
        ax5.set_xlabel('x')
        ax5.set_ylabel('u')
        ax5.set_title(f'Profile at y = {self.y[y_idx_05]:.2f}')
        ax5.legend()
        ax5.grid(True, alpha=0.3)

        # Profile at y = 0.75
        y_idx_075 = np.argmin(np.abs(self.y - 0.75))
        ax6 = fig.add_subplot(236)
        ax6.plot(self.x, self.u[y_idx_075, :], 'b-', linewidth=2, label='Numerical')
        ax6.plot(self.x, u_analytical[y_idx_075, :], 'r--', linewidth=2, label='Analytical')
        ax6.set_xlabel('x')
        ax6.set_ylabel('u')
        ax6.set_title(f'Profile at y = {self.y[y_idx_075]:.2f}')
        ax6.legend()
        ax6.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved to: {save_path}")

        if show_plot:
            plt.show()

        return fig

    def plot_grid(self, save_path: str = None, show_plot: bool = True,
                  title: str = "Grid Generation") -> plt.Figure:
        """
        Plot the computational grid (Mesh).
        
        Parameters
        ----------
        save_path : str, optional
            Path to save the figure
        show_plot : bool
            Whether to display the plot
        title : str
            Title of the plot
            
        Returns
        -------
        fig : matplotlib.figure.Figure
            The generated figure
        """
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        
        # Plot vertical lines (constant x)
        for i in range(self.N):
            ax.plot(self.X[:, i], self.Y[:, i], 'r-', linewidth=0.5, alpha=0.7)
            
        # Plot horizontal lines (constant y)
        for j in range(self.N):
            ax.plot(self.X[j, :], self.Y[j, :], 'r-', linewidth=0.5, alpha=0.7)
            
        ax.set_xlabel('x (Physical) / $\\xi$ (Computational)', fontsize=12)
        ax.set_ylabel('y (Physical) / $\\eta$ (Computational)', fontsize=12)
        ax.set_title(f'{title}\nGrid size: {self.N}x{self.N}', fontsize=14)
        ax.set_aspect('equal')
        
        # Set limits with small padding
        padding = 0.02
        ax.set_xlim(self.x_min - padding, self.x_max + padding)
        ax.set_ylim(self.y_min - padding, self.y_max + padding)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved to: {save_path}")
            
        if show_plot:
            plt.show()
            
        return fig

    def print_summary(self):
        """Print a summary of the solver state and solution."""
        print("\n" + "="*60)
        print("Burgers Solver Summary")
        print("="*60)
        print(f"Grid: {self.N} x {self.N}")
        print(f"Parameters: a={self.a}, b={self.b}, c={self.c}, d={self.d}")
        print(f"Grid spacing: dx={self.dx:.6f}, dy={self.dy:.6f}")
        print(f"Current dt: {self.dt:.6e}")
        print(f"CFL number: {self.cfl}")
        print("-"*60)
        print("Solution statistics:")
        print(f"  Min: {self.u.min():.6f}")
        print(f"  Max: {self.u.max():.6f}")
        print(f"  Mean: {self.u.mean():.6f}")
        print(f"  Std: {self.u.std():.6f}")
        if self.residual_history:
            print("-"*60)
            print("Convergence info:")
            print(f"  Iterations: {self.iteration_count}")
            print(f"  Final residual: {self.residual_history[-1]:.6e}")
        print("="*60 + "\n")


def run_case1_demo(N: int = 41, max_iterations: int = 10000,
                   tolerance: float = 1e-6) -> BurgersSolver2D:
    """
    Run demonstration for Case 1.

    Parameters
    ----------
    N : int
        Grid size
    max_iterations : int
        Maximum iterations for steady state solver
    tolerance : float
        Convergence tolerance

    Returns
    -------
    solver : BurgersSolver2D
        The configured and solved solver instance
    """
    print("\n" + "#"*70)
    print(" Case 1: a=1, b=1, c=1.5, d=-0.5")
    print("#"*70)

    # Create solver
    solver = BurgersSolver2D(N=N, a=1.0, b=1.0, c=1.5, d=-0.5, cfl=0.5)

    # Solve to steady state
    result = solver.solve_steady_state(
        max_iterations=max_iterations,
        tolerance=tolerance,
        print_interval=1000
    )

    # Print summary
    solver.print_summary()

    return solver


if __name__ == "__main__":
    # Run demo
    print("\nTesting BurgersSolver2D...")

    # Test with N=21
    solver_21 = run_case1_demo(N=21, max_iterations=5000, tolerance=1e-6)

    # Test with N=41
    solver_41 = run_case1_demo(N=41, max_iterations=10000, tolerance=1e-6)

    print("\nTest completed successfully!")
