"""
Finite Volume Method Solver for 2D Inviscid Burgers Equation with TVD Schemes

Implements:
    - First-Order Upwind Scheme
    - TVD Schemes with various flux limiters
    - Explicit Euler and Implicit TDMA time integration
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from typing import Tuple, Optional, Dict, List, Callable, Union
import time
from enum import Enum


class FluxLimiter(Enum):
    """Enumeration of available flux limiters for TVD schemes."""
    UPWIND = "upwind"
    MINMOD = "minmod"
    VAN_ALBADA = "van_albada"
    UMIST = "umist"
    SUPERBEE = "superbee"
    VAN_LEER = "van_leer"
    KOREN = "koren"
    QUICK = "quick"
    SMART = "smart"
    HCUSP = "hcusp"
    HQUICK = "hquick"


class SolverMethod(Enum):
    """Enumeration of available solver methods."""
    EXPLICIT_EULER = "explicit"
    IMPLICIT_TDMA = "implicit_tdma"


# =================================================================================
# FLUX LIMITER FUNCTIONS
# =================================================================================

def psi_upwind(r: np.ndarray) -> np.ndarray:
    """Upwind limiter: ψ(r) = 0"""
    return np.zeros_like(r)


def psi_minmod(r: np.ndarray) -> np.ndarray:
    """Min-Mod limiter"""
    return np.maximum(0, np.minimum(r, 1))


def psi_van_albada(r: np.ndarray) -> np.ndarray:
    """Van Albada limiter"""
    return (r + r**2) / (1 + r**2)


def psi_umist(r: np.ndarray) -> np.ndarray:
    """UMIST limiter"""
    return np.maximum(0, np.minimum(
        np.minimum(2*r, (1 + 3*r)/4),
        np.minimum((3 + r)/4, 2)
    ))


def psi_superbee(r: np.ndarray) -> np.ndarray:
    """Superbee limiter"""
    return np.maximum(0, np.maximum(
        np.minimum(2*r, 1),
        np.minimum(r, 2)
    ))


def psi_van_leer(r: np.ndarray) -> np.ndarray:
    """Van Leer limiter"""
    return (r + np.abs(r)) / (1 + np.abs(r))


def psi_koren(r: np.ndarray) -> np.ndarray:
    """Koren limiter"""
    return np.maximum(0, np.minimum(
        np.minimum(2*r, (2 + r)/3),
        2
    ))


def psi_quick(r: np.ndarray) -> np.ndarray:
    """TVD-limited QUICK"""
    return np.maximum(0, np.minimum(
        np.minimum(2*r, (3 + r)/4),
        2
    ))


def psi_smart(r: np.ndarray) -> np.ndarray:
    """SMART limiter"""
    return np.maximum(0, np.minimum(
        np.minimum(2*r, (1 + 3*r)/4),
        4
    ))


def psi_hcusp(r: np.ndarray) -> np.ndarray:
    """H-CUSP limiter"""
    with np.errstate(divide='ignore', invalid='ignore'):
        result = 1.5 * (r + np.abs(r)) / (r + 2)
        result = np.where(np.isfinite(result), result, 0)
    return np.maximum(0, result)


def psi_hquick(r: np.ndarray) -> np.ndarray:
    """H-QUICK limiter"""
    with np.errstate(divide='ignore', invalid='ignore'):
        result = 2 * (r + np.abs(r)) / (r + 3)
        result = np.where(np.isfinite(result), result, 0)
    return np.maximum(0, result)


def get_limiter_function(limiter: FluxLimiter) -> Callable:
    """
    Get the flux limiter function for the given limiter type.

    Parameters
    ----------
    limiter : FluxLimiter
        The type of flux limiter to use

    Returns
    -------
    psi : Callable
        The flux limiter function ψ(r)
    """
    limiter_functions = {
        FluxLimiter.UPWIND: psi_upwind,
        FluxLimiter.MINMOD: psi_minmod,
        FluxLimiter.VAN_ALBADA: psi_van_albada,
        FluxLimiter.UMIST: psi_umist,
        FluxLimiter.SUPERBEE: psi_superbee,
        FluxLimiter.VAN_LEER: psi_van_leer,
        FluxLimiter.KOREN: psi_koren,
        FluxLimiter.QUICK: psi_quick,
        FluxLimiter.SMART: psi_smart,
        FluxLimiter.HCUSP: psi_hcusp,
        FluxLimiter.HQUICK: psi_hquick,
    }
    return limiter_functions[limiter]


# =================================================================================
# TDMA SOLVER
# =================================================================================

def solve_tdma(a: np.ndarray, b: np.ndarray, c: np.ndarray,
               d: np.ndarray) -> np.ndarray:
    """Solve a tri-diagonal system using the Thomas Algorithm (TDMA)."""
    n = len(d)

    # Create working arrays
    c_prime = np.zeros(n)
    d_prime = np.zeros(n)
    x = np.zeros(n)

    # Forward sweep
    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]

    for i in range(1, n):
        denom = b[i] - a[i] * c_prime[i-1]
        if abs(denom) < 1e-15:
            denom = 1e-15  # Prevent division by zero
        c_prime[i] = c[i] / denom
        d_prime[i] = (d[i] - a[i] * d_prime[i-1]) / denom

    # Backward substitution
    x[n-1] = d_prime[n-1]
    for i in range(n-2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i+1]

    return x


# =================================================================================
# MAIN SOLVER CLASS
# =================================================================================

class BurgersSolverFVM_TVD:
    """
    2D Inviscid Burgers Equation Solver using Finite Volume Method
    with Upwind and TVD Schemes.

    ==========================================================================
    MATHEMATICAL DERIVATION - FINITE VOLUME METHOD
    ==========================================================================

    GOVERNING EQUATION (Conservative Form):
    ---------------------------------------
    Starting from the inviscid Burgers equation:

        ∂u/∂t + a*u*∂u/∂x + b*∂u/∂y = 0

    Rewriting in conservative form:

        ∂u/∂t + ∂E/∂x + ∂F/∂y = 0

    Where:
        E = a*u²/2  (flux in x-direction)
        F = b*u     (flux in y-direction)

    Or in vector notation:
        ∂u/∂t + ∇·H = 0,  where H = (E, F)

    FVM DISCRETIZATION:
    -------------------
    Integrating over a control volume V_c:

        ∫∫_{V_c} ∂u/∂t dV + ∫∫_{V_c} ∇·H dV = 0

    Using the divergence theorem for the second term:

        ∂/∂t(∫∫_{V_c} u dV) + ∮_{∂V_c} H·n dS = 0

    For a cell-centered finite volume with area ΔS = Δx·Δy:

        ΔS * (u_C^{n+1} - u_C^n)/Δt + ∑_{faces} H_f · S_f = 0

    Where the sum is over all faces (East, West, North, South).

    FACE FLUXES:
    ------------
    For a structured rectangular grid:

    East face (at i+1/2):   Flux_e = E_e · Δy = (a*u_e²/2) · Δy
    West face (at i-1/2):   Flux_w = -E_w · Δy = -(a*u_w²/2) · Δy
    North face (at j+1/2):  Flux_n = F_n · Δx = (b*u_n) · Δx
    South face (at j-1/2):  Flux_s = -F_s · Δx = -(b*u_s) · Δx

    (Note: Signs follow outward normal convention)

    UPWIND SCHEME:
    --------------
    For the face value approximation, the upwind scheme uses:

    If wave speed > 0 (information flows from left to right):
        u_e = u_C  (upwind value from center)
        u_w = u_W  (upwind value from west)

    If wave speed < 0 (information flows from right to left):
        u_e = u_E  (upwind value from east)
        u_w = u_C  (upwind value from center)

    For the x-direction (nonlinear term a*u*∂u/∂x):
        Wave speed = a * u_face ≈ a * (u_C + u_neighbor)/2

    For the y-direction (linear term b*∂u/∂y):
        Wave speed = b (constant)

    TVD SCHEME WITH FLUX LIMITERS:
    ------------------------------
    The TVD scheme interpolates face values using:

        u_f = u_C + (1/2)*ψ(r)*(u_D - u_C)

    Where:
        u_C = cell center value (upwind cell)
        u_D = downwind cell value
        r = ratio of consecutive gradients
        ψ(r) = flux limiter function

    For flow in +x direction at east face:
        r_e⁺ = (u_C² - u_W²)/(u_E² - u_C²)   [for conservative form]
        or
        r_e⁺ = (u_C - u_W)/(u_E - u_C)       [for non-conservative form]

    The flux limiter ψ(r) must satisfy:
        - ψ(1) = 1 (second-order for uniform flow)
        - 0 ≤ ψ(r) ≤ min(2r, 2) (TVD constraint)
        - ψ(r) = 0 for r ≤ 0 (upwind at extrema)

    DISCRETE EQUATION (Explicit):
    -----------------------------
    Rearranging:

        u_C^{n+1} = u_C^n - (Δt/ΔS) * [Flux_e - Flux_w + Flux_n - Flux_s]

    With appropriate upwind/TVD face value interpolations.

    DISCRETE EQUATION (Implicit with TDMA):
    ---------------------------------------
    For implicit time integration:

        a_C*u_C^{n+1} + a_E*u_E^{n+1} + a_W*u_W^{n+1} = b

    Where b contains contributions from north/south fluxes (treated explicitly)
    and the time derivative term.

    The coefficients are derived from linearizing the flux terms.

    ==========================================================================

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
    limiter : FluxLimiter, optional
        Flux limiter type for TVD scheme (default: UPWIND)
    """

    def __init__(self, N: int, a: float, b: float, c: float, d: float,
                 cfl: float = 0.5, limiter: FluxLimiter = FluxLimiter.UPWIND):
        """Initialize the Burgers equation solver."""

        # Store parameters
        self.N = N
        self.a = a
        self.b_coeff = b  # Using b_coeff to avoid conflict with array 'b' later
        self.c = c
        self.d = d
        self.cfl = cfl
        self.limiter = limiter
        self.psi = get_limiter_function(limiter)

        # Domain bounds
        self.x_min, self.x_max = 0.0, 1.0
        self.y_min, self.y_max = 0.0, 1.0

        # Grid spacing (for N points, there are N-1 intervals)
        self.dx = (self.x_max - self.x_min) / (N - 1)
        self.dy = (self.y_max - self.y_min) / (N - 1)

        # Cell area
        self.dS = self.dx * self.dy

        # Create grid (node-centered)
        self.x = np.linspace(self.x_min, self.x_max, N)
        self.y = np.linspace(self.y_min, self.y_max, N)
        self.X, self.Y = np.meshgrid(self.x, self.y)

        # Initialize solution array [j, i] indexing where j=row(y), i=col(x)
        self.u = np.zeros((N, N))

        # Time step
        self.dt = 0.0

        # Convergence history
        self.residual_history = []
        self.iteration_count = 0

        # Initialize with boundary and initial conditions
        self._apply_initial_conditions()
        self._apply_boundary_conditions()
        self._compute_time_step()

        print(f"BurgersSolverFVM_TVD initialized:")
        print(f"  Grid: {N} x {N}")
        print(f"  Parameters: a={a}, b={b}, c={c}, d={d}")
        print(f"  Grid spacing: dx={self.dx:.6f}, dy={self.dy:.6f}")
        print(f"  Flux limiter: {limiter.value}")
        print(f"  CFL number: {cfl}")
        print(f"  Initial dt: {self.dt:.6e}")

    def __init__(self, N: int, a: float, b: float, c: float, d: float,
                 cfl: float = 0.5,
                 limiter: FluxLimiter = FluxLimiter.VAN_ALBADA):
        self.N = N
        self.a = a
        self.b_coeff = b
        self.c = c
        self.d = d
        self.cfl = cfl
        self.limiter = limiter

        self.psi = get_limiter_function(limiter)

        self.x_min, self.x_max = 0.0, 1.0
        self.y_min, self.y_max = 0.0, 1.0

        self.dx = (self.x_max - self.x_min) / (N - 1)
        self.dy = (self.y_max - self.y_min) / (N - 1)

        self.dS = self.dx * self.dy

        self.x = np.linspace(self.x_min, self.x_max, N)
        self.y = np.linspace(self.y_min, self.y_max, N)
        self.X, self.Y = np.meshgrid(self.x, self.y)

        self.u = np.zeros((N, N))
        self.dt = 0.0

        self.residual_history = []
        self.iteration_count = 0

        self._apply_initial_conditions()
        self._apply_boundary_conditions()
        self._compute_time_step()

    def _apply_initial_conditions(self):
        """Apply initial conditions: u(x, 0) = c - (c - d) * x"""
        for i in range(self.N):
            for j in range(self.N):
                x = self.x[i]
                self.u[j, i] = self.c - (self.c - self.d) * x

    def _apply_boundary_conditions(self):
        """
        Apply boundary conditions:
            Left: u(0, y) = c
            Right: u(1, y) = d
            Bottom: u(x, 0) = c - (c-d)*x
            Top: du/dy = 0 (Neumann)
        """
        # Dirichlet BCs
        self.u[:, 0] = self.c
        self.u[:, -1] = self.d

        for i in range(self.N):
            self.u[0, i] = self.c - (self.c - self.d) * self.x[i]

    def _compute_time_step(self):
        """Compute time step based on CFL condition."""
        u_max = np.max(np.abs(self.u))
        if u_max < 1e-10:
            u_max = max(abs(self.c), abs(self.d))

        speed_x = np.abs(self.a) * u_max
        speed_y = np.abs(self.b_coeff)

        if speed_x < 1e-10: speed_x = 1e-10
        if speed_y < 1e-10: speed_y = 1e-10

        dt_x = self.dx / speed_x
        dt_y = self.dy / speed_y

        self.dt = self.cfl * min(dt_x, dt_y)

    def _compute_gradient_ratio_x(self, u: np.ndarray, i: int, j: int,
                                   positive_flow: bool) -> Tuple[float, float]:
        """
        Compute gradient ratios for TVD scheme in x-direction.

        TVD GRADIENT RATIO COMPUTATION:
        -------------------------------
        For face interpolation using TVD, we need the ratio r of consecutive
        gradients to determine the limiter value.

        For positive flow direction (a*u > 0) at east face:
            Upwind cell: C = (i,j)
            Downwind cell: D = (i+1,j)
            Far upwind cell: U = (i-1,j)

            r_e⁺ = (u_C - u_U)/(u_D - u_C)
                 = (u[j,i] - u[j,i-1])/(u[j,i+1] - u[j,i])

        For negative flow direction (a*u < 0) at east face:
            Upwind cell: C = (i+1,j)
            Downwind cell: D = (i,j)
            Far upwind cell: U = (i+2,j)

            r_e⁻ = (u_C - u_U)/(u_D - u_C)
                 = (u[j,i+1] - u[j,i+2])/(u[j,i] - u[j,i+1])

        Returns
        -------
        r_e, r_w : float
            Gradient ratios for east and west face interpolations
        """
        N = self.N
        eps = 1e-10  # Small number to avoid division by zero

        if positive_flow:
            # East face (i+1/2): flow from left
            # r_e⁺ = (u_C - u_W)/(u_E - u_C)
            if i >= 1 and i < N - 1:
                delta_d = u[j, i+1] - u[j, i]
                delta_u = u[j, i] - u[j, i-1]
                r_e = delta_u / (delta_d + np.sign(delta_d + eps) * eps)
            else:
                r_e = 0

            # West face (i-1/2): flow from left
            # r_w⁺ = (u_W - u_WW)/(u_C - u_W)
            if i >= 2:
                delta_d = u[j, i] - u[j, i-1]
                delta_u = u[j, i-1] - u[j, i-2]
                r_w = delta_u / (delta_d + np.sign(delta_d + eps) * eps)
            else:
                r_w = 0
        else:
            # East face (i+1/2): flow from right
            # r_e⁻ = (u_E - u_EE)/(u_C - u_E)
            if i < N - 2:
                delta_d = u[j, i] - u[j, i+1]
                delta_u = u[j, i+1] - u[j, i+2]
                r_e = delta_u / (delta_d + np.sign(delta_d + eps) * eps)
            else:
                r_e = 0

            # West face (i-1/2): flow from right
            # r_w⁻ = (u_C - u_E)/(u_W - u_C)
            if i >= 1 and i < N - 1:
                delta_d = u[j, i-1] - u[j, i]
                delta_u = u[j, i] - u[j, i+1]
                r_w = delta_u / (delta_d + np.sign(delta_d + eps) * eps)
            else:
                r_w = 0

        return r_e, r_w

    def _compute_gradient_ratio_y(self, u: np.ndarray, i: int, j: int,
                                   positive_flow: bool) -> Tuple[float, float]:
        """
        Compute gradient ratios for TVD scheme in y-direction.

        Similar to x-direction but for north/south faces.

        For positive flow (b > 0) at north face:
            r_n⁺ = (u_C - u_S)/(u_N - u_C)

        For negative flow (b < 0) at north face:
            r_n⁻ = (u_N - u_NN)/(u_C - u_N)
        """
        N = self.N
        eps = 1e-10

        if positive_flow:
            # North face (j+1/2): flow from below
            if j >= 1 and j < N - 1:
                delta_d = u[j+1, i] - u[j, i]
                delta_u = u[j, i] - u[j-1, i]
                r_n = delta_u / (delta_d + np.sign(delta_d + eps) * eps)
            else:
                r_n = 0

            # South face (j-1/2): flow from below
            if j >= 2:
                delta_d = u[j, i] - u[j-1, i]
                delta_u = u[j-1, i] - u[j-2, i]
                r_s = delta_u / (delta_d + np.sign(delta_d + eps) * eps)
            else:
                r_s = 0
        else:
            # North face (j+1/2): flow from above
            if j < N - 2:
                delta_d = u[j, i] - u[j+1, i]
                delta_u = u[j+1, i] - u[j+2, i]
                r_n = delta_u / (delta_d + np.sign(delta_d + eps) * eps)
            else:
                r_n = 0

            # South face (j-1/2): flow from above
            if j >= 1 and j < N - 1:
                delta_d = u[j-1, i] - u[j, i]
                delta_u = u[j, i] - u[j+1, i]
                r_s = delta_u / (delta_d + np.sign(delta_d + eps) * eps)
            else:
                r_s = 0

        return r_n, r_s

    def _compute_face_value_tvd(self, u_C: float, u_D: float, r: float) -> float:
        """
        Compute face value using TVD interpolation.

        TVD FACE VALUE FORMULA:
        -----------------------
        u_f = u_C + (1/2)*ψ(r)*(u_D - u_C)

        Where:
            u_C = upwind cell value
            u_D = downwind cell value
            r = gradient ratio
            ψ = flux limiter function

        For upwind scheme: ψ = 0, giving u_f = u_C
        For central difference: ψ = 1, giving u_f = (u_C + u_D)/2
        """
        psi_r = self.psi(np.array([r]))[0]
        return u_C + 0.5 * psi_r * (u_D - u_C)

    def _compute_fvm_flux_divergence_upwind(self, u: np.ndarray) -> np.ndarray:
        """Compute flux divergence using EXACT method from solver.py (reference)."""
        N = self.N
        flux_div = np.zeros_like(u)
        cell_area = self.dx * self.dy

        k_e, k_w = self.dy, -self.dy
        l_n, l_s = self.dx, -self.dx

        for j in range(1, N - 1):
            for i in range(1, N - 1):
                u_C = u[j, i]
                u_E = u[j, i + 1]
                u_W = u[j, i - 1]
                u_N = u[j + 1, i]
                u_S = u[j - 1, i]

                # East face
                u_e = 0.5 * (u_C + u_E)
                phi_e = u_C if self.a * u_e >= 0 else u_E
                flux_e = k_e * (self.a / 2.0) * phi_e**2

                # West face
                u_w = 0.5 * (u_C + u_W)
                phi_w = u_W if self.a * u_w >= 0 else u_C
                flux_w = k_w * (self.a / 2.0) * phi_w**2

                # North face
                phi_n = u_C if self.b_coeff >= 0 else u_N
                flux_n = l_n * self.b_coeff * phi_n

                # South face
                phi_s = u_S if self.b_coeff >= 0 else u_C
                flux_s = l_s * self.b_coeff * phi_s

                flux_div[j, i] = (flux_e + flux_w + flux_n + flux_s) / cell_area

        return flux_div

    def _compute_flux_x_tvd(self, u: np.ndarray) -> np.ndarray:
        """Compute x-direction flux using TVD scheme."""
        N = self.N
        flux = np.zeros_like(u)

        for j in range(1, N - 1):
            for i in range(1, N - 1):
                # East face (i+1/2)
                u_face_e = 0.5 * (u[j, i] + u[j, i+1])
                wave_speed_e = self.a * u_face_e

                if wave_speed_e >= 0:
                    r_e, _ = self._compute_gradient_ratio_x(u, i, j, True)
                    u_e = self._compute_face_value_tvd(u[j, i], u[j, i+1], r_e)
                else:
                    r_e, _ = self._compute_gradient_ratio_x(u, i, j, False)
                    u_e = self._compute_face_value_tvd(u[j, i+1], u[j, i], r_e)

                # West face (i-1/2)
                u_face_w = 0.5 * (u[j, i-1] + u[j, i])
                wave_speed_w = self.a * u_face_w

                if wave_speed_w >= 0:
                    _, r_w = self._compute_gradient_ratio_x(u, i, j, True)
                    u_w = self._compute_face_value_tvd(u[j, i-1], u[j, i], r_w)
                else:
                    _, r_w = self._compute_gradient_ratio_x(u, i, j, False)
                    u_w = self._compute_face_value_tvd(u[j, i], u[j, i-1], r_w)

                # Net flux in x-direction
                E_e = 0.5 * self.a * u_e**2
                E_w = 0.5 * self.a * u_w**2
                flux[j, i] = (E_e - E_w) / self.dx

        return flux

    def _compute_flux_y_tvd(self, u: np.ndarray) -> np.ndarray:
        """Compute y-direction flux using TVD scheme."""
        N = self.N
        flux = np.zeros_like(u)

        for j in range(1, N - 1):
            for i in range(1, N - 1):
                # Wave speed in y-direction (constant)
                wave_speed = self.b_coeff

                # North face (j+1/2)
                if wave_speed >= 0:
                    r_n, _ = self._compute_gradient_ratio_y(u, i, j, True)
                    u_n = self._compute_face_value_tvd(u[j, i], u[j+1, i], r_n)
                else:
                    r_n, _ = self._compute_gradient_ratio_y(u, i, j, False)
                    u_n = self._compute_face_value_tvd(u[j+1, i], u[j, i], r_n)

                # South face (j-1/2)
                if wave_speed >= 0:
                    _, r_s = self._compute_gradient_ratio_y(u, i, j, True)
                    u_s = self._compute_face_value_tvd(u[j-1, i], u[j, i], r_s)
                else:
                    _, r_s = self._compute_gradient_ratio_y(u, i, j, False)
                    u_s = self._compute_face_value_tvd(u[j, i], u[j-1, i], r_s)

                # Net flux in y-direction
                F_n = self.b_coeff * u_n
                F_s = self.b_coeff * u_s
                flux[j, i] = (F_n - F_s) / self.dy

        return flux

    def _compute_flux_x_upwind(self, u: np.ndarray) -> np.ndarray:
        """Compute x-direction flux using First-Order Upwind FVM scheme."""
        N = self.N
        flux = np.zeros_like(u)

        # Surface areas for FVM
        k_e = self.dy
        k_w = -self.dy
        cell_area = self.dx * self.dy

        # Loop over interior cells
        for j in range(1, N - 1):
            for i in range(1, N - 1):
                u_C = u[j, i]
                u_E = u[j, i + 1]
                u_W = u[j, i - 1]

                # East face (i+1/2)
                u_e = 0.5 * (u_C + u_E)
                phi_e = u_C if self.a * u_e >= 0 else u_E
                flux_e = k_e * (self.a / 2.0) * phi_e**2

                # West face (i-1/2)
                u_w = 0.5 * (u_C + u_W)
                phi_w = u_W if self.a * u_w >= 0 else u_C
                flux_w = k_w * (self.a / 2.0) * phi_w**2

                # Net flux divergence in x-direction
                flux[j, i] = (flux_e + flux_w) / cell_area

        return flux

    def _compute_flux_y_upwind(self, u: np.ndarray) -> np.ndarray:
        """Compute y-direction flux using First-Order Upwind FVM scheme."""
        N = self.N
        flux = np.zeros_like(u)

        # Surface areas for FVM
        l_n = self.dx
        l_s = -self.dx
        cell_area = self.dx * self.dy

        for j in range(1, N - 1):
            for i in range(1, N - 1):
                u_C = u[j, i]
                u_N = u[j + 1, i]
                u_S = u[j - 1, i]

                # North face (j+1/2)
                phi_n = u_C if self.b_coeff >= 0 else u_N
                flux_n = l_n * self.b_coeff * phi_n

                # South face (j-1/2)
                phi_s = u_S if self.b_coeff >= 0 else u_C
                flux_s = l_s * self.b_coeff * phi_s

                # Net flux divergence in y-direction
                flux[j, i] = (flux_n + flux_s) / cell_area

        return flux

    def _compute_flux_x_upwind_vectorized(self, u: np.ndarray) -> np.ndarray:
        """Vectorized x-direction flux computation using upwind scheme."""
        flux = np.zeros_like(u)

        # Interior slice indices
        j_int = slice(1, self.N - 1)
        i_int = slice(1, self.N - 1)

        # Cell center values for wave speed determination
        u_center = u[j_int, i_int]
        wave_speed = self.a * u_center

        # Backward difference: (u_C - u_W)/dx
        du_backward = (u[j_int, i_int] - u[j_int, :-2]) / self.dx

        # Forward difference: (u_E - u_C)/dx
        du_forward = (u[j_int, 2:] - u[j_int, i_int]) / self.dx

        # Upwind selection based on wave speed
        flux[j_int, i_int] = np.where(
            wave_speed >= 0,
            self.a * u_center * du_backward,
            self.a * u_center * du_forward
        )

        return flux

    def _compute_flux_y_upwind_vectorized(self, u: np.ndarray) -> np.ndarray:
        """Vectorized y-direction flux computation using upwind scheme."""
        flux = np.zeros_like(u)

        j_int = slice(1, self.N - 1)
        i_int = slice(1, self.N - 1)

        if self.b_coeff >= 0:
            flux[j_int, i_int] = self.b_coeff * (
                u[j_int, i_int] - u[:-2, i_int]
            ) / self.dy
        else:
            flux[j_int, i_int] = self.b_coeff * (
                u[2:, i_int] - u[j_int, i_int]
            ) / self.dy

        return flux

    def time_step_explicit(self, use_vectorized: bool = True) -> float:
        """
        Perform one explicit Euler time step.

        Parameters
        ----------
        use_vectorized : bool
            Use vectorized (faster) flux computation when available

        Returns
        -------
        residual : float
            L2 norm of the change in solution
        """
        u_old = self.u.copy()
        self._compute_time_step()

        # Compute fluxes based on selected scheme
        if self.limiter == FluxLimiter.UPWIND:
            # Use EXACT FVM flux divergence from solver.py (PROVEN CORRECT)
            flux_div = self._compute_fvm_flux_divergence_upwind(self.u)
        else:
            # Use TVD scheme with selected flux limiter
            flux_x = self._compute_flux_x_tvd(self.u)
            flux_y = self._compute_flux_y_tvd(self.u)
            flux_div = flux_x + flux_y

        # Explicit Euler update
        self.u = u_old - self.dt * flux_div

        # Apply boundary conditions
        self._apply_boundary_conditions()

        # Apply Neumann BC at top explicitly (zero gradient)
        self.u[-1, 1:-1] = self.u[-2, 1:-1]

        # Compute residual (L2 norm of change)
        diff = self.u - u_old
        residual = np.sqrt(np.sum(diff**2)) / self.N

        return residual

    def _build_tdma_coefficients(self, j: int, u: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Build TDMA coefficient arrays for row j.
        Returns: (a_lower, a_diag, a_upper, b_rhs)
        """
        N = self.N

        # Initialize coefficient arrays for interior points only
        n_interior = N - 2
        a_lower = np.zeros(n_interior)  # a_W (lower diagonal)
        a_diag = np.zeros(n_interior)   # a_C (main diagonal)
        a_upper = np.zeros(n_interior)  # a_E (upper diagonal)
        b_rhs = np.zeros(n_interior)    # Right-hand side

        # S/Δt coefficient
        S_dt = self.dS / self.dt

        for k, i in enumerate(range(1, N - 1)):  # Interior columns
            u_C = u[j, i]

            # Wave speed at cell center for linearization
            wave_speed = self.a * u_C

            # Face coefficients
            k_e = self.dy  # East face area
            k_w = -self.dy  # West face area (negative for outward normal)

            # Upwind direction
            if wave_speed >= 0:
                # Flow from left to right
                # u_e = u_C, u_w = u_W
                coeff_C_from_e = k_e * self.a * u_C
                coeff_E_from_e = 0
                coeff_C_from_w = 0
                coeff_W_from_w = k_w * self.a * u[j, i-1]
            else:
                # Flow from right to left
                # u_e = u_E, u_w = u_C
                coeff_C_from_e = 0
                coeff_E_from_e = k_e * self.a * u[j, i+1]
                coeff_C_from_w = k_w * self.a * u_C
                coeff_W_from_w = 0

            # Build coefficients
            a_diag[k] = S_dt + coeff_C_from_e - coeff_C_from_w

            # East neighbor (upper diagonal)
            if i < N - 2:  # Not at right interior boundary
                a_upper[k] = coeff_E_from_e

            # West neighbor (lower diagonal)
            if i > 1:  # Not at left interior boundary
                a_lower[k] = -coeff_W_from_w

            # RHS: time term minus y-flux (treated explicitly)
            # Y-flux
            if self.b_coeff >= 0:
                flux_y = self.b_coeff * (u[j, i] - u[j-1, i]) / self.dy
            else:
                flux_y = self.b_coeff * (u[j+1, i] - u[j, i]) / self.dy

            b_rhs[k] = S_dt * u_C - self.dS * flux_y

            # Boundary contributions
            if i == 1:
                # Left boundary: add flux from west BC
                u_bc_w = self.c
                flux_bc_w = k_w * self.a * 0.5 * u_bc_w**2
                b_rhs[k] -= flux_bc_w

            if i == N - 2:
                # Right boundary: add flux from east BC
                u_bc_e = self.d
                flux_bc_e = k_e * self.a * 0.5 * u_bc_e**2
                b_rhs[k] -= flux_bc_e

        return a_lower, a_diag, a_upper, b_rhs

    def time_step_implicit_tdma(self) -> float:
        """Perform one implicit time step using TDMA (line-by-line)."""
        u_old = self.u.copy()
        self._compute_time_step()

        # Solve for each row (line-by-line TDMA)
        for j in range(1, self.N - 1):
            a_lower, a_diag, a_upper, b_rhs = self._build_tdma_coefficients(j, u_old)

            # Solve tri-diagonal system
            u_interior = solve_tdma(a_lower, a_diag, a_upper, b_rhs)

            # Update solution for this row
            self.u[j, 1:-1] = u_interior

        # Apply boundary conditions
        self._apply_boundary_conditions()
        self.u[-1, 1:-1] = self.u[-2, 1:-1]  # Neumann at top

        # Compute residual
        diff = self.u - u_old
        residual = np.sqrt(np.sum(diff**2)) / self.N

        return residual

    def time_step(self, method: SolverMethod = SolverMethod.EXPLICIT_EULER) -> float:
        """
        Perform one time step using the specified method.

        Parameters
        ----------
        method : SolverMethod
            The solver method to use

        Returns
        -------
        residual : float
            L2 norm of the change in solution
        """
        if method == SolverMethod.EXPLICIT_EULER:
            return self.time_step_explicit()
        elif method == SolverMethod.IMPLICIT_TDMA:
            return self.time_step_implicit_tdma()
        else:
            raise ValueError(f"Unknown solver method: {method}")

    def solve_steady_state(self, max_iterations: int = 10000,
                           tolerance: float = 1e-6,
                           print_interval: int = 500,
                           method: SolverMethod = SolverMethod.EXPLICIT_EULER) -> Dict:
        """
        Solve until steady state is reached.
        """
        print(f"Solving steady-state (Limiter: {self.limiter.value}, Method: {method.value})...")
        
        start_time = time.time()
        self.residual_history = []
        converged = False

        for iteration in range(1, max_iterations + 1):
            residual = self.time_step(method)
            self.residual_history.append(residual)
            self.iteration_count = iteration

            if iteration % print_interval == 0:
                print(f"  Iter {iteration}: residual = {residual:.6e}")

            if residual < tolerance:
                converged = True
                print(f"  Converged at iter {iteration}, residual: {residual:.6e}")
                break

        elapsed_time = time.time() - start_time
        if not converged:
             print(f"  Warning: Max iterations ({max_iterations}) reached. Final residual: {residual:.6e}")

        return {
            'converged': converged,
            'iterations': self.iteration_count,
            'final_residual': self.residual_history[-1] if self.residual_history else 0,
            'residual_history': self.residual_history.copy(),
            'time_elapsed': elapsed_time,
            'limiter': self.limiter.value,
            'method': method.value
        }

    def get_solution(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Get the current solution."""
        return self.X, self.Y, self.u.copy()

    def compute_error(self, analytical_solution: callable) -> Dict:
        """Compute error metrics against an analytical solution."""
        u_analytical = analytical_solution(self.X, self.Y)
        diff = self.u - u_analytical

        L1 = np.mean(np.abs(diff))
        L2 = np.sqrt(np.mean(diff**2))
        Linf = np.max(np.abs(diff))
        total_error = np.sum(np.abs(diff))

        u_analytical_norm = np.sqrt(np.mean(u_analytical**2))
        relative_L2 = L2 / u_analytical_norm if u_analytical_norm > 1e-10 else L2

        return {
            'L1': L1,
            'L2': L2,
            'Linf': Linf,
            'relative_L2': relative_L2,
            'total_error': total_error,
            'averaged_error': total_error / (self.N * self.N)
        }

    def plot_solution(self, save_path: str = None, show_plot: bool = True,
                      title_suffix: str = "") -> plt.Figure:
        """Plot the numerical solution."""
        fig = plt.figure(figsize=(16, 5))

        # Subplot 1: Filled contour
        ax1 = fig.add_subplot(131)
        contourf = ax1.contourf(self.X, self.Y, self.u, levels=20, cmap='RdYlBu_r')
        ax1.set_xlabel('x', fontsize=12)
        ax1.set_ylabel('y', fontsize=12)
        ax1.set_title(f'Numerical Solution\n({self.limiter.value}, N={self.N})'
                      f'\n{title_suffix}', fontsize=11)
        ax1.set_aspect('equal')
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

    def plot_comparison(self, analytical_solution: callable,
                        save_path: str = None, show_plot: bool = True) -> plt.Figure:
        """Plot comparison between numerical and analytical solutions."""
        u_analytical = analytical_solution(self.X, self.Y)
        error = self.u - u_analytical

        fig = plt.figure(figsize=(18, 5))

        # Analytical
        ax1 = fig.add_subplot(131)
        cf1 = ax1.contourf(self.X, self.Y, u_analytical, levels=20, cmap='RdYlBu_r')
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_title('Analytical Solution')
        ax1.set_aspect('equal')
        plt.colorbar(cf1, ax=ax1)

        # Numerical
        ax2 = fig.add_subplot(132)
        cf2 = ax2.contourf(self.X, self.Y, self.u, levels=20, cmap='RdYlBu_r')
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        ax2.set_title(f'Numerical ({self.limiter.value})')
        ax2.set_aspect('equal')
        plt.colorbar(cf2, ax=ax2)

        # Error
        ax3 = fig.add_subplot(133)
        cf3 = ax3.contourf(self.X, self.Y, error, levels=20, cmap='coolwarm')
        ax3.set_xlabel('x')
        ax3.set_ylabel('y')
        ax3.set_title(f'Error (N={self.N})')
        ax3.set_aspect('equal')
        plt.colorbar(cf3, ax=ax3)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')

        if show_plot:
            plt.show()

        return fig

    def print_summary(self):
        """Print solver summary."""
        print("\n" + "="*70)
        print("Burgers Solver Summary (FVM-TVD)")
        print("="*70)
        print(f"Grid: {self.N} x {self.N}")
        print(f"Parameters: a={self.a}, b={self.b_coeff}, c={self.c}, d={self.d}")
        print(f"Grid spacing: dx={self.dx:.6f}, dy={self.dy:.6f}")
        print(f"Flux limiter: {self.limiter.value}")
        print(f"Current dt: {self.dt:.6e}")
        print("-"*70)
        print("Solution statistics:")
        print(f"  Min: {self.u.min():.6f}")
        print(f"  Max: {self.u.max():.6f}")
        print(f"  Mean: {self.u.mean():.6f}")
        if self.residual_history:
            print(f"  Final residual: {self.residual_history[-1]:.6e}")
            print(f"  Iterations: {self.iteration_count}")
        print("="*70 + "\n")


# =================================================================================
# DEMONSTRATION AND TESTING
# =================================================================================

def run_limiter_comparison(N: int = 21, max_iterations: int = 10000,
                           tolerance: float = 1e-6) -> Dict:
    """
    Run comparison of different flux limiters.

    Parameters
    ----------
    N : int
        Grid size
    max_iterations : int
        Maximum iterations
    tolerance : float
        Convergence tolerance

    Returns
    -------
    results : dict
        Dictionary of results for each limiter
    """
    from src.analytical.case1_solution import analytical_solution_case1

    limiters = [
        FluxLimiter.UPWIND,
        FluxLimiter.MINMOD,
        FluxLimiter.VAN_ALBADA,
        FluxLimiter.UMIST,
        FluxLimiter.SUPERBEE,
        FluxLimiter.VAN_LEER,
        FluxLimiter.KOREN,
    ]

    results = {}

    print("\n" + "#"*70)
    print(f" Flux Limiter Comparison (N={N})")
    print("#"*70)

    for limiter in limiters:
        print(f"\n--- Testing {limiter.value} ---")
        solver = BurgersSolverFVM_TVD(
            N=N, a=1.0, b=1.0, c=1.5, d=-0.5,
            cfl=0.5, limiter=limiter
        )

        result = solver.solve_steady_state(
            max_iterations=max_iterations,
            tolerance=tolerance,
            print_interval=2000
        )

        errors = solver.compute_error(analytical_solution_case1)
        result['errors'] = errors
        result['solver'] = solver

        results[limiter.value] = result

        print(f"  Total error: {errors['total_error']:.4f}")
        print(f"  Averaged error: {errors['averaged_error']:.6f}")

    return results


def print_comparison_table(results: Dict):
    """Print a comparison table of limiter results."""
    print("\n" + "="*80)
    print("FLUX LIMITER COMPARISON TABLE")
    print("="*80)
    print(f"{'Limiter':<15} {'Iterations':>12} {'Total Error':>15} {'Avg Error':>12} {'L2 Error':>12}")
    print("-"*80)

    for name, result in results.items():
        print(f"{name:<15} {result['iterations']:>12d} "
              f"{result['errors']['total_error']:>15.4f} "
              f"{result['errors']['averaged_error']:>12.6f} "
              f"{result['errors']['L2']:>12.6f}")

    print("="*80)


if __name__ == "__main__":
    print("\n" + "="*70)
    print(" FVM-TVD Solver for 2D Inviscid Burgers Equation")
    print("="*70)

    # Test with N=21 using different limiters
    print("\n--- Testing with N=21 ---\n")

    # UMIST limiter test
    solver = BurgersSolverFVM_TVD(
        N=21, a=1.0, b=1.0, c=1.5, d=-0.5,
        cfl=0.5, limiter=FluxLimiter.UMIST
    )

    result = solver.solve_steady_state(
        max_iterations=10000,
        tolerance=1e-6,
        print_interval=1000
    )

    solver.print_summary()

    print("\nTest completed successfully!")
