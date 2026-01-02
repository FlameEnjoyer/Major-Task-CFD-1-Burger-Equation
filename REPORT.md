# Answers to Project Tasks (CFD 1 2025)

This document provides the specific answers and explanations for the 8 tasks outlined in the "Statement of Final task of CFD1 2025.pdf".

---

## 1. Analytical Verification
**Objective:** Implement and plot the piecewise analytical solution for Case 1 ($a=1, b=1, c=1.5, d=-0.5$).

**Solution:**
The analytical solution is derived using the method of characteristics. For the given boundary conditions, it contains shock waves (discontinuities) where characteristics intersect.

*   **Implementation**: See `src/analytical/case1_solution.py`.
*   **Results**: The plot `plots/analytical/case1_contour.png` visualizes this solution. You can clearly see the piecewise constant regions ($u=1.5$ and $u=-0.5$) separated by expansion fans or shocks.

## 2. Numerical Derivation (FVM)
**Objective:** Derive the discrete equations using the Finite Volume Method (FVM).

**Derivation:**
Starting with the conservative form of the 2D Burgers equation:
$$ \frac{\partial u}{\partial t} + \frac{\partial F(u)}{\partial x} + \frac{\partial G(u)}{\partial y} = 0 $$
Where flux vectors are $F(u) = \frac{a u^2}{2}$ and $G(u) = b u$.

Integrating over a control volume $V_{i,j}$ and applying the Divergence Theorem:
$$ \int_{V} \frac{\partial u}{\partial t} dV + \oint_{S} \mathbf{H} \cdot \mathbf{n} dS = 0 $$
Where $\mathbf{H} = [F, G]^T$.

Discretizing on a uniform grid ($\Delta x, \Delta y$):
$$ \frac{u_{i,j}^{n+1} - u_{i,j}^n}{\Delta t} \Delta x \Delta y + \left[ (F_{i+1/2,j} - F_{i-1/2,j}) \Delta y + (G_{i,j+1/2} - G_{i,j-1/2}) \Delta x \right] = 0 $$

Rearranging for the explicit update rule:
$$ u_{i,j}^{n+1} = u_{i,j}^n - \frac{\Delta t}{\Delta x} (F_{i+1/2} - F_{i-1/2}) - \frac{\Delta t}{\Delta y} (G_{j+1/2} - G_{j-1/2}) $$

## 3. Algorithm Design (Upwind Scheme)
**Objective:** Develop an algorithm using a first-order Upwind Scheme.

**Explanation:**
To calculate the interface fluxes ($F_{i+1/2}$), we cannot simply average the values because the equation is hyperbolic (information travels in specific directions). We use the **Upwind Scheme** based on the local wave speed ($a \cdot u$).

**Algorithm Logic for $F_{i+1/2}$:**
*   Check wave speed $c = a \cdot u_{i+1/2}$.
*   **Case 1 ($c > 0$):** Flow moves Left $\to$ Right.
    $$ F_{i+1/2} \approx a \cdot u_i \cdot u_i = a (u_i)^2 $$
    *(Implementation uses linearized form $a u \frac{\partial u}{\partial x} \approx a u_i \frac{u_i - u_{i-1}}{\Delta x}$)*
*   **Case 2 ($c < 0$):** Flow moves Right $\to$ Left.
    $$ F_{i+1/2} \approx a \cdot u_{i+1} \cdot u_{i+1} $$

**Stability (CFL Condition):**
The time step $\Delta t$ is controlled dynamically to ensure stability:
$$ \Delta t \le \text{CFL} \cdot \min \left( \frac{\Delta x}{|a u|_{max}}, \frac{\Delta y}{|b|} \right) $$
We used $\text{CFL} = 0.5$.

## 4. Programming (Grid Grids)
**Objective:** Solve on structured grids of $21 \times 21$ and $41 \times 41$.

**Implementation:**
The code supports arbitrary $N$.
*   **Result**: The outputs `solution_N21.png` and `solution_N41.png` show the steady-state results.
*   **Observation**: The grid for $N=41$ (finer) produces sharper shock resolution compared to $N=21$ (coarser), which is more diffusive.

## 5. Convergence
**Objective:** Track and plot convergence of the solution over time.

**Methodology:**
Convergence is defined as the system reaching a Steady State ($\partial u / \partial t \approx 0$).
*   **Metric:** L2 Norm of the Residual ($R = \sqrt{\sum (u^{n+1} - u^n)^2}$).
*   **Criterion:** Stop when $R < 10^{-6}$.
*   **Result**: See `convergence_#.png`. The plots show the residual decreasing exponentially (linear on log-scale) until hitting the tolerance floor.

## 6. Validation
**Objective:** Compare numerical results against the analytical solution.

**Comparison:**
We computed the error $E = |u_{num} - u_{exact}|$.
*   **L2 Error (N=21):** ~3.5e-2
*   **L2 Error (N=41):** ~2.5e-2 (Error decreases as grid is refined).
*   **Visual Check:** Comparisons are plotted in `comparison_N41.png`. The numerical solution matches the analytical "plateaus" perfectly ($u=1.5, -0.5$). The main error is located exactly at the shock lines, where the numerical scheme "smears" the discontinuity over 2-3 cells.

## 7. Parametric Study
**Objective:** Run the solver for three additional cases with varying coefficients $a, b, c, d$.

**Cases Run:**
1.  **Case 1 (Reference):** $a=1, b=1 \to$ Diagonal shock.
2.  **Case 2 ($a=1.5, b=0.5$):** Stronger horizontal advection. Discontinuities are steeper/more vertical.
3.  **Case 3 ($a=1.5, b=2.0$):** Strong vertical advection. Discontinuities are more horizontal.
4.  **Case 4 ($d=-0.6$):** Stronger jump condition at the boundary.

**Results:**
All comparative plots are saved in `plots/multicase/`. The solver successfully converged for all parameter variations, demonstrating the robustness of the implementation.
