import yaml
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Required for 3D plotting
import numpy as np
import os
import sys

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.numerical.solver import BurgersSolver2D

def load_cases(filename='cases.yaml'):
    with open(filename, 'r') as f:
        return yaml.safe_load(f)

def run_simulation(case_name, params, grid_sizes=[21, 41, 201], output_dir='results'):
    os.makedirs(output_dir, exist_ok=True)
    
    case_history = {}
    solvers = {}
    
    print(f"RUNNING {case_name}")
    print(f"Parameters: {params}")

    for N in grid_sizes:
        print(f"\n  Running Grid Size N={N}...")
        

        solver = BurgersSolver2D(
            N=N, 
            a=params['a'], 
            b=params['b'], 
            c=params['c'], 
            d=params['d'],
            cfl=0.5
        )
        

        result = solver.solve_steady_state(
            max_iterations=30000, 
            tolerance=1e-6, 
            print_interval=5000,
            use_implicit=False,
            use_fvm=True,
            use_vectorized=True
        )
        
        solvers[N] = solver
        
        # Store essential history
        case_history[N] = {
            'converged': bool(result['converged']),
            'iterations': int(result['iterations']),
            'final_residual': float(result['final_residual']),
            'time_elapsed': float(result['time_elapsed']),
            'residual_history': [float(r) for r in result['residual_history']] # Store for yaml if needed, but might be large
        }

    # Generate Plots
    plot_fields(solvers, case_name, params, output_dir)
    plot_convergence(solvers, case_name, params, case_history, output_dir)
    
    yaml_history = {}
    for N in grid_sizes:
        yaml_history[f'N_{N}'] = {
            'converged': case_history[N]['converged'],
            'iterations': case_history[N]['iterations'],
            'final_residual': case_history[N]['final_residual'],
            'time_elapsed': case_history[N]['time_elapsed']
        }
        
    return yaml_history

def plot_fields(solvers, case_name, params, output_dir):
    """
    Figure 1: 6 subplots
    Row 1: 2D field for N=21, 41, 201
    Row 2: 3D field for N=21, 41, 201
    """
    grid_sizes = sorted(solvers.keys())
    
    # Setup Figure: 2 rows (2D, 3D), 3 columns (N sizes)
    fig = plt.figure(figsize=(18, 10))
    
    # Find global min/max for shared colorbar
    vmin = min([s.u.min() for s in solvers.values()])
    vmax = max([s.u.max() for s in solvers.values()])
    
    for i, N in enumerate(grid_sizes):
        solver = solvers[N]
        
        # --- Row 1: 2D Contour ---
        ax_2d = fig.add_subplot(2, 3, i + 1)
        contour = ax_2d.contourf(solver.X, solver.Y, solver.u, levels=25, cmap='RdYlBu_r', vmin=vmin, vmax=vmax)
        ax_2d.set_title(f'2D Field (N={N})', fontsize=12, fontweight='bold')
        ax_2d.set_xlabel('x')
        ax_2d.set_ylabel('y')
        ax_2d.set_aspect('equal')

        # --- Row 2: 3D Surface ---
        ax_3d = fig.add_subplot(2, 3, i + 4, projection='3d')
        surf = ax_3d.plot_surface(solver.X, solver.Y, solver.u, cmap='RdYlBu_r', 
                                  vmin=vmin, vmax=vmax, edgecolor='none', alpha=0.9)
        ax_3d.set_title(f'3D Field (N={N})', fontsize=12, fontweight='bold')
        ax_3d.set_xlabel('x')
        ax_3d.set_ylabel('y')
        ax_3d.set_zlabel('u')
        ax_3d.view_init(elev=30, azim=45)
        
    # Add shared colorbar
    fig.subplots_adjust(right=0.9, top=0.90, hspace=0.3)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    fig.colorbar(contour, cax=cbar_ax, label='u(x,y)')
    
    # Title with params
    param_str = f"a={params['a']}, b={params['b']}, c={params['c']}, d={params['d']}"
    fig.suptitle(f'{case_name}: Solution Fields (Pure Upwind)\n{param_str}', fontsize=16, fontweight='bold')
    
    save_path = os.path.join(output_dir, f'{case_name}_Fields_2D_3D.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved fields plot to {save_path}")

def plot_convergence(solvers, case_name, params, history, output_dir):
    """
    Figure 2: 3 subplots of convergence curve for N=21, 41, 201
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    grid_sizes = sorted(solvers.keys())
    
    for ax, N in zip(axes, grid_sizes):
        res_hist = history[N]['residual_history']
        iters = range(1, len(res_hist) + 1)
        
        ax.semilogy(iters, res_hist, 'b-', linewidth=1.5)
        ax.set_title(f'N = {N}', fontsize=12, fontweight='bold')
        ax.set_xlabel('Iteration')
        ax.set_ylabel('Residual (L2)')
        ax.grid(True, which="both", ls="-", alpha=0.3)
        
        # Annotate final
        final_res = history[N]['final_residual']
        ax.text(0.95, 0.95, f'Final: {final_res:.2e}', 
                transform=ax.transAxes, ha='right', va='top', 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Adjust layout to prevent title overlap
    fig.subplots_adjust(top=0.85)

    param_str = f"a={params['a']}, b={params['b']}, c={params['c']}, d={params['d']}"
    fig.suptitle(f'{case_name}: Residual Convergence (Pure Upwind)\n{param_str}', fontsize=14, fontweight='bold')
    
    save_path = os.path.join(output_dir, f'{case_name}_Convergence.png')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved Convergence plot to {save_path}")

def main():
    cases = load_cases('cases.yaml')
    overall_history = {}
    
    for case_name, params in cases.items():
        case_summary = run_simulation(case_name, params, output_dir='results/cases_output')
        overall_history[case_name] = {
            'parameters': params,
            'results': case_summary
        }
        
    # Save overall history to yaml
    summary_path = 'results/cases_output/simulation_history.yaml'
    with open(summary_path, 'w') as f:
        yaml.dump(overall_history, f, sort_keys=False)
        
    print(f"ALL CASES COMPLETED. Summary saved to {summary_path}")

if __name__ == "__main__":
    main()
