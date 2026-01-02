# 2D Inviscid Burgers Equation CFD Project

A collaborative computational fluid dynamics (CFD) project implementing numerical and analytical solutions for the 2D inviscid Burgers equation.

## Project Overview

This project solves the 2D inviscid Burgers equation using:
- **Analytical solutions**: Direct mathematical solutions
- **Numerical methods**: Finite difference, finite volume, or finite element approaches

The equation being solved is:
```
∂u/∂t + u·∇u = 0  (in 2D domain)
```

## Project Structure

```
├── src/
│   ├── analytical/          # Analytical solution implementations
│   └── numerical/           # Numerical method implementations
├── plots/
│   ├── analytical/          # Plots from analytical solutions
│   └── numerical/           # Plots from numerical solutions
├── data/                    # Simulation data and results
├── requirements.txt         # Python dependencies
├── README.md               # This file
└── .gitignore             # Git ignore rules
```

## Installation

### Prerequisites
- Python 3.8 or higher
- Git (command-line) or [GitHub Desktop](https://desktop.github.com) (recommended for beginners)

### Setup Steps

1. **Clone the repository** (after it's created on GitHub):
```bash
git clone https://github.com/your-username/burgers-equation-cfd.git
cd burgers-equation-cfd
```

2. **Initialize the project structure**:
```bash
python init_project.py
```

3. **Create a virtual environment** (recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

4. **Install dependencies**:
```bash
pip install -r requirements.txt
```

## Collaborative Workflow

**For detailed GitHub Desktop and Git instructions, see [GIT_WORKFLOW.md](GIT_WORKFLOW.md)**

### Branch Strategy

Two main development branches:
- **`analytical-branch`**: For analytical solution development
- **`numerical-branch`**: For numerical method development
- **`main`**: Stable, integrated version

### Getting Started with Git (Command-Line)

1. **Initial GitHub Setup** (first time only):
```bash
# Initialize git (if not already initialized)
git init

# Add all files
git add .

# Create initial commit
git commit -m "Initial project setup"

# Add remote repository (replace with your GitHub URL)
git remote add origin https://github.com/your-username/burgers-equation-cfd.git

# Push to GitHub
git push -u origin main
```

2. **Creating branches** (each developer):
```bash
# Create and switch to your working branch
git checkout -b analytical-branch  # or numerical-branch

# Push the branch to remote
git push -u origin analytical-branch
```

3. **Daily workflow**:
```bash
# Check current branch
git branch

# Make changes to your files
# ... edit code ...

# Stage and commit your changes
git add .
git commit -m "Describe your changes here"

# Push to your branch
git push origin analytical-branch  # or numerical-branch

# Pull latest changes from your branch
git pull origin analytical-branch
```

4. **Merging back to main** (when ready):
```bash
# Switch to main branch
git checkout main

# Pull latest main
git pull origin main

# Merge your branch
git merge analytical-branch  # or numerical-branch

# Push to main
git push origin main
```

## Dependencies

- **NumPy**: Numerical computing
- **Matplotlib**: Visualization and plotting
- **SciPy**: Scientific computing utilities
- **Jupyter** (optional): Interactive notebooks for analysis

See `requirements.txt` for versions.

## Usage Examples

### Running Analytical Solution
```python
from src.analytical.solver import solve_burgers_analytical
import matplotlib.pyplot as plt

# Define domain
x = np.linspace(0, 1, 100)
y = np.linspace(0, 1, 100)
t = 0.5

# Solve
u, v = solve_burgers_analytical(x, y, t)

# Plot
plt.contourf(x, y, u)
plt.colorbar()
plt.show()
```

### Running Numerical Solution
```python
from src.numerical.solver import solve_burgers_numerical
import matplotlib.pyplot as plt

# Parameters
domain = {'x': (0, 1, 100), 'y': (0, 1, 100)}
t_final = 0.5
dt = 0.001

# Solve
u, v, t_array = solve_burgers_numerical(domain, t_final, dt)

# Plot
plt.contourf(u[:,:,-1])
plt.colorbar()
plt.show()
```

## File Descriptions

### Source Files
- `src/analytical/solver.py`: Analytical solution implementation
- `src/numerical/solver.py`: Numerical method implementation
- `src/utils.py`: Shared utility functions

### Output Files
- `plots/analytical/`: Generated analytical solution plots
- `plots/numerical/`: Generated numerical solution plots
- `data/`: Stored simulation results (CSV, HDF5, etc.)

## Git Workflow Best Practices

**See [GIT_WORKFLOW.md](GIT_WORKFLOW.md) for complete GitHub Desktop and command-line workflows.**

1. **Before starting work**: Always pull latest changes
   - GitHub Desktop: Click "Fetch origin" → "Pull origin"
   - Command-line: `git pull origin analytical-branch`

2. **Commit frequently** with meaningful messages
   - GitHub Desktop: Use the commit form at bottom-left
   - Command-line: `git commit -m "Implement finite difference scheme"`

3. **Avoid conflicts**: Keep changes localized to your files
   - Analytical work → `src/analytical/`
   - Numerical work → `src/numerical/`

4. **Communicate**: Let your partner know before major merges

5. **Review code**: Before merging to main, review each other's work

## Troubleshooting

### Merge Conflicts
If conflicts occur:
```bash
# View conflicting files
git status

# Open conflicted files and resolve manually
# After resolving:
git add .
git commit -m "Resolve merge conflicts"
git push origin main
```

### Accidentally Committed to Wrong Branch
```bash
# Undo last commit (keep changes)
git reset --soft HEAD~1

# Switch to correct branch
git checkout correct-branch

# Commit again
git commit -m "Your message"
```

## Contributing

1. Create a feature branch from your main branch
2. Make changes and commit regularly
3. Push to remote and create a pull request
4. After review, merge to main
5. Delete the feature branch

## License

[Specify your license here]

## Contact & Contributors

- Partner 1: [Name/Email]
- Partner 2: [Name/Email]

## References

[Add relevant papers or resources about Burgers equation]
