#!/usr/bin/env python3
"""
Project initialization script for 2D Inviscid Burgers Equation CFD
Creates the necessary directory structure and files
"""

import os
from pathlib import Path

def create_project_structure():
    """Create the project directory structure"""

    # Define base directory paths
    base_dir = Path.cwd()
    directories = [
        'src',
        'plots',
        'data',
        'src/analytical',
        'src/numerical',
        'plots/analytical',
        'plots/numerical',
    ]

    # Create directories
    print("Creating project directories...")
    for directory in directories:
        dir_path = base_dir / directory
        dir_path.mkdir(parents=True, exist_ok=True)
        print(f"✓ Created: {directory}/")

    # Create __init__.py files for Python packages
    print("\nCreating Python package files...")
    init_files = [
        'src/__init__.py',
        'src/analytical/__init__.py',
        'src/numerical/__init__.py',
    ]

    for init_file in init_files:
        file_path = base_dir / init_file
        if not file_path.exists():
            file_path.touch()
            print(f"✓ Created: {init_file}")

    # Create placeholder files
    print("\nCreating placeholder files...")

    # .gitkeep files to track empty directories
    gitkeep_files = [
        'plots/analytical/.gitkeep',
        'plots/numerical/.gitkeep',
        'data/.gitkeep',
    ]

    for gitkeep_file in gitkeep_files:
        file_path = base_dir / gitkeep_file
        if not file_path.exists():
            file_path.touch()
            print(f"✓ Created: {gitkeep_file}")

    print("\n✓ Project structure initialized successfully!")
    print("\nDirectory structure:")
    print("""
    project-root/
    ├── src/
    │   ├── __init__.py
    │   ├── analytical/
    │   │   └── __init__.py
    │   └── numerical/
    │       └── __init__.py
    ├── plots/
    │   ├── analytical/
    │   ├── numerical/
    │   └── (generated plots saved here)
    ├── data/
    │   └── (simulation data and results)
    ├── requirements.txt
    ├── README.md
    └── .gitignore
    """)

if __name__ == '__main__':
    create_project_structure()
