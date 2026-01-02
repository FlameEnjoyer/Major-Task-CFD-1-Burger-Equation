# Git Workflow Guide: GitHub Desktop & Command-Line

**Complete guide for collaborative CFD project using GitHub Desktop (primary) with command-line reference**

---

## Table of Contents

### Part A: GitHub Desktop Workflow (Recommended)
1. [Installation & Setup](#1-installation--setup)
2. [Creating the Repository (Person 1)](#2-creating-the-repository-person-1)
3. [Cloning Repository (Person 2)](#3-cloning-repository-person-2)
4. [Branching Strategy](#4-branching-strategy)
5. [Daily Workflow](#5-daily-workflow)
6. [Merging Changes](#6-merging-changes)
7. [Handling Merge Conflicts](#7-handling-merge-conflicts)
8. [Best Practices](#8-best-practices)

### Part B: Command-Line Reference
9. [Git Commands Quick Reference](#9-git-commands-quick-reference)

### Part C: Troubleshooting
10. [Common Issues & Solutions](#10-common-issues--solutions)

---

## Part A: GitHub Desktop Workflow

# 1. Installation & Setup

## 1.1 Download and Install GitHub Desktop

### For Windows:
1. Visit [desktop.github.com](https://desktop.github.com)
2. Click **"Download for Windows"** (64-bit)
3. Run the downloaded `GitHubDesktopSetup-x64.exe` installer
4. The installer will automatically install and launch GitHub Desktop

### For macOS:
1. Visit [desktop.github.com](https://desktop.github.com)
2. Click **"Download for macOS"**
3. Open the downloaded `.zip` file
4. Drag **GitHub Desktop** to your Applications folder
5. Launch from Applications

## 1.2 Sign In to GitHub

**First-time launch screen:**

1. GitHub Desktop will show a welcome screen
2. Click **"Sign in to GitHub.com"**
3. A browser window opens → Enter your GitHub username and password
4. Click **"Authorize desktop"** to connect GitHub Desktop to your account
5. The browser will show "Success! You may now close this window"
6. Return to GitHub Desktop

**If you don't have a GitHub account:**
1. Click **"Create your free account"** on the welcome screen
2. Fill in username, email, and password on GitHub.com
3. Verify your email address
4. Return to GitHub Desktop and sign in

## 1.3 Configure Git

After signing in, GitHub Desktop will ask you to configure Git:

1. **Name**: Enter your full name (e.g., "John Doe")
   - This appears in commit history
2. **Email**: Use the email associated with your GitHub account
   - Recommendation: Use your GitHub-provided `noreply` email for privacy
   - Find it at: GitHub.com → Settings → Emails → Keep my email private
3. Click **"Finish"** or **"Continue"**

You're now ready to use GitHub Desktop!

---

# 2. Creating the Repository (Person 1)

**This section is for the person who will create and own the repository.**

You have two options:
- **Option A**: Create repository on GitHub.com first, then clone locally *(Recommended)*
- **Option B**: Create repository locally in GitHub Desktop, then publish

## Option A: Create on GitHub.com (Recommended)

### Step 1: Create Repository on GitHub Website

1. Open your browser and go to [GitHub.com](https://github.com)
2. Log in to your account
3. Click the **"+"** icon in the top-right corner
4. Select **"New repository"**

**Repository settings:**
- **Repository name**: `burgers-equation-cfd` (or your preferred name)
- **Description**: "2D Inviscid Burgers Equation CFD Solver"
- **Visibility**:
  - Choose **"Private"** if only you and your partner should access it
  - Choose **"Public"** if you want it publicly visible
- **Initialize repository**:
  - ❌ **Do NOT check** "Add a README file"
  - ❌ **Do NOT check** "Add .gitignore"
  - ❌ **Do NOT check** "Choose a license"
  - (You already have these files locally)
5. Click **"Create repository"**

### Step 2: Add Your Local Files to the Repository

**In GitHub Desktop:**

1. Click **File** → **Add Local Repository**
2. Click **"Choose..."** and navigate to your project folder:
   ```
   c:\Users\robys\OneDrive\Tugas Akhir (TA)\Thesis\Matkul\CFD\Tugas Besar
   ```
3. Click **"Select Folder"**
4. If GitHub Desktop says "This directory does not appear to be a Git repository":
   - Click **"Create a repository"**
   - The repository will be initialized

### Step 3: Publish to GitHub

1. In GitHub Desktop, you'll see **"Publish repository"** button in the top toolbar
2. Click **"Publish repository"**
3. A dialog appears with settings:
   - **Name**: `burgers-equation-cfd` (should match your GitHub repo)
   - **Description**: "2D Inviscid Burgers Equation CFD Solver"
   - **Keep this code private**: Check if you created a private repo
   - **Organization**: None (unless you're using a GitHub organization)
4. Click **"Publish repository"**

**GitHub Desktop will:**
- Create an initial commit with all your files
- Push everything to GitHub.com
- Connect your local repository to the remote

### Step 4: Verify Upload

1. Go to GitHub.com and navigate to your repository
2. You should see all your files:
   - `init_project.py`
   - `README.md`
   - `.gitignore`
   - `requirements.txt`
   - `GIT_WORKFLOW.md`
   - `References/` folder

## Option B: Create Repository Locally First

1. Open GitHub Desktop
2. Click **File** → **New Repository**
3. Fill in the form:
   - **Name**: `burgers-equation-cfd`
   - **Description**: "2D Inviscid Burgers Equation CFD Solver"
   - **Local path**: Choose the **parent folder** of your project
     - Example: `c:\Users\robys\OneDrive\Tugas Akhir (TA)\Thesis\Matkul\CFD\`
     - GitHub Desktop will create a `burgers-equation-cfd` folder inside
   - **Git ignore**: Select "Python"
   - **License**: None (or choose one if you prefer)
4. Click **"Create repository"**
5. Copy your existing files into the new folder
6. In GitHub Desktop, click **"Publish repository"** (top toolbar)
7. Choose visibility (private/public) and click **"Publish repository"**

---

## 2.4 Add Your Partner as Collaborator

**On GitHub.com:**

1. Go to your repository page
2. Click **"Settings"** tab (top of the page)
3. In the left sidebar, click **"Collaborators"** (under "Access")
4. You may need to confirm your GitHub password
5. Click **"Add people"** button
6. Enter your partner's GitHub username or email
7. Click **"Add [username] to this repository"**

**Your partner will receive:**
- An email invitation
- A notification on GitHub.com

**Your partner must:**
- Check their email or GitHub notifications
- Click **"Accept invitation"**
- They now have push/pull access to the repository

---

# 3. Cloning Repository (Person 2)

**This section is for the collaborator who will work on the project.**

## 3.1 Accept Invitation

1. Check your email for invitation from GitHub
2. Click **"View invitation"** in the email
   - Or go to GitHub.com → Check notifications (bell icon)
3. Click **"Accept invitation"**
4. You now have access to the repository

## 3.2 Clone Repository in GitHub Desktop

1. Open **GitHub Desktop**
2. Sign in to your GitHub account (if not already signed in)
3. Click **File** → **Clone Repository**
   - Or use the shortcut: `Ctrl+Shift+O` (Windows) / `Cmd+Shift+O` (Mac)

**In the Clone Repository dialog:**

### Option 1: From GitHub.com List
1. Click the **"GitHub.com"** tab
2. You should see `your-partner/burgers-equation-cfd` in the list
3. Select the repository
4. **Local path**: Choose where to save it on your computer
   - Example: `C:\Users\YourName\Projects\burgers-equation-cfd`
5. Click **"Clone"**

### Option 2: By URL
1. Click the **"URL"** tab
2. Enter the repository URL:
   ```
   https://github.com/your-partner-username/burgers-equation-cfd.git
   ```
3. **Local path**: Choose where to save it
4. Click **"Clone"**

**GitHub Desktop will:**
- Download all files from GitHub
- Set up the connection to the remote repository
- Open the repository in GitHub Desktop

## 3.3 Set Up Your Development Environment

### Open in File Explorer
1. In GitHub Desktop, click **Repository** → **Show in Explorer** (Windows)
   - Or **Show in Finder** (Mac)

### Install Python Dependencies
1. Open **Command Prompt** (Windows) or **Terminal** (Mac)
2. Navigate to the project folder:
   ```bash
   cd path\to\burgers-equation-cfd
   ```
3. Create virtual environment (recommended):
   ```bash
   python -m venv venv
   venv\Scripts\activate       # Windows
   source venv/bin/activate    # Mac/Linux
   ```
4. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

### Initialize Project Structure
```bash
python init_project.py
```

This creates the `src/`, `plots/`, and `data/` directories.

**You're now ready to collaborate!**

---

# 4. Branching Strategy

## 4.1 Why Use Branches?

**Branches allow you and your partner to work simultaneously without conflicts.**

Our strategy:
- **`main`**: Stable, tested code (both analytical + numerical merged)
- **`analytical-branch`**: Person 1 works on analytical solutions
- **`numerical-branch`**: Person 2 works on numerical methods

**Benefits:**
- No file conflicts (different folders: `src/analytical/` vs `src/numerical/`)
- Independent testing before merging
- Clear separation of work

## 4.2 Understanding Branches in GitHub Desktop

**Current Branch Indicator:**
- Look at the top center of GitHub Desktop
- You'll see a dropdown showing **"Current Branch: main"** with a branch icon
- Click it to see all branches and switch between them

## 4.3 Creating Your Branch

### Person 1: Create Analytical Branch

1. In GitHub Desktop, ensure you're on the **main** branch
   - Check the "Current Branch" dropdown at the top
2. Click **"Current Branch"** dropdown
3. Click **"New Branch"** button
4. **Name**: Type `analytical-branch`
5. **Create branch based on**: Select **main**
6. Click **"Create Branch"**

**You're now on `analytical-branch` (you'll see it in the dropdown)**

7. Click **"Publish branch"** in the top toolbar
   - This uploads the branch to GitHub so your partner can see it

### Person 2: Create Numerical Branch

Follow the same steps but name it `numerical-branch`:

1. Current Branch dropdown → **"New Branch"**
2. Name: `numerical-branch`
3. Create based on: **main**
4. Click **"Create Branch"**
5. Click **"Publish branch"**

## 4.4 Switching Between Branches

**To switch branches:**
1. Click the **"Current Branch"** dropdown
2. You'll see a list of all branches:
   - `main`
   - `analytical-branch`
   - `numerical-branch`
3. Click the branch you want to switch to

**GitHub Desktop will:**
- Update your working files to match that branch
- Show any uncommitted changes you have

**Visual indicator:**
- The current branch is shown with a checkmark ✓

## 4.5 Understanding Branch Visualization

**In GitHub Desktop:**
- Go to **History** tab (left sidebar)
- You'll see a graph showing:
  - Commit timeline
  - Where branches diverged
  - Which commits belong to which branch

**Branch colors:**
- Each branch gets a different color in the graph
- This helps you visualize when branches split and merge

---

# 5. Daily Workflow

## 5.1 Morning Routine: Pull Latest Changes

**Before starting work each day:**

1. Open GitHub Desktop
2. Check you're on your working branch (`analytical-branch` or `numerical-branch`)
   - See "Current Branch" dropdown at top
3. Click **"Fetch origin"** button in the top toolbar
   - This checks for new changes from GitHub
4. If changes are available, click **"Pull origin"**
   - This downloads the changes to your local computer

**Why do this?**
- Your partner may have made changes
- Main branch might have updates
- Prevents conflicts later

## 5.2 Making Changes to Your Code

1. Open your code editor (VS Code, PyCharm, etc.)
2. Work on your files in your designated folder:
   - **Person 1 (Analytical)**: Edit files in `src/analytical/`
   - **Person 2 (Numerical)**: Edit files in `src/numerical/`
3. Save your changes

**Example changes:**
- Create `src/analytical/solver.py`
- Implement analytical solution functions
- Create test scripts
- Generate plots

## 5.3 Viewing Changes in GitHub Desktop

**Switch to GitHub Desktop to see what changed:**

**Left Sidebar - "Changes" Tab:**
- Lists all modified, new, or deleted files
- Files with:
  - Green **+** icon = New file
  - Yellow **M** icon = Modified file
  - Red **-** icon = Deleted file

**Right Panel - Diff Viewer:**
- Click any file in the left sidebar
- The right panel shows the "diff" (differences):
  - **Green lines** = Added lines (with + symbol)
  - **Red lines** = Removed lines (with - symbol)
  - **White lines** = Unchanged context

**Example:**
```diff
  import numpy as np
  import matplotlib.pyplot as plt
+
+ def solve_burgers_analytical(x, y, t):
+     """Solve 2D Burgers equation analytically"""
+     u = np.sin(np.pi * x) * np.cos(np.pi * y)
+     v = -np.cos(np.pi * x) * np.sin(np.pi * y)
+     return u, v
```

**Checkboxes next to files:**
- Checked = Will be included in the commit
- Unchecked = Won't be included
- Use this to create focused commits (one feature at a time)

## 5.4 Committing Your Changes

**A commit is a snapshot of your work with a description.**

### Writing a Good Commit

**In GitHub Desktop:**

1. Review changed files in the left sidebar
2. Check/uncheck files to include in this commit
3. At the bottom left, you'll see two text fields:

**Summary Field (Required):**
- Short description (50 characters or less)
- Examples:
  ```
  Add analytical solver for Burgers equation
  Implement finite difference scheme
  Fix boundary condition bug
  Add visualization functions
  Update mesh generation algorithm
  ```

**Description Field (Optional):**
- Detailed explanation of what and why
- Use bullet points for multiple changes
- Example:
  ```
  - Implemented method of characteristics solver
  - Added initial condition functions
  - Created test cases for validation
  - Verified against known solutions
  ```

4. Click **"Commit to analytical-branch"** button (blue button at bottom)
   - Button text shows which branch you're committing to

**After committing:**
- Your changes are saved locally
- The "Changes" list clears
- The commit appears in the "History" tab

### Commit Best Practices

**Good commit messages:**
```
✓ Add analytical solution using method of characteristics
✓ Fix: Correct initial condition in numerical solver
✓ Refactor: Simplify mesh generation code
✓ Update: Improve plotting function with labels
```

**Bad commit messages:**
```
✗ Update
✗ Fix stuff
✗ Changes
✗ asdf
✗ Final version (really)
```

**How often to commit:**
- After completing a logical unit of work
- Every 30-60 minutes during active development
- Before taking a break
- Before trying something experimental

## 5.5 Pushing Changes to GitHub

**After committing locally, you need to push to GitHub:**

1. Look at the top toolbar in GitHub Desktop
2. You'll see a button that says:
   - **"Push origin"** with a number badge (e.g., ↑1, ↑3)
   - The number = how many commits haven't been pushed
3. Click **"Push origin"**

**GitHub Desktop will:**
- Upload your commits to GitHub.com
- Update the remote branch
- Sync your work so your partner can see it

**Important:**
- Push at least once per day (end of work session)
- Push before switching to a different computer
- Push before asking partner to review your code

### Understanding Push/Fetch/Pull

**Three operations:**

| Operation | Button | What it does |
|-----------|--------|--------------|
| **Fetch** | "Fetch origin" | Check if GitHub has new changes (doesn't modify your files) |
| **Pull** | "Pull origin" | Download changes from GitHub and merge into your local branch |
| **Push** | "Push origin" | Upload your local commits to GitHub |

**Visual indicator:**
- **↑ with number**: You have X commits to push
- **↓ with number**: GitHub has X commits you haven't pulled
- **↑↓ with numbers**: Both exist

## 5.6 Pulling Partner's Changes (Advanced)

**Scenario**: Your partner pushed changes to `main` that you want in your branch.

### Option 1: Merge main into your branch

1. Make sure all your changes are committed
2. Switch to **main** branch
   - Current Branch dropdown → Select **main**
3. Click **"Pull origin"** to get latest main
4. Switch back to your working branch
   - Current Branch dropdown → Select **analytical-branch**
5. Click **Branch** menu → **"Merge into current branch..."**
6. Select **main** from the list
7. Click **"Merge main into analytical-branch"**

### Option 2: Rebase (advanced users)

1. Branch → Rebase current branch
2. Select main
3. Resolve conflicts if any

---

# 6. Merging Changes

**When your work is complete and tested, merge it into `main`.**

## 6.1 Before Merging: Checklist

- [ ] All changes are committed
- [ ] Code is tested and working
- [ ] Your branch is pushed to GitHub
- [ ] You've communicated with your partner

## 6.2 Creating a Pull Request (Recommended Method)

**Pull Requests (PRs) allow code review before merging.**

### In GitHub Desktop:

1. Ensure you're on your working branch (`analytical-branch`)
2. Make sure everything is committed and pushed
3. Click **Branch** menu → **"Create Pull Request"**
   - Or click the **"Preview Pull Request"** button in the toolbar
   - Or press `Ctrl+R` (Windows) / `Cmd+R` (Mac)

**GitHub Desktop will:**
- Open your default browser
- Navigate to GitHub.com
- Show the "Open a pull request" page

### On GitHub.com:

**The pull request form shows:**

1. **Base branch**: `main` (where you want to merge TO)
2. **Compare branch**: `analytical-branch` (what you want to merge FROM)
3. **Title**: Auto-filled from your recent commit
   - Edit to be descriptive: "Add analytical solution implementation"
4. **Description**: Write a summary
   ```markdown
   ## What this PR does
   - Implements analytical solver using method of characteristics
   - Adds boundary condition functions
   - Creates visualization tools for analytical results

   ## Testing
   - Verified against exact solutions
   - Tested with various initial conditions
   - Generated validation plots

   ## Notes
   - Ready for review
   - No known issues
   ```

5. **Reviewers**: Click the gear icon and select your partner
6. Click **"Create pull request"** button

### Partner Review Process

**Your partner will:**
1. Get a notification email
2. Go to the Pull Request on GitHub.com
3. Click **"Files changed"** tab to review code
4. Leave comments or approve
5. Click **"Review changes"** → **"Approve"** → **"Submit review"**

### Merging the PR

**After approval:**

1. On the PR page, click **"Merge pull request"** button
2. Confirm by clicking **"Confirm merge"**
3. Optionally, click **"Delete branch"** (keeps repo clean)

**GitHub will:**
- Merge your code into `main`
- Close the pull request
- Update the repository

### After Merge: Update Your Local Repository

**In GitHub Desktop:**

1. Switch to **main** branch
2. Click **"Pull origin"** to get the merged code
3. Your `main` branch now includes your changes

## 6.3 Direct Merge (Alternative Method)

**If you don't want to use pull requests:**

1. In GitHub Desktop, switch to **main** branch
2. Click **"Pull origin"** to ensure main is up to date
3. Click **Branch** → **"Merge into current branch..."**
4. Select your working branch (e.g., `analytical-branch`)
5. Click **"Merge analytical-branch into main"**
6. Click **"Push origin"** to upload the merged main

---

# 7. Handling Merge Conflicts

**Conflicts happen when you and your partner edit the same lines in the same file.**

## 7.1 When Do Conflicts Occur?

**Common scenarios:**
- Both edited `README.md`
- Both modified shared utility functions
- Both changed the same configuration file

**How to avoid:**
- Work in separate folders (`src/analytical/` vs `src/numerical/`)
- Communicate before editing shared files
- Merge frequently to stay in sync

## 7.2 Identifying Conflicts in GitHub Desktop

**When you try to merge and there's a conflict:**

1. GitHub Desktop shows a warning:
   ```
   ⚠️ Resolve conflicts before merging
   ```
2. The conflicted files are listed with a warning icon
3. The **"Merge"** button is disabled

## 7.3 Resolving Conflicts

### Step 1: Open the Conflicted File

**In GitHub Desktop:**
1. Right-click the conflicted file
2. Select **"Open in [Your Editor]"**
   - Visual Studio Code
   - Sublime Text
   - Notepad
   - etc.

### Step 2: Understand Conflict Markers

**The file will contain special markers:**

```python
def initialize_mesh(nx, ny):
<<<<<<< HEAD (Current Change - main branch)
    x = np.linspace(0, 1, 100)
    y = np.linspace(0, 1, 100)
=======
    x = np.linspace(0, 1, 50)
    y = np.linspace(0, 1, 50)
>>>>>>> analytical-branch (Incoming Change)
    return x, y
```

**Markers explained:**
- `<<<<<<< HEAD`: Start of current branch's version
- `=======`: Separator between versions
- `>>>>>>> analytical-branch`: End of incoming branch's version

### Step 3: Choose the Correct Version

**You have three options:**

**Option 1: Keep current change (HEAD)**
```python
def initialize_mesh(nx, ny):
    x = np.linspace(0, 1, 100)
    y = np.linspace(0, 1, 100)
    return x, y
```

**Option 2: Keep incoming change**
```python
def initialize_mesh(nx, ny):
    x = np.linspace(0, 1, 50)
    y = np.linspace(0, 1, 50)
    return x, y
```

**Option 3: Combine both (manual edit)**
```python
def initialize_mesh(nx, ny):
    # Use 100 points for better resolution
    x = np.linspace(0, 1, 100)
    y = np.linspace(0, 1, 100)
    return x, y
```

**Important:**
- Delete all conflict markers (`<<<<<<<`, `=======`, `>>>>>>>`)
- Keep only the final code you want

### Step 4: Save and Mark as Resolved

1. Save the file in your editor
2. Return to GitHub Desktop
3. The file will now show a checkmark (resolved)
4. If multiple files have conflicts, repeat for each

### Step 5: Complete the Merge

1. All conflicts must be resolved
2. GitHub Desktop shows: "All conflicts resolved"
3. Write a commit message:
   ```
   Merge analytical-branch into main

   Resolved conflicts in mesh initialization
   ```
4. Click **"Commit merge"**
5. Click **"Push origin"** to upload

## 7.4 Using Visual Merge Tools

**GitHub Desktop can use external merge tools:**

1. **File** → **Options** (Windows) / **Preferences** (Mac)
2. Click **"Integrations"** tab
3. Choose your preferred merge tool:
   - Visual Studio Code
   - Atom
   - Sublime Text
   - etc.

**Then when conflicts occur:**
1. Right-click conflicted file
2. Select **"Open in external editor"**
3. Use the visual conflict resolution UI

---

# 8. Best Practices

## 8.1 Commit Frequency and Messages

**How often to commit:**
- ✓ After completing a feature
- ✓ Before taking a break
- ✓ Every 30-60 minutes during active work
- ✓ Before trying something risky
- ✗ Not every single line change
- ✗ Not once a week

**Commit message structure:**

```
Short summary (50 chars or less)

Detailed explanation if needed:
- What changed
- Why it changed
- Any side effects or notes

Fixes #123 (if fixing an issue)
```

**Examples:**

Good:
```
Add method of characteristics solver

- Implement MOC algorithm for Burgers equation
- Add initial condition handling
- Create validation test cases
```

Better:
```
Implement finite difference scheme for numerical solver

- Add forward-time, centered-space discretization
- Implement CFL condition check
- Add boundary condition handling for Dirichlet and Neumann
- Tested with known analytical solutions (error < 1e-6)
```

## 8.2 Branch Hygiene

**Keep branches clean:**

1. **Don't commit large binary files**
   - Plots (PNG, PDF) → Use `.gitignore`
   - Data files (CSV, HDF5) → Use `.gitignore`
   - Videos, archives → Use `.gitignore`

2. **Review changes before committing**
   - Look at the diff in GitHub Desktop
   - Uncheck files you don't want to commit
   - Ensure no debug code (e.g., `print()` statements)

3. **Delete merged branches**
   - After merging to main, delete feature branches
   - Keeps branch list clean
   - GitHub Desktop: Branch → Delete → Select branch

4. **Keep main stable**
   - Only merge tested, working code
   - Main should always be runnable

## 8.3 Communication with Partner

**Before major changes:**
- "I'm going to refactor the mesh generation code"
- "I'm adding a new dependency (Numba)"
- "I'm about to merge analytical-branch to main"

**Daily updates:**
- "Committed analytical solver, pushed to analytical-branch"
- "Found a bug in boundary conditions, fixing today"
- "Finished my part, ready for merge review"

**Use GitHub features:**
- Issues: Track bugs and features
- Discussions: Ask questions
- Pull Request comments: Code review feedback

## 8.4 File Organization

**Maintain the structure:**

```
burgers-equation-cfd/
├── src/
│   ├── __init__.py
│   ├── analytical/           ← Person 1 works here
│   │   ├── __init__.py
│   │   ├── solver.py
│   │   ├── initial_conditions.py
│   │   └── boundary_conditions.py
│   └── numerical/            ← Person 2 works here
│       ├── __init__.py
│       ├── solver.py
│       ├── mesh.py
│       └── time_stepping.py
├── plots/
│   ├── analytical/           ← .gitignored (not committed)
│   └── numerical/
├── data/                     ← .gitignored (not committed)
├── tests/                    ← Both can add tests
│   ├── test_analytical.py
│   └── test_numerical.py
├── README.md                 ← Shared documentation
├── requirements.txt          ← Coordinate dependency changes
└── .gitignore
```

**Minimizes conflicts because:**
- Each person has dedicated folders
- Shared files (README, requirements) are coordinated

## 8.5 Backup and Safety

**GitHub Desktop safety features:**

1. **Undo commit**: Repository → Undo most recent commit
2. **Revert commit**: Right-click commit in History → Revert
3. **Discard changes**: Right-click file → Discard changes

**Important:**
- Discarding changes is **permanent** - use carefully
- Undoing a pushed commit requires force-push (avoid)

**Backup strategy:**
- Pushed commits are safe (on GitHub servers)
- Unpushed commits are only on your computer
- Push daily to avoid losing work

---

## Part B: Command-Line Reference

# 9. Git Commands Quick Reference

**For users who prefer terminal or want to understand what's happening under the hood.**

## 9.1 Mapping GitHub Desktop Actions to Commands

| GitHub Desktop Action | Git Command |
|----------------------|-------------|
| Fetch origin | `git fetch origin` |
| Pull origin | `git pull origin main` |
| Push origin | `git push origin analytical-branch` |
| Create branch | `git checkout -b analytical-branch` |
| Switch branch | `git checkout main` |
| Commit | `git commit -m "message"` |
| Stage all files | `git add .` |
| Stage specific file | `git add src/analytical/solver.py` |
| View status | `git status` |
| View history | `git log --oneline -10` |
| View changes | `git diff` |
| Merge branch | `git merge analytical-branch` |
| Clone repository | `git clone https://github.com/user/repo.git` |

## 9.2 Initial Setup (Command-Line)

```bash
# Configure Git (first time only)
git config --global user.name "Your Name"
git config --global user.email "your-email@example.com"

# Initialize repository
cd path/to/project
git init

# Add remote
git remote add origin https://github.com/username/burgers-equation-cfd.git

# Initial commit and push
git add .
git commit -m "Initial project setup"
git push -u origin main
```

## 9.3 Daily Workflow (Command-Line)

```bash
# Morning: Update your branch
git checkout analytical-branch
git pull origin analytical-branch

# Make changes to files...
# Check what changed
git status
git diff

# Stage and commit
git add src/analytical/solver.py
git commit -m "Add analytical solver implementation"

# Push to GitHub
git push origin analytical-branch

# End of day: Ensure everything is pushed
git status
git push origin analytical-branch
```

## 9.4 Branching (Command-Line)

```bash
# Create and switch to new branch
git checkout -b analytical-branch

# Push branch to GitHub
git push -u origin analytical-branch

# List all branches
git branch -a

# Switch branches
git checkout main
git checkout analytical-branch

# Delete local branch
git branch -d analytical-branch

# Delete remote branch
git push origin --delete analytical-branch
```

## 9.5 Merging (Command-Line)

```bash
# Merge analytical-branch into main
git checkout main
git pull origin main
git merge analytical-branch
git push origin main

# Merge with pull request (using GitHub CLI)
gh pr create --title "Add analytical solver" --body "Description"
gh pr merge 1
```

## 9.6 Handling Conflicts (Command-Line)

```bash
# When merge has conflicts
git status                    # See conflicted files

# Edit conflicted files manually
# Remove <<<<<<, =======, >>>>>> markers

# After resolving
git add .
git commit -m "Resolve merge conflicts"
git push origin main
```

## 9.7 Useful Advanced Commands

```bash
# View commit history with graph
git log --graph --oneline --all

# See what changed in last commit
git show

# See changes since main branch
git diff main...analytical-branch

# Undo last commit (keep changes)
git reset --soft HEAD~1

# Undo last commit (discard changes)
git reset --hard HEAD~1

# Stash changes temporarily
git stash
git stash pop

# See who changed each line of a file
git blame src/analytical/solver.py

# Search commit messages
git log --grep="bug fix"
```

---

## Part C: Troubleshooting

# 10. Common Issues & Solutions

## 10.1 Authentication Issues

### Problem: "Authentication failed" when pushing

**Solution 1: Use Personal Access Token**
1. Go to GitHub.com → Settings → Developer settings → Personal access tokens
2. Click "Generate new token (classic)"
3. Select scopes: `repo` (full control)
4. Copy the token
5. In GitHub Desktop: File → Options → Accounts → Sign out → Sign in again
6. Use token as password

**Solution 2: Re-authenticate GitHub Desktop**
1. File → Options (Windows) / Preferences (Mac)
2. Accounts tab → Sign out
3. Sign in again

### Problem: "Permission denied" when accessing repository

**Solution:**
- Ensure you've accepted the collaboration invitation
- Check: GitHub.com → Repository → Settings → Collaborators
- Your username should be listed

## 10.2 Sync Issues

### Problem: "Your branch is behind origin/analytical-branch"

**Meaning:** GitHub has commits you don't have locally

**Solution in GitHub Desktop:**
1. Click "Pull origin" to download changes

**Solution in command-line:**
```bash
git pull origin analytical-branch
```

### Problem: "Your branch is ahead of origin/analytical-branch"

**Meaning:** You have local commits not pushed to GitHub

**Solution in GitHub Desktop:**
1. Click "Push origin" to upload commits

### Problem: "Failed to push - rejected"

**Meaning:** GitHub has changes you don't have

**Solution:**
1. Click "Pull origin" first (merge remote changes)
2. Resolve any conflicts
3. Then click "Push origin"

## 10.3 Branch Issues

### Problem: "I committed to the wrong branch"

**Solution in GitHub Desktop:**
1. Right-click the commit in History
2. Select "Undo commit"
3. Switch to the correct branch
4. Recommit the changes

**Solution in command-line:**
```bash
# Undo commit (keep changes)
git reset --soft HEAD~1

# Switch to correct branch
git checkout correct-branch

# Commit again
git commit -m "Your message"
```

### Problem: "I can't see my partner's branch"

**Solution:**
1. Click "Fetch origin" to update branch list
2. Current Branch dropdown should now show all branches
3. If still missing, ask partner to confirm they pushed their branch

### Problem: "Changes I made are gone after switching branches"

**Explanation:** Each branch has its own version of files

**Solution:**
- This is normal Git behavior
- Switch back to your branch to see your changes
- If uncommitted changes: Commit them before switching
- Or use "Stash" feature to temporarily save them

## 10.4 Merge and Conflict Issues

### Problem: "Merge conflict in README.md"

**Solution:**
1. Open the file in your editor
2. Find the conflict markers (`<<<<<<<`, `=======`, `>>>>>>>`)
3. Decide which version to keep
4. Delete the conflict markers
5. Save the file
6. In GitHub Desktop, file will show as resolved
7. Click "Commit merge"

### Problem: "Too many conflicts, want to start over"

**Solution in GitHub Desktop:**
1. Repository → Abort merge
2. Try merging again, or ask partner for help

**Solution in command-line:**
```bash
git merge --abort
```

### Problem: "I merged the wrong branch"

**Solution:**
```bash
# Undo merge (if not yet pushed)
git reset --hard HEAD~1

# If already pushed (dangerous - coordinate with partner)
git revert -m 1 HEAD
git push origin main
```

## 10.5 File and Change Issues

### Problem: "Changed files aren't showing in GitHub Desktop"

**Checklist:**
- [ ] Are you on the correct branch?
- [ ] Did you save the file in your editor?
- [ ] Is the file inside the repository folder?
- [ ] Is the file ignored in `.gitignore`?

**Solution:**
1. Repository → Refresh (F5)
2. Check `.gitignore` - remove the file pattern if it shouldn't be ignored

### Problem: "I want to uncommit a file but keep the changes"

**Solution in GitHub Desktop:**
1. Repository → Undo most recent commit
2. The changes reappear in the Changes tab
3. Uncheck the file you don't want
4. Recommit

### Problem: "Large file rejected by GitHub"

**Error:** "File is larger than 100 MB"

**Solution:**
1. Add the file to `.gitignore`
2. Use Git LFS (Large File Storage):
   ```bash
   git lfs install
   git lfs track "*.csv"
   git add .gitattributes
   ```
3. Or store large files elsewhere (Google Drive, OneDrive)

## 10.6 General Troubleshooting

### Problem: "GitHub Desktop is frozen or not responding"

**Solution:**
1. Close GitHub Desktop
2. Reopen it
3. If still frozen: Restart computer
4. Last resort: Reinstall GitHub Desktop (your repositories are safe)

### Problem: "I don't understand what went wrong"

**Debugging steps:**
1. **Check status**: What does GitHub Desktop show?
2. **Check history**: Repository → History - what's the last commit?
3. **Check branch**: What branch are you on?
4. **Ask partner**: Did they make changes recently?
5. **Google the error**: Copy exact error message to Google
6. **Ask for help**: GitHub Community Forum, Stack Overflow

### Problem: "Everything is broken, I give up"

**Nuclear option (fresh start):**

1. **Backup your work:**
   - Copy all your code files to a safe location outside the repository

2. **Delete the local repository:**
   - In GitHub Desktop: Repository → Remove → Delete

3. **Clone again:**
   - File → Clone Repository
   - Select your repository

4. **Restore your changes:**
   - Copy your saved code files back
   - Commit and push

**This loses local commit history but preserves your code.**

---

## Emergency Commands

### Undo Last Commit (Keep Changes)
```bash
git reset --soft HEAD~1
```

### Undo Last Commit (Discard Changes)
```bash
git reset --hard HEAD~1
```

### Discard All Local Changes
```bash
git reset --hard HEAD
```

### See What You Recently Did
```bash
git reflog
```

### Recover Deleted Commit
```bash
git reflog              # Find commit hash
git checkout <hash>     # View it
git cherry-pick <hash>  # Apply it to current branch
```

---

## Getting Help

### In GitHub Desktop
- Help menu → Show Logs
- Help menu → Report Issue

### Online Resources
- [GitHub Desktop Documentation](https://docs.github.com/en/desktop)
- [Git Documentation](https://git-scm.com/doc)
- [GitHub Community Forum](https://github.community/)

### Command Help
```bash
git help <command>
# Example:
git help merge
git help commit
```

---

## Summary: Quick Daily Checklist

**Every morning:**
- [ ] Open GitHub Desktop
- [ ] Fetch origin
- [ ] Pull origin (if changes available)

**During work:**
- [ ] Make changes to code
- [ ] Review changes in GitHub Desktop
- [ ] Write good commit message
- [ ] Commit
- [ ] Push (at least once a day)

**Before merging:**
- [ ] Test your code
- [ ] Communicate with partner
- [ ] Create pull request
- [ ] Review and merge

**After merging:**
- [ ] Switch to main
- [ ] Pull origin
- [ ] Delete merged branch (optional)

---

**You're now ready to collaborate effectively on your CFD project!**

For questions or issues not covered here, check the troubleshooting section or ask your partner for help. Git is powerful and forgiving - you can almost always recover from mistakes.
