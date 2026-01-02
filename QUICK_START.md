# Quick Start Checklist

**Fast-track guide for setting up the 2D Inviscid Burgers Equation CFD project**

---

## Person 1: Repository Creator

### Initial Setup (One-time)

- [ ] **Install GitHub Desktop**
  - Download from [desktop.github.com](https://desktop.github.com)
  - Sign in with your GitHub account
  - Configure name and email

- [ ] **Run project initialization script**
  ```bash
  python init_project.py
  ```

- [ ] **Create repository on GitHub.com**
  - Go to GitHub.com → New repository
  - Name: `burgers-equation-cfd`
  - Visibility: Private or Public
  - Do NOT initialize with README (you have one locally)

- [ ] **Publish project to GitHub**
  - In GitHub Desktop: File → Add Local Repository
  - Select your project folder
  - Click "Publish repository"
  - Verify files appear on GitHub.com

- [ ] **Invite your partner**
  - GitHub.com → Repository → Settings → Collaborators
  - Click "Add people"
  - Enter partner's GitHub username or email
  - Send invitation

### Development Setup

- [ ] **Create your working branch**
  - In GitHub Desktop: Current Branch → New Branch
  - Name: `analytical-branch`
  - Based on: `main`
  - Click "Publish branch"

- [ ] **Make your first commit**
  - Create `src/analytical/solver.py` (or similar file)
  - In GitHub Desktop, review changes
  - Write commit message: "Initial analytical solver structure"
  - Click "Commit to analytical-branch"
  - Click "Push origin"

- [ ] **Verify everything is working**
  - Check GitHub.com to see your branch and commit
  - Share repository link with partner

---

## Person 2: Collaborator

### Initial Setup (One-time)

- [ ] **Install GitHub Desktop**
  - Download from [desktop.github.com](https://desktop.github.com)
  - Sign in with your GitHub account
  - Configure name and email

- [ ] **Accept collaboration invitation**
  - Check email or GitHub.com notifications
  - Click "Accept invitation"
  - Confirm you have access to repository

- [ ] **Clone the repository**
  - In GitHub Desktop: File → Clone Repository
  - Select `your-partner/burgers-equation-cfd` from list
  - Choose local path (where to save on your computer)
  - Click "Clone"

- [ ] **Set up Python environment**
  - Open Command Prompt / Terminal
  - Navigate to project folder
  ```bash
  cd path/to/burgers-equation-cfd

  # Create virtual environment
  python -m venv venv

  # Activate virtual environment
  venv\Scripts\activate          # Windows
  source venv/bin/activate       # Mac/Linux

  # Install dependencies
  pip install -r requirements.txt
  ```

- [ ] **Initialize project structure**
  ```bash
  python init_project.py
  ```

### Development Setup

- [ ] **Create your working branch**
  - In GitHub Desktop: Current Branch → New Branch
  - Name: `numerical-branch`
  - Based on: `main`
  - Click "Publish branch"

- [ ] **Make your first commit**
  - Create `src/numerical/solver.py` (or similar file)
  - In GitHub Desktop, review changes
  - Write commit message: "Initial numerical solver structure"
  - Click "Commit to numerical-branch"
  - Click "Push origin"

- [ ] **Verify everything is working**
  - Check GitHub.com to see your branch and commit
  - Confirm you can see partner's branch

---

## Both Partners: Daily Workflow

### Every Morning

- [ ] Open GitHub Desktop
- [ ] Click "Fetch origin" (check for updates)
- [ ] Click "Pull origin" (if updates available)
- [ ] Verify you're on your working branch

### During Development

- [ ] Make code changes in your designated folder:
  - Person 1: `src/analytical/`
  - Person 2: `src/numerical/`
- [ ] Review changes in GitHub Desktop
- [ ] Commit with meaningful messages
- [ ] Push at least once per day

### Before Merging to Main

- [ ] Test your code thoroughly
- [ ] Commit all changes
- [ ] Push to your branch
- [ ] Communicate with partner
- [ ] Create Pull Request on GitHub.com
- [ ] Wait for partner review
- [ ] Merge after approval
- [ ] Pull latest main branch
- [ ] Continue development

---

## Verification Tests

### Test 1: Can you see each other's branches?

**Both partners:**
- [ ] Open GitHub Desktop
- [ ] Click "Fetch origin"
- [ ] Click "Current Branch" dropdown
- [ ] Verify you can see:
  - `main`
  - `analytical-branch`
  - `numerical-branch`

### Test 2: Can you push and pull?

**Person 1:**
- [ ] Create a test file in `src/analytical/test.txt`
- [ ] Commit and push
- [ ] Verify on GitHub.com

**Person 2:**
- [ ] Switch to `analytical-branch` (to view partner's work)
- [ ] Click "Pull origin"
- [ ] Verify test file appears locally
- [ ] Switch back to `numerical-branch`

### Test 3: Can you merge?

**After both have made some commits:**
- [ ] Both push latest changes
- [ ] Person 1 creates PR: `analytical-branch` → `main`
- [ ] Person 2 reviews and approves
- [ ] Person 1 merges PR
- [ ] Both pull latest `main` branch
- [ ] Verify merged code works

---

## Common First-Time Issues

### "I can't find GitHub Desktop"

**Solution:**
- Windows: Check Start Menu → Search "GitHub Desktop"
- Mac: Check Applications folder

### "Authentication failed"

**Solution:**
- Sign out of GitHub Desktop
- Sign back in
- May need Personal Access Token (see [GIT_WORKFLOW.md](GIT_WORKFLOW.md#101-authentication-issues))

### "I don't see my partner's branch"

**Solution:**
- Click "Fetch origin" to update branch list
- Ask partner to confirm they clicked "Publish branch"

### "Python packages won't install"

**Solution:**
- Ensure virtual environment is activated
- Check Python version: `python --version` (should be 3.8+)
- Try: `pip install --upgrade pip` then retry

### "init_project.py doesn't work"

**Solution:**
- Ensure you're in the project root directory
- Check file exists: `ls init_project.py` or `dir init_project.py`
- Run with: `python init_project.py`

### "I committed to the wrong branch"

**Solution:**
- Right-click the commit in History
- Select "Undo commit"
- Switch to correct branch
- Commit again

---

## Project Structure After Setup

```
burgers-equation-cfd/
├── src/
│   ├── __init__.py
│   ├── analytical/              ← Person 1 works here
│   │   ├── __init__.py
│   │   └── solver.py
│   └── numerical/               ← Person 2 works here
│       ├── __init__.py
│       └── solver.py
├── plots/
│   ├── analytical/
│   └── numerical/
├── data/
├── tests/
├── References/                  ← Keep reference PDFs here
├── README.md
├── GIT_WORKFLOW.md             ← Detailed GitHub Desktop guide
├── QUICK_START.md              ← This file
├── requirements.txt
├── .gitignore
└── init_project.py
```

---

## Next Steps

Once both partners have completed all checklists:

1. **Read [GIT_WORKFLOW.md](GIT_WORKFLOW.md)** for detailed GitHub Desktop instructions
2. **Review [README.md](README.md)** for project overview and usage examples
3. **Start implementing**:
   - Person 1: Analytical solutions in `src/analytical/`
   - Person 2: Numerical methods in `src/numerical/`
4. **Communicate regularly** via:
   - GitHub Issues (for bugs and features)
   - Pull Request comments (for code review)
   - Direct messaging (for coordination)

---

## Quick Reference Commands

### Python Virtual Environment

```bash
# Activate
venv\Scripts\activate              # Windows
source venv/bin/activate           # Mac/Linux

# Deactivate
deactivate
```

### GitHub Desktop Shortcuts

- `Ctrl+Shift+O` / `Cmd+Shift+O`: Clone Repository
- `Ctrl+Shift+N` / `Cmd+Shift+N`: New Repository
- `Ctrl+T` / `Cmd+T`: Switch Branch
- `Ctrl+R` / `Cmd+R`: Create Pull Request
- `F5`: Refresh

### Git Commands (Alternative to GitHub Desktop)

```bash
git status                    # Check current state
git add .                     # Stage all changes
git commit -m "message"       # Commit
git push origin branch-name   # Push
git pull origin branch-name   # Pull
```

---

## Getting Help

**Detailed guides:**
- [GIT_WORKFLOW.md](GIT_WORKFLOW.md) - Complete GitHub Desktop workflow
- [README.md](README.md) - Project documentation

**Online resources:**
- [GitHub Desktop Docs](https://docs.github.com/en/desktop)
- [Git Basics](https://git-scm.com/book/en/v2/Getting-Started-Git-Basics)

**Troubleshooting:**
- See [GIT_WORKFLOW.md - Part C: Troubleshooting](GIT_WORKFLOW.md#10-common-issues--solutions)

---

**You're ready to collaborate! Good luck with your CFD project!** 🚀
