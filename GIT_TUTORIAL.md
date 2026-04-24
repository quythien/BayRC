# Git Tutorial for BayRC

A practical reference covering everything you need for daily work on this project.

---

## 1. The mental model

Git tracks **snapshots** of your project over time. Every `git commit` creates a permanent
snapshot you can always return to. Think of it as an unlimited undo history that also lets
multiple people work in parallel.

Three zones to know:

```
Working directory  →  Staging area  →  Repository (commits)
  (your files)       (git add)         (git commit)
```

---

## 2. Checking what you have

```bash
# What files have changed?
git status

# What exactly changed in each file?
git diff

# What changed in staged files (already git add'd)?
git diff --staged

# History of commits
git log --oneline -10

# Who changed what line in a file?
git blame R/mcmc.R
```

---

## 3. Making a commit

Always stage specific files — never `git add .` blindly (risk committing `.RData`, secrets, etc.).

```bash
# Stage specific files
git add R/mcmc.R R/internal.R DESCRIPTION

# Stage parts of a file interactively (pick individual hunks)
git add -p R/internal.R

# Commit with a message
git commit -m "Fix RJMCMC acceptance ratio: log.phi.prior = log(1/P)"

# Check the commit was created
git log --oneline -3
```

**Good commit messages** follow this pattern:
```
Short summary (50 chars max, imperative mood)

Optional body explaining WHY, not WHAT. The diff shows what changed;
the message explains the reasoning. Reference issue numbers if relevant.
```

---

## 4. Branches — working on features safely

```bash
# See all branches
git branch

# Create and switch to a new branch
git checkout -b fix/convergence-diagnostics

# Switch between branches
git checkout main
git checkout fix/convergence-diagnostics

# Merge a branch back to main (fast-forward when possible)
git checkout main
git merge fix/convergence-diagnostics

# Delete a branch after merging
git branch -d fix/convergence-diagnostics
```

**Branch naming convention for this project:**
- `fix/short-description` — bug fixes
- `feat/short-description` — new features
- `paper/short-description` — paper-only changes

---

## 5. Viewing history and comparing

```bash
# See what changed between two commits
git diff HEAD~3 HEAD -- R/mcmc.R

# See a specific commit
git show 5ae6edd

# Search commit messages
git log --oneline --grep="MCMC"

# See when a line was introduced
git log -S "set.seed" R/mcmc.R

# Graph of branches
git log --oneline --graph --all
```

---

## 6. Undoing things

```bash
# Unstage a file (keep changes in working directory)
git restore --staged R/mcmc.R

# Discard working directory changes in a file (IRREVERSIBLE)
git restore R/mcmc.R

# Undo the last commit but keep changes staged
git reset --soft HEAD~1

# Create a new commit that reverses a previous one (safe for shared history)
git revert 5ae6edd
```

**Rule of thumb:** `git restore` and `git reset --soft` are safe.
`git reset --hard` and `git push --force` are destructive — ask first.

---

## 7. Stashing — saving work without committing

```bash
# Save current changes temporarily (e.g., need to switch branch quickly)
git stash

# List stashes
git stash list

# Bring stash back
git stash pop

# Name a stash
git stash push -m "half-finished ESS diagnostics"
```

---

## 8. Working with remotes

```bash
# See configured remotes
git remote -v

# Push your commits to the remote
git push origin main

# Pull latest changes from remote
git pull origin main

# Push a new branch to remote
git push -u origin fix/convergence-diagnostics
```

---

## 9. The `.gitignore` in this project

The `.gitignore` uses a **whitelist** approach: everything is ignored by default (`*`)
and specific paths are explicitly allowed (`!R/`, `!DESCRIPTION`, etc.).

```
*                   ← ignore everything by default
!R/                 ← but allow R/ source
!R/**               ← and everything inside it
!DESCRIPTION        ← allow DESCRIPTION
!NEWS.md            ← allow NEWS.md
paper/              ← paper is LOCAL ONLY — not tracked in git
```

**To add a new file type that should be tracked:**
```bash
echo "!vignettes/my_new_file.Rmd" >> .gitignore
git add .gitignore vignettes/my_new_file.Rmd
```

**To check why a file is being ignored:**
```bash
git check-ignore -v paper/BayRC.tex
```

---

## 10. Typical daily workflow for this project

```bash
# 1. Start of day: make sure you have the latest
git pull origin main

# 2. Create a branch for what you're working on
git checkout -b fix/phase-shift-threshold

# 3. Make changes, check what changed
git status
git diff R/internal.R

# 4. Stage and commit
git add R/internal.R paper/BayRC.tex
git commit -m "Fix phase_infer default shift: 4h -> 2h (match paper ±2h claim)"

# 5. When done, merge back and push
git checkout main
git merge fix/phase-shift-threshold
git push origin main

# 6. Clean up branch
git branch -d fix/phase-shift-threshold
```

---

## 11. Common scenarios in R package development

**After running `devtools::document()` (regenerates NAMESPACE and man/):**
```bash
git add NAMESPACE man/
git commit -m "Regenerate documentation after API changes"
```

**After bumping version:**
```bash
# Edit DESCRIPTION Version: field, then:
git add DESCRIPTION NEWS.md
git commit -m "Bump version to 0.2.1"
git tag v0.2.1          # tag the release
git push origin v0.2.1  # push the tag
```

**Checking what changed since last version:**
```bash
git log v0.1.0..HEAD --oneline
git diff v0.1.0..HEAD -- R/mcmc.R
```

---

## 12. Quick reference card

| What you want | Command |
|---|---|
| See what changed | `git status` / `git diff` |
| Stage a file | `git add filename` |
| Commit | `git commit -m "message"` |
| See history | `git log --oneline -10` |
| New branch | `git checkout -b branch-name` |
| Switch branch | `git checkout main` |
| Merge branch | `git merge branch-name` |
| Undo staged | `git restore --staged filename` |
| Undo file changes | `git restore filename` |
| Save work temporarily | `git stash` / `git stash pop` |
| Push to remote | `git push origin main` |
| Pull from remote | `git pull origin main` |
| See a specific commit | `git show <hash>` |
| Search history | `git log --grep="keyword"` |
