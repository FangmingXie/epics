# CLAUDE.md

## Project

Epigenomics analysis. General-purpose pipeline for processing and analyzing epigenomic data.

## Repository Structure

```
epics/
├── CLAUDE.md
├── README.md
├── LICENSE
├── .gitignore
├── scripts/                                     # pipeline shell/python scripts
├── notebooks/                                   # exploratory Jupyter notebooks
├── plan/                                        # finalized implementation plans (see automation rules)
├── env/                                         # saved conda/pip package lists
├── logs/                                        # run logs (gitignored)
└── data/                                        # input data and results (gitignored)
```

## Git Branches

- `main`: active development branch

## Git Configuration

- user.name: FangmingXie
- user.email: fmxie1993@gmail.com

## Environment

- Use this conda env to run this project: `epics`

## .gitignore Notes

- `data/`, `logs/`, and `*.log` are gitignored

## Coding Styles

- Define all file paths (input and output files) in the beginning of each script as much as possible. Capitalize the variables that store these file paths.

**Simplify Relentlessly**: Remove complexity aggressively - the simplest design that works is usually best

#### Fail-Fast, No Fallbacks
- **No Silent Fallbacks**: Code must fail immediately when expected conditions aren't met.
- **Explicit Error Messages**: When something goes wrong, stop execution with clear error messages explaining what failed and what was expected.

### ⚠️ **IMPORTANT: Rewrite Project - Breaking Changes Encouraged**

**This project is under active development**, not a stable codebase with external dependencies. This means:

- **Breaking changes are encouraged** when they follow best practices
- **No backward compatibility constraints** - optimize for clean design
- **Clean organization** - each script has a single, clear purpose

## Claude Code Automation Rules
- When operating in Plan Mode, ALWAYS save the finalized implementation plan as a distinct markdown file under the `plan/` folder before concluding the turn.
- Never execute modifications while Plan Mode is toggled active.

## ⚠️ **IMPORTANT: Installing a package
- do not attempt to install new packages without asking permission explicitly. 
- before installing anything new, save a copy of the current package list under `env/`. 
- always try using conda first, and pip later. 
