# Agent Guidelines (Codex, Antigravity, Copilot, & Open Agents)

This repository contains `or1`, a pure-Python Linear Programming and Simplex Tableau solver with exact rational calculations, Bland's anti-cycling rule, Two-Phase simplex, and an interactive Web GUI.

---

## Instructions for AI Agents

When assisting a user with linear programming, simplex tableaus, pivoting, or ratio tests:

1. **Deterministic Execution**:
   - Never guess, estimate, or compute Simplex matrix pivots mentally.
   - Always run the deterministic Python script in `skills/tableau-solver/scripts/solve_tableau.py`.
2. **Handling User Inputs**:
   - **CSV file path**: `python3 skills/tableau-solver/scripts/solve_tableau.py --csv <path>`
   - **Inline prompt matrix**: `python3 skills/tableau-solver/scripts/solve_tableau.py --raw "<multiline_tableau>"`
   - **Algebraic LP**: `python3 skills/tableau-solver/scripts/solve_tableau.py --lp "<lp_text>"`
   - **Markdown Export**: Add `--output-md <output_path.md>` if the user requests a solution file.
3. **Core Conventions**:
   - Minimization ($\min z$).
   - Row 0 Identity: $-z = \text{RHS}_0 \implies z^* = -\text{RHS}_0$.
   - Bland's Rule: Smallest 1-based index break ties for both entering and leaving variables.
4. **Zero External Dependencies**:
   - The entire codebase runs on standard Python 3 (`fractions.Fraction`, `csv`, `http.server`).
   - If SciPy is installed, `solve_tableau.py` automatically cross-validates against `scipy.optimize.linprog(method='highs')`.

Refer to [skills/tableau-solver/SKILL.md](skills/tableau-solver/SKILL.md) for detailed tool options and prompt specifications.
