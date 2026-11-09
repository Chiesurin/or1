---
name: tableau-solver
description: >-
  Solves Simplex and Two-Phase tableaus step-by-step using a deterministic Python engine.
  Use when the user provides a Simplex tableau (via CSV, file path, or inline text in the prompt)
  or asks to solve a Linear Program through tableaus, expliciting all iterations, row operations,
  ratio tests, Bland's rule, and generating solution markdown files without external LLM APIs.
---

# Simplex Tableau Solver Skill

This skill enables an LLM agent to reliably solve any Simplex tableau or Linear Program by executing a pure-Python exact calculation engine.

> [!IMPORTANT]
> **No Mental Math / No LLM Arithmetic**: LLMs should never attempt to calculate matrix pivots, fractional quotients, or Gaussian row operations internally. Always invoke the local Python solver script [solve_tableau.py](./scripts/solve_tableau.py) to perform exact calculations and format the output.
> **Python Optimization Library (SciPy)**: Integrates directly with SciPy's linear programming engine (`scipy.optimize.linprog`, method `highs`). If SciPy is installed, it computes and cross-verifies the optimal solution. If SciPy is not yet installed, install via `pip3 install --user scipy numpy`.
> **Zero LLM API**: All calculations are executed locally on the user's machine. No OpenAI, Anthropic, or external API keys are required.

---

## Capabilities & Mathematical Conventions

1. **Universal Minimization ($\min z$)**:
   - The engine operates under universal minimization.
   - **Row 0 Identity**: $-z = \text{RHS}_0 \implies z^* = -\text{RHS}_0$.
2. **Exact Rational Arithmetic**:
   - Every number is represented as an exact irreducible fraction (`Fraction(a, b)`). Zero floating-point rounding errors.
3. **Bland's Smallest-Subscript Anti-Cycling Rule**:
   - **Entering Variable**: Smallest variable index $j$ among candidate columns with strictly negative reduced costs ($c_j < 0$).
   - **Leaving Variable**: Smallest variable index $i$ among candidate rows achieving the minimum ratio $\theta = \text{RHS}_i / a_{ij}$ ($a_{ij} > 0$).
4. **Two-Phase Simplex**:
   - Automatic when negative RHS or equality/$\ge$ constraints exist, or explicitly requested via `--two-phase`.
   - Phase 1 drives artificial variables $y_1, \dots, y_m$ to 0. Phase 2 restores and optimizes the original objective.
5. **Dual Simplex Algorithm (`--dual-simplex`)**:
   - Solves dual-feasible tableaus (all reduced costs in Row 0 satisfy $c_j \ge 0$) with negative RHS entries ($b_i < 0$).
   - Pivots by selecting leaving variable $R_l$ with most negative RHS, then performing the dual minimum ratio test $\min_{j: a_{lj} < 0} (c_j / |a_{lj}|)$ to pick the entering variable.
6. **Dual Problem Formulation & Duality Theorems (`--show-dual`)**:
   - Automatically formulates the mathematical dual Linear Program from primal LP formulations (min $\leftrightarrow$ max, transposed constraint matrix $A^T$, dual variable bounds according to primal constraint types).
   - Computes exact shadow prices (dual variables $y^* = c_B^T B^{-1}$) from Row 0 reduced costs and basis inversion.
   - Verifies **Strong Duality Theorem** ($z^*_{\text{primal}} == w^*_{\text{dual}}$).
   - Verifies **Complementary Slackness Theorem** ($y_i^* \cdot s_i = 0$ for primal slacks, $x_j^* \cdot \bar{c}_j = 0$ for dual slacks).
7. **Full Calculation Trace**:
   - Explicitly records every step:
     - Optimality test (Row 0 inspection)
     - Bland candidate selection
     - Ratio test table ($\theta$ calculations)
     - Pivot element $(p, q)$
     - Elementary row operations ($R_p \leftarrow \frac{1}{a_{pq}} R_p$, $R_i \leftarrow R_i - a_{iq} R_p$)
     - Resulting tableau and basic solution

---

## Workflow for LLM Agents

When a user requests to solve a tableau or LP:

### Step 1: Identify Input Source

- **Case A: CSV file path provided** (e.g. `tableau.csv`, `csv/tableau.csv`):
  Run:
  ```bash
  python3 skills/tableau-solver/scripts/solve_tableau.py --csv <path_to_csv>
  ```
- **Case B: Tableau text pasted inside the prompt** (comma, semicolon, or whitespace separated):
  Run:
  ```bash
  python3 skills/tableau-solver/scripts/solve_tableau.py --raw "<pasted_tableau_text>"
  ```
- **Case C: Standard Linear Program text** (e.g. `min 4*x1 + x2 s.t. ...`):
  Run:
  ```bash
  python3 skills/tableau-solver/scripts/solve_tableau.py --lp "<lp_text>"
  ```
- **Case D: Two-Phase Simplex requested**:
  Add `--two-phase` flag to any command above.

### Step 2: Generate Solution Markdown File

When the user asks for a `.md` file containing the solution (or when creating an artifact for the user):
Add `--output-md <output_path.md>` to the command:
```bash
python3 skills/tableau-solver/scripts/solve_tableau.py --csv <file.csv> --output-md solution.md
```
Or for inline prompts:
```bash
python3 skills/tableau-solver/scripts/solve_tableau.py --raw "<tableau_text>" --output-md solution.md
```

### Step 3: Present the Results to the User

Present the output cleanly to the user:
1. **Summary**: State problem dimension, initial basic feasible solution, and phase mode.
2. **Walkthrough of Every Iteration**:
   - Show the tableau table.
   - State entering variable chosen by Bland's rule.
   - Show the ratio test table with all evaluated fractions.
   - State leaving variable and pivot element.
   - Explicitly write the elementary row operation formulas ($R_i \leftarrow R_i - c \cdot R_p$).
3. **Final Result**:
   - Status (`OPTIMAL`, `INFEASIBLE`, or `UNBOUNDED`).
   - Optimal objective value: $\min z^* = -\text{RHS}_0$.
   - Optimal variable vector: $x^* = (x_1^*, x_2^*, \dots)$.
   - Clickable link to the generated `.md` file.

---

## CLI Options Reference

| Flag | Description |
| :--- | :--- |
| `--csv <path>` | Load tableau from CSV file (auto-searches `csv/`, workspace root, `examples/`). |
| `--raw "<text>"` | Parse raw multiline tableau text directly from command line / prompt. |
| `--stdin` | Read raw tableau text from standard input pipe. |
| `--lp "<text>"` | Parse standard Linear Program formulation text (min/max, constraints). |
| `--two-phase` | Force Two-Phase Simplex even if an initial slack basis exists. |
| `--dual-simplex` | Solve using the Dual Simplex algorithm (requires dual-feasible tableau $c_j \ge 0$). |
| `--show-dual` | Formulate and display the mathematical dual Linear Program. |
| `--output-md <path>` | Export complete, publication-grade solution report to a Markdown `.md` file. |
| `--format {text,markdown,json}` | Output format in terminal stdout (default: `text`). |
| `--verify-scipy` | Cross-validate optimal objective value with `scipy.optimize.linprog` if available. |

---

## Directory Reference

- [solve_tableau.py](./scripts/solve_tableau.py): Core standalone Python solver engine.
- [simplex_cheatsheet.md](./references/simplex_cheatsheet.md): Mathematical derivations, Row 0 identity, Bland's rule.
- [input_formats.md](./references/input_formats.md): Accepted CSV and prompt formatting rules.
- [sample_tableau.csv](./examples/sample_tableau.csv): Example standard CSV tableau.
- [sample_solution.md](./examples/sample_solution.md): Example generated solution report.
