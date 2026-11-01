# Simplex Method with Bland's Rule in Python (CLI & Web GUI)

A complete, exact implementation of the **Simplex Algorithm** in pure Python using **Bland's Smallest-Subscript Rule** to prevent cycling and guarantee finite convergence to optimality.

Includes both a **Terminal CLI** and a modern, responsive, dependency-free **Web GUI Visualizer** served via Python's standard `http.server`.

---

## 🤖 AI / Agent Quick Reference Summary

> **Context for AI Agents & Automated Tools:**
> - **Language & Runtime:** Python 3 (3.8+). Zero external third-party dependencies (`pip install` not required; uses only standard library: `fractions`, `http.server`, `urllib`, `json`, `csv`, `glob`, `socket`, `webbrowser`).
> - **Arithmetic Engine:** Exact rational arithmetic via `fractions.Fraction`. All calculations (pivot operations, reduced costs, ratio tests) are exact with no floating-point rounding errors.
> - **Primary Source File:** [`source/gui.py`](source/gui.py) contains both the Simplex engine, CLI runner, and the self-contained HTML/CSS/JS single-page web GUI.
> - **CSV Tableaus Directory:** [`csv/`](csv/) holds standard CSV-formatted tableaus.
> - **Optimization Direction:** The solver universally minimizes $z$ ($\min z$). Maximization problems $\max c^T x$ are converted via $\min -c^T x$.
> - **Row 0 Convention:** Tableau Row 0 is labeled `-z` (or `-w` during Phase 1). The internal equation is $-z + \sum \bar{c}_j x_j = \text{RHS}_0$, meaning the true objective value is $z = -\text{RHS}_0$.
> - **Tableau Layout:** Column headers are ordered as `Basis, RHS, x1, x2, ...` with **no dedicated $z$ column** and RHS on the left.
> - **Anti-Cycling Guarantee:** Employs Bland's smallest-index rule for both entering and leaving variable selection to avoid degenerate cycling (verified against Beale's 1955 cycling problem).

---

## 📁 Repository Structure

```
.
├── .gitignore              # Ignores .DS_Store, __pycache__, bytecode
├── LICENSE                 # Project license
├── README.md               # Complete documentation and AI reference guide
├── csv/                    # Simplex initial tableau CSV files
│   ├── example_standard.csv # Standard 2D textbook LP problem
│   └── tableau.csv         # Beale's classic cycling problem (1955)
└── source/                 # Python source files
    └── gui.py              # Core solver engine, REST API, CLI, and Web GUI
```

---

## 📐 Mathematical Formulation & Conventions

### 1. Optimization Goal (Universal Minimization)
All problems are formulated in standard minimization form:
$$\min_{\mathbf{x}} \; \mathbf{c}^T \mathbf{x} \quad \text{subject to} \quad A\mathbf{x} = \mathbf{b}, \quad \mathbf{x} \ge 0$$
- If a maximization objective is provided, it is negated: $\max f(\mathbf{x}) \iff \min -f(\mathbf{x})$.
- All constraint right-hand sides are normalized so that $\mathbf{b} \ge 0$. Any row with $\text{RHS} < 0$ is multiplied by $-1$ with its inequality reversed before tableau construction.

### 2. Tableau Layout & RHS Placement
In this implementation, the tableau columns are arranged as:
$$\begin{array}{|c|c|cccc|}
\hline
\textbf{Basis} & \textbf{RHS} & x_1 & x_2 & \dots & x_n \\
\hline
-z & \text{RHS}_0 & \bar{c}_1 & \bar{c}_2 & \dots & \bar{c}_n \\
\hline
x_{B(1)} & \text{RHS}_1 & a_{1,1} & a_{1,2} & \dots & a_{1,n} \\
\vdots & \vdots & \vdots & \vdots & \ddots & \vdots \\
x_{B(m)} & \text{RHS}_m & a_{m,1} & a_{m,2} & \dots & a_{m,n} \\
\hline
\end{array}$$

- **Column 1 (`Basis`)**: Name of the basic variable for that row (`-z` or `-w` for Row 0; $x_i$, $s_i$, or $y_i$ for constraint rows).
- **Column 2 (`RHS`)**: Constant term of the row equation.
  - For Row 0: represents the raw identity $-z = \text{RHS}_0 \implies z = -\text{RHS}_0$.
  - For Row $i \ge 1$: indicates the value of basic variable $x_{B(i)} = \text{RHS}_i$. Non-basic variables are $0$.
- **No $z$ column**: The identity column for $z$ is eliminated to save space and match textbook practice.

### 3. Bland's Smallest-Subscript Anti-Cycling Rule
Dantzig's standard "most negative reduced cost" rule can cycle indefinitely on degenerate linear programs. Bland's rule guarantees finite termination without cycling:

1. **Entering Variable Selection**:
   - Compute candidate non-basic variables with strictly negative Row 0 reduced costs: $\mathcal{C} = \{ j \mid \bar{c}_j < 0 \}$.
   - If $\mathcal{C} = \emptyset$, current tableau is **optimal**.
   - Otherwise, select the candidate variable with the **smallest subscript index**:
     $$j^* = \arg\min_{j \in \mathcal{C}} \; \text{index}(x_j)$$

2. **Leaving Variable Selection (Minimum $\theta$ Ratio Test)**:
   - For all constraint rows $i \in \{1, \dots, m\}$ where the pivot column entry is strictly positive ($a_{i, j^*} > 0$), calculate:
     $$\theta_i = \frac{\text{RHS}_i}{a_{i, j^*}}$$
   - If no $a_{i, j^*} > 0$ exists, the linear program is **unbounded**.
   - Compute the minimum ratio: $\theta^* = \min \{ \theta_i \mid a_{i, j^*} > 0 \}$.
   - If there is a tie ($\theta_i = \theta^*$ for multiple rows), choose the basic variable with the **smallest subscript index**:
     $$i^* = \arg\min_{\{i \mid \theta_i = \theta^*\}} \; \text{index}(x_{B(i)})$$

3. **Subscript Index Ordering**:
   Variables are assigned a strict, 1-based sequential index:
   $$\text{Decision Variables } (x_1, \dots, x_n) < \text{Slack Variables } (x_{n+1}, \dots) < \text{Surplus } (e_1, \dots) < \text{Artificial } (y_1, \dots)$$

### 4. Two-Phase Simplex
When constraints contain $\ge$ or $=$ relations (or when `--force-two-phase` is requested), basic slack variables cannot provide an initial feasible basis:
- **Phase 1**:
  - Introduces non-negative artificial variables $y_1, \dots, y_k$.
  - Row 0 minimizes the sum of artificials: $-w = -\sum y_i$.
  - Canonicalizes Row 0 against the initial artificial basis so reduced costs of $y_i$ become $0$.
  - Solves using Bland's rule. If $\min W > 0$, the problem is **infeasible**. If $\min W = 0$, a feasible basis has been found.
  - If any $y_i$ remains in the basis at value $0$, it is pivoted out with a non-artificial non-basic variable.
- **Phase 2**:
  - Drops all artificial columns $y_i$.
  - Restores the original objective function Row 0 ($-z$).
  - Canonicalizes Row 0 with respect to the Phase 1 final basis.
  - Continues Bland's simplex iterations until optimality.

---

## 📊 CSV Tableau Format & Specifications

CSV files stored in [`csv/`](csv/) represent initial tableaus.

### Example: `csv/tableau.csv` (Beale's Problem)
```csv
Basis,RHS,x1,x2,x3,x4,x5,x6,x7
-z,0,-10,57,9,24,0,0,0
x5,0,1/2,-11/2,-5/2,9,1,0,0
x6,0,1/2,-3/2,-1/2,1,0,1,0
x7,1,1,0,0,0,0,0,1
```

### Format Specifications:
1. **Headers (Row 1)**: `Basis, RHS, var1, var2, ...`
2. **Row 0 (Row 2)**: Objective row labeled `-z` (or `z`, `-w`, `w`).
3. **Rows $1 \dots m$**: Constraint rows labeled with basic variables.
4. **Delimiters**: Auto-detects comma (`,`), semicolon (`;`), or tab (`\t`).
5. **Number Parsing**: Exact integer, standard fraction (`1/2`, `-11/2`), decimal with point (`0.5`), or decimal with European comma (`0,5`).
6. **Automatic Normalization**:
   - Legacy slack columns labeled `s1, s2, ...` are automatically converted to sequential $x$ variables (`x5, x6, ...`).
   - If basic variables are omitted in column 0, they are inferred automatically via identity unit vectors.
   - Older CSV tableaus with RHS on the far right or containing an explicit $z$ column are auto-converted to the standard format.

---

## 🚀 Execution & Usage

### 1. Launch Interactive Web GUI
Run from repository root:
```bash
python3 source/gui.py
```
Or from within `source/`:
```bash
cd source && python3 gui.py
```
- Starts local HTTP server (automatically picks a free port starting at `5050`).
- Automatically launches your default web browser at `http://127.0.0.1:<port>`.

### 2. Command Line Interface (CLI Mode)
Solve any CSV tableau directly in the terminal without opening a browser:
```bash
# Relative to csv/ directory
python3 source/gui.py csv/tableau.csv --cli

# Or simply specify filename (auto-resolved from csv/)
python3 source/gui.py tableau.csv --cli

# Solve standard textbook 2D LP
python3 source/gui.py example_standard.csv --cli
```

---

## 🌐 Embedded REST API Reference

When `source/gui.py` runs, it exposes a lightweight JSON HTTP API. This allows AI subagents, headless scripts, or automated tests to interact with the solver programmatically without a browser.

| Endpoint | Method | Parameters / Body | Response Description |
|---|---|---|---|
| `/` | `GET` | None | Returns the full single-page web GUI application (`text/html`). |
| `/api/csv-files` | `GET` | None | Returns JSON array of available CSV filenames found in `csv/` and root: `["tableau.csv", "example_standard.csv"]`. |
| `/api/load-csv` | `GET` | `?file=<filename>` | Returns parsed `{ "headers": [...], "rows": [[...]] }`. Automatically resolves path in `csv/`. |
| `/api/solve-file` | `GET` | `?file=<filename>&force_two_phase=<bool>` | Solves the specified CSV tableau. Returns `{ "status": "OPTIMAL", "steps": [...], "final_solution": {...} }`. |
| `/api/demo` | `GET` | `?name=<beale\|twophase\|standard2d>` | Solves one of the built-in reference benchmarks and returns full trajectory. |
| `/api/solve-tableau` | `POST` | JSON: `{ "headers": [...], "rows": [[...]], "force_two_phase": bool }` | Solves an in-memory edited tableau. |
| `/api/solve-lp` | `POST` | JSON: `{ "obj_type": "min", "obj_coeffs": [...], "constraints": [...], "force_two_phase": bool }` | Constructs standard form, adds artificials if needed, and solves. |
| `/api/save-csv` | `POST` | JSON: `{ "filename": "name.csv", "content": "..." }` | Saves or overwrites a CSV tableau into the `csv/` directory. |

### Step Object Schema in Response (`steps` array)
Each item in `steps` represents an elementary sub-step:
```json
{
  "iteration": 0,
  "substep_idx": 1,
  "substep_total": 4,
  "substep_type": "pivot_selection",
  "title": "Iteration 0 • Step 1: Select Pivot",
  "substep_desc": "Entering: 'x1' ... Leaving: 'x6' ...",
  "formula": "Pivot: a(R2, x1) = 1/2, min θ = 0",
  "pivot": { "row": 2, "col": 0, "val": "1/2" },
  "active_rows": [2],
  "active_cols": [0],
  "var_names": ["x1", "x2", "x3", "x4", "x5", "x6", "x7"],
  "basic_vars": ["x5", "x6", "x7"],
  "obj_name": "-z",
  "rows": [
    { "basic": "-z", "rhs": "0", "cells": ["-10", "57", "9", "24", "0", "0", "0"] },
    { "basic": "x5", "rhs": "0", "cells": ["1/2", "-11/2", "-5/2", "9", "1", "0", "0"] }
  ],
  "status": "CONTINUE",
  "candidates": [
    { "name": "x1", "index": 1, "col": 0, "coeff": "-10" }
  ],
  "ratios": [
    { "row": 1, "basic": "x5", "index": 5, "theta": "0", "valid": true, "rhs": "0", "coeff": "1/2" },
    { "row": 2, "basic": "x6", "index": 6, "theta": "0", "valid": true, "rhs": "0", "coeff": "1/2" }
  ],
  "min_ratio": "0",
  "tied": true,
  "bland_tie_note": "Tie detected between 2 rows with minimum θ = 0. Bland's Rule selects basic variable with smallest index: 'x5' (Index 5)."
}
```

---

## 🐍 Programmatic Python API

Any script or test can import and use the engine directly without network overhead:

```python
import sys
sys.path.insert(0, "source")
import gui

# 1. Load and solve a CSV tableau directly:
tableau, obj_type = gui.load_tableau_from_csv("tableau.csv") # Auto-resolves to csv/tableau.csv
result = gui.record_tableau_steps_twophase(tableau, obj_type="min")

print("Status:", result["status"])
print("Optimal Objective:", result["final_solution"]["obj_val"])
print("Solution Vector:", result["final_solution"]["solution_vector"])

# 2. Formulate and solve an arbitrary LP from equations:
lp = gui.LinearProgram(
    obj_type="min",
    obj_coeffs=[gui.parse_fraction("-3"), gui.parse_fraction("-5")],
    constraints=[
        ([gui.parse_fraction("1"), gui.parse_fraction("0")], "<=", gui.parse_fraction("4")),
        ([gui.parse_fraction("0"), gui.parse_fraction("2")], "<=", gui.parse_fraction("12")),
        ([gui.parse_fraction("3"), gui.parse_fraction("2")], "<=", gui.parse_fraction("18")),
    ],
    var_names=["x1", "x2"],
    force_two_phase=False
)
lp_result = gui.record_lp_steps(lp)
print("LP Optimal Status:", lp_result["status"])
print("LP Min z*:", lp_result["final_solution"]["obj_val"]) # -36
```

---

## 🏆 Reference Benchmark Verification Cases

### 1. Beale's Cycling Benchmark (1955)
- **Problem**:
  $$\min \; z = -10x_1 + 57x_2 + 9x_3 + 24x_4$$
  Subject to:
  $$\begin{aligned}
  \frac{1}{2}x_1 - \frac{11}{2}x_2 - \frac{5}{2}x_3 + 9x_4 &\le 0 \\
  \frac{1}{2}x_1 - \frac{3}{2}x_2 - \frac{1}{2}x_3 + x_4 &\le 0 \\
  x_1 &\le 1 \\
  x_1, x_2, x_3, x_4 &\ge 0
  \end{aligned}$$
- **Dantzig Rule Behavior**: Cycles indefinitely among 6 bases at objective value $0$.
- **Bland's Rule Behavior**: Resolves ties deterministically, avoiding cycling and reaching optimality in 7 iterations.
- **Optimal Solution**: $x_1 = 1$, $x_2 = 0$, $x_3 = 1$, $x_4 = 0$, $x_5 = 2$, $x_6 = 0$, $x_7 = 0$, with $\min z^* = -1$.

### 2. Standard 2D Textbook LP
- **Problem**:
  $$\min \; z = -3x_1 - 5x_2 \quad \text{s.t.} \quad x_1 \le 4, \quad 2x_2 \le 12, \quad 3x_1 + 2x_2 \le 18, \quad x_1, x_2 \ge 0$$
- **Optimal Solution**: $x_1 = 2$, $x_2 = 6$, slack $x_3 = 2$, $x_4 = 0$, $x_5 = 0$, with $\min z^* = -36$.

### 3. Two-Phase LP with Artificial Variables
- **Problem**:
  $$\min \; z = 4x_1 + x_2 \quad \text{s.t.} \quad 3x_1 + x_2 = 3, \quad 4x_1 + 3x_2 \ge 6, \quad x_1 + 2x_2 \le 4, \quad x_1, x_2 \ge 0$$
- **Phase 1**: Minimizes artificial variables $y_1, y_2$ to $0$.
- **Phase 2**: Optimizes original objective from feasible basis.
- **Optimal Solution**: $x_1 = 2/5 = 0.4$, $x_2 = 9/5 = 1.8$, with $\min z^* = 17/5 = 3.4$.

---

## 📜 License

MIT License. See [LICENSE](LICENSE) for full text.
