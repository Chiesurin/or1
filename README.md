# OR1: Operations Research Linear & Integer Programming Toolkit

**OR1** (*Operations Research 1*) is an open-source, exact-arithmetic toolkit for solving, verifying, and visualizing **Linear Programming (LP)** and **Integer Linear Programming (ILP)** problems.

It combines a self-contained **Web GUI Visualizer & REST API**, a **Terminal CLI**, and two deterministic **AI Agent Skills** built around exact rational arithmetic (`fractions.Fraction`), **Bland's anti-cycling rule**, **Two-Phase Simplex**, **LP Duality & Dual Simplex**, **Total Unimodularity (TUM)** analysis, **Chvátal-Gomory / GMI cutting planes**, **Branch and Bound**, and **Branch and Cut**.

---

## Core Capabilities

- **Exact Rational Arithmetic**: Uses Python's built-in `fractions.Fraction` across every pivot, ratio test, dual inversion, and cutting plane. Eliminates floating-point drift (e.g., `103/26`, `-11/2`, `17/5`).
- **Anti-Cycling Simplex & Two-Phase Engine**: Implements Bland's smallest-subscript rule for both entering and leaving variable ties (verified on Beale's 1955 cycling benchmark), plus automatic or forced Two-Phase Simplex with artificial variables `y1, y2, ...`.
- **Duality Theory & Dual Simplex**:
  - Automatic formulation of the mathematical dual problem ($\min \leftrightarrow \max$, transposed $A^T$, exact sign/free bounds).
  - Exact shadow price computation via basis inversion ($B^T y^* = c_B$) for both `--lp` and `--csv`/`--raw` tableaus.
  - Explicit verification of **Strong Duality** ($z^* = w^*$) and **Complementary Slackness**.
  - Full step-by-step **Dual Simplex** execution for dual-feasible tableaus ($c_j \ge 0$) with negative RHS entries.
- **Integer & Mixed-Integer Linear Programming (ILP / MILP)**:
  - Pre-solve **Total Unimodularity (TUM)** inspection: checks entry necessity ($a_{ij} \in \{-1, 0, 1\}$), Poincaré / Heller-Tompkins network criteria, Consecutive Ones (interval matrix) property, and Bareiss exact submatrix determinants.
  - **Hoffman-Kruskal Theorem** evaluation to identify Ideal Polyhedra ($P = P_I$) before introducing cuts or branches.
  - **Chvátal-Gomory & Gomory Mixed-Integer (GMI) Cutting Planes** derived both in tableau slack form and re-projected into original decision variables, followed by Dual Simplex re-optimization.
  - **Branch and Bound** supporting `best-bound`, `dfs`, and `bfs` node selection, `most-fractional` and `first-fractional` branching, and exact pruning (`INFEASIBLE`, `PRUNED_BOUND`, `INTEGER_INCUMBENT`, `INTEGER_FATHOMED`).
  - **Branch and Cut** combining Two-Phase LP node relaxation, per-node CG/GMI cut generation with Dual Simplex re-optimization, and subtree branching.
  - **Interactive 2D HTML/SVG Plotter** showing the continuous LP polyhedron, integer lattice points / MILP slices, branching bounds, cutting hyperplanes, node tree navigation, continuous optimum $x^*_{\text{LP}}$, and integer optimum $x^*_{\text{ILP}}$.
- **Optional SciPy Cross-Verification**: Runs with zero mandatory third-party packages; when `scipy` and `numpy` are installed, solutions are automatically cross-checked against HiGHS (`scipy.optimize.linprog` and `scipy.optimize.milp`).
- **Publication-Ready Exports**: Generates standalone Markdown walkthroughs (`--output-md`) and interactive single-file HTML visualizations (`--visualize-2d`).

---

## Requirements

- **Runtime**: Python 3.8 or newer.
- **Mandatory Dependencies**: None (Python Standard Library only).
- **Optional Dependencies** (for HiGHS cross-verification):

```bash
pip3 install --user scipy numpy
```

---

## Quick Start

```bash
# 1. Launch the Interactive Web GUI (opens http://127.0.0.1:5050)
python3 source/gui.py

# 2. Solve a CSV tableau directly in the terminal
python3 source/gui.py csv/tableau.csv --cli

# 3. Solve an LP with Dual formulation, Shadow Prices, and SciPy verification
python3 skills/tableau-solver/scripts/solve_tableau.py \
  --lp "min 4*x1 + x2 s.t. 3*x1 + x2 = 3; 4*x1 + 3*x2 >= 6; x1 + 2*x2 <= 4; x1, x2 >= 0" \
  --show-dual --verify-scipy

# 4. Solve a 2D ILP with TUM check, Gomory cuts, and generate an interactive 2D HTML plot
python3 skills/ilp-solver/scripts/solve_ilp.py \
  --lp "max 2*x1 + 3*x2 s.t. 2*x1 + 5*x2 <= 16; 6*x1 + 2*x2 <= 27; x1, x2 >= 0" \
  --visualize-2d cuts_2d.html --verify-scipy
```

---

## Repository Structure

```text
.
├── .agents/
│   └── skills/                         # Auto-discovery skill symlinks for AI coding agents
│       ├── ilp-solver -> ../../skills/ilp-solver
│       └── tableau-solver -> ../../skills/tableau-solver
├── .gitattributes
├── .gitignore
├── INSTRUCTIONS.md                     # Detailed user guide and learning walkthrough
├── LICENSE                             # MIT License
├── README.md                           # Project overview, specs, CLI & API reference
├── csv/                                # Sample initial Simplex tableau datasets
│   ├── example_standard.csv            # Standard 2D minimization tableau (opt z* = -33)
│   └── tableau.csv                     # Beale's 1955 cycling benchmark (opt z* = -1)
├── skills/
│   ├── ilp-solver/                     # ILP/MILP TUM, Cutting Planes, B&B & B&C skill
│   │   ├── SKILL.md                    # Agent specification and workflow
│   │   ├── examples/
│   │   │   ├── bnb_2d.lp               # Classic 2D Branch & Bound / Branch & Cut ILP
│   │   │   ├── gomory_2d.lp            # Classical 2D non-TUM ILP example
│   │   │   ├── gomory_2d_viz.html      # Interactive 2D cut visualization HTML output
│   │   │   ├── milp_2d.lp              # 2D Mixed-Integer LP (integer + continuous)
│   │   │   └── tum_network.lp          # Ideal polyhedron (TUM) benchmark example
│   │   ├── references/
│   │   │   ├── branch_and_bound_cut.md # B&B, Branch & Cut, and GMI cuts reference
│   │   │   ├── gomory_cuts.md          # Chvátal-Gomory & Dual Simplex derivation notes
│   │   │   └── total_unimodularity.md  # TUM criteria & Hoffman-Kruskal reference
│   │   └── scripts/
│   │       └── solve_ilp.py            # Deterministic ILP/MILP, TUM, B&B, B&C & 2D engine
│   └── tableau-solver/                 # LP Simplex, Two-Phase & Duality agent skill
│       ├── SKILL.md                    # Agent specification and workflow
│       ├── examples/
│       │   ├── input_inline_prompt.txt # Sample inline prompt input
│       │   ├── sample_solution.md      # Sample exported Markdown walkthrough
│       │   ├── sample_tableau.csv      # 2D textbook tableau (opt z* = -36)
│       │   └── test_run.sh             # End-to-end verification shell script
│       ├── references/
│       │   ├── input_formats.md        # Accepted CSV and LP syntax specifications
│       │   └── simplex_cheatsheet.md   # Mathematical conventions & pivoting rules
│       └── scripts/
│           └── solve_tableau.py        # Deterministic Simplex, Two-Phase & Duality engine
└── source/
    └── gui.py                          # Standalone Web GUI server, REST API & CSV CLI
```

---

## Execution Modes

### 1. Interactive Web GUI & CSV CLI (`source/gui.py`)

Launch the local web application:

```bash
python3 source/gui.py
```

- Automatically selects the first available local port starting at `5050` (`http://127.0.0.1:5050`) and binds strictly to loopback `127.0.0.1`.
- Provides color-coded tableau rendering (basic vs. non-basic columns, artificial variables `y_i`, Row 0 reduced costs, RHS), sub-step playback (pivot selection, row normalization, Gauss-Jordan elimination), interactive CSV editor backed by `csv/`, manual LP entry form, built-in benchmarks (`beale`, `twophase`, `standard2d`), and automated SciPy HiGHS correctness checking.

To solve any CSV tableau directly in the terminal without starting the HTTP server:

```bash
python3 source/gui.py csv/tableau.csv --cli
python3 source/gui.py example_standard.csv --cli
```

---

### 2. LP & Duality Solver CLI (`skills/tableau-solver/scripts/solve_tableau.py`)

Designed for both direct terminal use and zero-mental-math execution by LLM agents:

```bash
# Solve from CSV (auto-searches workspace root, csv/, and examples/)
python3 skills/tableau-solver/scripts/solve_tableau.py --csv tableau.csv --verify-scipy

# Solve from raw CSV/text inline
python3 skills/tableau-solver/scripts/solve_tableau.py \
  --raw "Basis,RHS,x1,x2\n-z,0,-3,-5\nx3,4,1,0\nx4,12,0,2\nx5,18,3,2"

# Solve algebraic LP, display Dual LP, and export full Markdown report
python3 skills/tableau-solver/scripts/solve_tableau.py \
  --lp "min 4*x1 + x2 s.t. 3*x1 + x2 = 3; 4*x1 + 3*x2 >= 6; x1 + 2*x2 <= 4; x1, x2 >= 0" \
  --show-dual --output-md solution.md --verify-scipy

# Run via Dual Simplex (requires initial dual-feasible tableau with c_j >= 0)
python3 skills/tableau-solver/scripts/solve_tableau.py \
  --lp "min 2*x1 + 3*x2 s.t. x1 + x2 >= 2; 2*x1 + x2 >= 3; x1, x2 >= 0" \
  --dual-simplex

# Pipe from stdin and emit structured JSON
cat csv/tableau.csv | python3 skills/tableau-solver/scripts/solve_tableau.py --stdin --format json
```

| Flag | Description |
|---|---|
| `--csv <path>` | Load initial tableau from a CSV file. |
| `--raw "<text>"` | Parse multiline or `\n`-separated tableau text directly from CLI. |
| `--stdin` | Read tableau text from standard input. |
| `--lp "<text>"` | Parse algebraic LP text (`min`/`max`, constraints separated by newlines or `;`). |
| `--two-phase` | Force Two-Phase Simplex even when a slack basis already exists. |
| `--dual-simplex` | Run Dual Simplex (for tableaus with $c_j \ge 0$ and negative RHS). |
| `--show-dual` | Print the mathematical dual LP formulation (when using `--lp`). |
| `--output-md <path>` | Export a complete Markdown walkthrough report. |
| `--format {text,markdown,json}` | Select stdout format (default: `text`). |
| `--verify-scipy` | Cross-verify optimal value and status against `scipy.optimize.linprog`. |

---

### 3. Integer & Mixed-Integer LP Solver (`skills/ilp-solver/scripts/solve_ilp.py`)

Analyzes Total Unimodularity and solves pure ILPs and MILPs using **Chvátal-Gomory / GMI Cutting Planes**, **Branch and Bound**, or **Branch and Cut**:

```bash
# 1. Pure Cutting Planes on a built-in example ('gomory2d', 'tum', 'bnb', 'milp')
python3 skills/ilp-solver/scripts/solve_ilp.py --example gomory2d --verify-scipy

# 2. Branch and Bound with Best-Bound search and 2D tree visualizer
python3 skills/ilp-solver/scripts/solve_ilp.py \
  --example bnb --method branch-and-bound \
  --visualize-2d bnb_2d_viz.html --verify-scipy

# 3. Branch and Cut (CG/GMI cuts + Dual Simplex + Branching)
python3 skills/ilp-solver/scripts/solve_ilp.py \
  --example bnb --method branch-and-cut --verify-scipy

# 4. Mixed-Integer LP solved with Branch and Cut (GMI cuts)
python3 skills/ilp-solver/scripts/solve_ilp.py \
  --example milp --method branch-and-cut --output-md milp_report.md --verify-scipy

# 5. Solve from an .lp file and export interactive 2D cut visualizer
python3 skills/ilp-solver/scripts/solve_ilp.py \
  --lp skills/ilp-solver/examples/gomory_2d.lp \
  --visualize-2d skills/ilp-solver/examples/gomory_2d_viz.html \
  --output-md ilp_report.md --verify-scipy
```

| Flag | Description |
|---|---|
| `--lp <text_or_path>` | Algebraic formulation string or path to `.lp` file (supports `integer`, `binary`, `continuous`). |
| `--example {gomory2d,tum,bnb,milp}` | Run a built-in benchmark problem. |
| `--method {cutting-plane,cg,branch-and-bound,bnb,branch-and-cut,bnc}` | Solving algorithm (default: `cutting-plane`). |
| `--strategy {best-bound,dfs,bfs}` | Node selection rule for Branch & Bound / Branch & Cut (default: `best-bound`). |
| `--branch-rule {most-fractional,first-fractional}` | Variable selection rule for branching (default: `most-fractional`). |
| `--cut-type {auto,cg,gmi}` | Cut generator mode (default: `auto`: CG when integer-safe, GMI otherwise). |
| `--max-cuts <N>` | Maximum cutting plane iterations in `cutting-plane` mode (default: `15`). |
| `--cuts-per-node <N>` | Max cut rounds per tree node in `branch-and-cut` (default: `2`). |
| `--cut-depth <N>` | Max tree depth for cut generation in `branch-and-cut` (default: `5`). |
| `--max-nodes <N>` | Maximum explored nodes in tree search (default: `200`). |
| `--visualize-2d <path.html>` | Export standalone 2D HTML/SVG visualizer across all three methods (when $n = 2$). |
| `--output-md <path.md>` | Export TUM analysis, tree/cut trace, and optimal solution to Markdown. |
| `--format {text,json}` | Select terminal stdout format (default: `text`). |
| `--verify-scipy` | Cross-verify optimal objective and status against `scipy.optimize.milp` (HiGHS). |

---

## Mathematical Conventions

### 1. Universal Minimization & Row 0 Identity

In the Simplex tableau engines, problems are solved in standard minimization form:

$$
\min_{\mathbf{x}} \; \mathbf{c}^T \mathbf{x} \quad \text{subject to} \quad A\mathbf{x} = \mathbf{b}, \quad \mathbf{x} \ge 0
$$

- Maximization problems $\max \mathbf{c}^T \mathbf{x}$ are transformed internally via $\min (-\mathbf{c}^T \mathbf{x})$.
- Row 0 is labeled `-z` (or `-w` in Phase 1) and encodes:

$$
-z + \sum_{j \in N} \bar{c}_j x_j = \text{RHS}_0 \implies z = -\text{RHS}_0
$$

- Constraint rows with $\text{RHS} < 0$ are multiplied by $-1$ (reversing inequality sense) prior to primal/two-phase setup, except when running `--dual-simplex` directly on a dual-feasible tableau.

### 2. Tableau Layout

Columns follow left-hand RHS placement without a separate $z$ column:

```text
+-------+-------+------+------+-----+------+
| Basis |  RHS  |  x1  |  x2  | ... |  xn  |
+-------+-------+------+------+-----+------+
|  -z   | RHS_0 | c̄_1  | c̄_2  | ... | c̄_n  |
| xB(1) | RHS_1 | a_11 | a_12 | ... | a_1n |
|  ...  |  ...  | ...  | ...  | ... | ...  |
| xB(m) | RHS_m | a_m1 | a_m2 | ... | a_mn |
+-------+-------+------+------+-----+------+
```

### 3. Bland's Anti-Cycling Rule

All variables are assigned a strict 1-based subscript index:

$$
\text{Decision } (x_1, \dots, x_n) < \text{Slack } (x_{n+1}, \dots) < \text{Surplus } (e_1, \dots) < \text{Artificial } (y_1, \dots)
$$

1. **Entering Variable**: Among non-basic variables with $\bar{c}_j < 0$, choose the one with the smallest index. If none exist, the tableau is optimal.
2. **Leaving Variable**: Among rows with $a_{ij^*} > 0$, compute $\theta_i = \text{RHS}_i / a_{ij^*}$. If no entry is positive, the LP is unbounded. On a minimum-ratio tie, choose the basic variable with the smallest index.

### 4. Duality, Shadow Prices & Dual Simplex

- Given optimal basis $B$, the exact dual solution is computed as $y^* = B^{-T} c_B$, satisfying Strong Duality ($z^* = w^*$) and Complementary Slackness ($y_i^* s_i^* = 0$, $x_j^* \bar{c}_j = 0$).
- Under **Dual Simplex** ($\bar{c}_j \ge 0$ for all $j$), the leaving row $l$ is chosen by most negative RHS (ties broken by Bland's index), and the entering column minimizes $\theta_j = \bar{c}_j / |a_{lj}|$ over $a_{lj} < 0$.

### 5. Total Unimodularity & Chvátal-Gomory Cuts

- **Total Unimodularity (TUM)**: Before solving the ILP cutting loop, `solve_ilp.py` tests whether $A$ is TUM and $b \in \mathbb{Z}^m$. By the **Hoffman-Kruskal Theorem**, when $A$ is TUM and $b$ is integral, $P = \{x \ge 0 : Ax \le b\}$ is an **Ideal Polyhedron** ($P = P_I = \text{conv}(P \cap \mathbb{Z}^n)$), so the continuous LP optimum is automatically integral with zero cuts required.
- **Chvátal-Gomory Cut Derivation**: From any fractional basic row $x_B + \sum_{j \in N} \bar{a}_{ij} x_j = \bar{b}_i$ with $f_0 = \bar{b}_i - \lfloor \bar{b}_i \rfloor \in (0, 1)$ and $f_j = \bar{a}_{ij} - \lfloor \bar{a}_{ij} \rfloor \in [0, 1)$, the valid cut is:

$$
\sum_{j \in N} f_j x_j \ge f_0 \iff -s_{\text{cut}} + \sum_{j \in N} (-f_j) x_j = -f_0
$$

- Substituting slack definitions yields the equivalent cut in original decision variables, while appending the `-f_0` row preserves dual feasibility so **Dual Simplex** immediately re-optimizes.

---

## CSV Tableau Format

Files in [`csv/`](csv/) store initial tableaus with auto-detected formatting:

- **Delimiters**: Comma (`,`), semicolon (`;`), or tab (`\t`).
- **Numeric Values**: Integers, rationals (`1/2`, `-11/2`), dot decimals (`0.5`), and European comma decimals (`0,5`).
- **Auto-Normalization**:
  - Slack labels `s1, s2, ...` are normalized into sequential `x` variables (`x5, x6, ...`).
  - Missing basic variable labels are inferred from unit identity columns.
  - Legacy tableaus with right-hand RHS or explicit `z` columns are automatically normalized on load.

Example (`csv/tableau.csv` — Beale's cycling benchmark):

```csv
Basis,RHS,x1,x2,x3,x4,x5,x6,x7
-z,0,-10,57,9,24,0,0,0
x5,0,1/2,-11/2,-5/2,9,1,0,0
x6,0,1/2,-3/2,-1/2,1,0,1,0
x7,1,1,0,0,0,0,0,1
```

---

## Verified Benchmark Datasets

| Dataset / Problem | Solver | Optimal Objective | Optimal Solution | Notes |
|---|---|---|---|---|
| [csv/tableau.csv](csv/tableau.csv) | Simplex (Bland) | $\min z^* = -1$ | $x_1=1, x_2=0, x_3=1, x_4=0, x_5=2$ | Beale (1955) cycling test; 7 pivots; dual $y^*=(0, -18, -1)$ |
| [csv/example_standard.csv](csv/example_standard.csv) | Simplex (Bland) | $\min z^* = -33$ | $x_1=3, x_2=12, x_3=0, x_4=0, x_5=3$ | Standard 2D LP; dual $y^*=(-5/4, -1/4, 0)$ |
| [skills/tableau-solver/examples/sample_tableau.csv](skills/tableau-solver/examples/sample_tableau.csv) | Simplex (Bland) | $\min z^* = -36$ | $x_1=2, x_2=6, x_3=2, x_4=0, x_5=0$ | Classic textbook LP ($\min -3x_1 - 5x_2$) |
| Two-Phase Benchmark (`--demo twophase`) | Two-Phase Simplex | $\min z^* = 17/5$ | $x_1=2/5, x_2=9/5, e_1=1$ | Mixed $=, \ge, \le$ constraints; dual $y^*=(7/5, 0, -1/5)$ |
| [skills/ilp-solver/examples/gomory_2d.lp](skills/ilp-solver/examples/gomory_2d.lp) | Chvátal-Gomory ILP | $\max z^*_{\text{ILP}} = 12$ | $x_1=3, x_2=2$ | Continuous LP $z_{\text{LP}} = 166/13$ at $(103/26, 21/13)$; 3 CG cuts (or 3 nodes, 2 cuts under B&C) |
| [skills/ilp-solver/examples/bnb_2d.lp](skills/ilp-solver/examples/bnb_2d.lp) | Branch & Bound / Branch & Cut | $\max z^*_{\text{ILP}} = 40$ | $x_1=0, x_2=5$ | Root $z_{\text{LP}} = 165/4$; 9 nodes under B&B vs. 1 node + 1 CG cut under B&C |
| [skills/ilp-solver/examples/milp_2d.lp](skills/ilp-solver/examples/milp_2d.lp) | MILP (B&B / GMI B&C) | $\max z^*_{\text{MILP}} = 55/3$ | $x_1=3, x_2=7/3$ | Mixed $x_1 \in \mathbb{Z}, x_2 \in \mathbb{R}_{\ge 0}$; 3 nodes under B&B vs. 1 node + 1 GMI cut under B&C |
| [skills/ilp-solver/examples/tum_network.lp](skills/ilp-solver/examples/tum_network.lp) | TUM + LP Relaxation | $\max z^*_{\text{ILP}} = 5$ | $x_1=3, x_2=2$ | Verified TUM & Ideal Polyhedron ($P = P_I$); 0 cuts needed |

---

## Embedded REST API (`source/gui.py`)

When the Web GUI runs on `http://127.0.0.1:<port>`, it exposes a JSON API restricted to loopback connections. CSV file operations are confined strictly to plain `.csv` basenames inside `csv/` or the repository root.

| Endpoint | Method | Parameters / JSON Body | Description |
|---|---|---|---|
| `/` | `GET` | None | Serves the single-page Web GUI HTML/CSS/JS application. |
| `/api/csv-files` | `GET` | None | Lists available `.csv` basenames in `csv/` and root. |
| `/api/load-csv` | `GET` | `?file=<name.csv>` | Returns parsed `{ "headers": [...], "rows": [[...]] }`. |
| `/api/solve-file` | `GET` | `?file=<name.csv>&force_two_phase=<bool>` | Solves a CSV tableau; returns `{ "status", "steps", "final_solution", "scipy", "correctness" }`. |
| `/api/demo` | `GET` | `?name=<beale\|twophase\|standard2d>` | Solves a built-in benchmark and returns its step trajectory. |
| `/api/solve-tableau` | `POST` | `{ "headers": [...], "rows": [[...]], "force_two_phase": bool }` | Solves an in-memory edited tableau. |
| `/api/solve-lp` | `POST` | `{ "obj_type": "min", "obj_coeffs": [...], "constraints": [...], "force_two_phase": bool }` | Converts an LP to standard form and returns full trajectory and SciPy check. |
| `/api/save-csv` | `POST` | `{ "filename": "<name.csv>", "content": "..." }` | Saves a `.csv` tableau into `csv/`. |

---

## Programmatic Python Usage

```python
import sys
sys.path.insert(0, "source")
import gui

# 1. Load and solve a CSV tableau directly
tableau, obj_type = gui.load_tableau_from_csv("tableau.csv")
result = gui.record_tableau_steps_twophase(tableau, obj_type="min")
print("Status:", result["status"])
print("Optimal Objective:", result["final_solution"]["obj_val"])
print("SciPy Verification:", result["correctness"]["status"])

# 2. Construct and solve a Linear Program in Python
lp = gui.LinearProgram(
    obj_type="min",
    obj_coeffs=[gui.parse_fraction("-3"), gui.parse_fraction("-2")],
    constraints=[
        ([gui.parse_fraction("2"), gui.parse_fraction("1")], "<=", gui.parse_fraction("18")),
        ([gui.parse_fraction("2"), gui.parse_fraction("3")], "<=", gui.parse_fraction("42")),
        ([gui.parse_fraction("3"), gui.parse_fraction("1")], "<=", gui.parse_fraction("24")),
    ],
    var_names=["x1", "x2"],
)
lp_res = gui.record_lp_steps(lp)
print("Min z*:", lp_res["final_solution"]["raw_obj_val"])  # -33
```

---

## AI Coding Agent Integration

Relative symlinks in `.agents/skills/` expose [`skills/tableau-solver/SKILL.md`](skills/tableau-solver/SKILL.md) and [`skills/ilp-solver/SKILL.md`](skills/ilp-solver/SKILL.md) to compatible AI coding assistants (such as Antigravity, Claude Code, Codex, and Cursor). Agents automatically invoke the Python scripts instead of performing approximate mental arithmetic.

You can also run the bundled skill test suite with:

```bash
bash skills/tableau-solver/examples/test_run.sh
```

---

## License

Licensed under the [MIT License](LICENSE).
