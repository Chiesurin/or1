---
name: ilp-solver
description: Solves Integer Linear Programs (ILP) and Mixed-Integer Linear Programs (MILP) using the Chvátal-Gomory / GMI Cutting Plane Algorithm, Branch and Bound, and Branch and Cut, performs Total Unimodularity (TUM) and Ideal Polyhedron analysis, and generates interactive 2D graphical visualizations. Use when the user asks to solve an ILP or MILP, use Branch and Bound or Branch and Cut, analyze Total Unimodularity, derive Gomory or Chvátal cuts, inspect the integer hull, or visualize trees and cutting planes on 2D polyhedra.
---

# Integer & Mixed-Integer Linear Programming (ILP / MILP) Solver Skill

This skill provides deterministic exact-arithmetic solving (`fractions.Fraction`) and visual explanations for Integer and Mixed-Integer Linear Programming problems. It implements **Total Unimodularity (TUM)** pre-solve analysis, **Chvátal-Gomory & Gomory Mixed-Integer (GMI) Cutting Planes**, **Branch and Bound**, **Branch and Cut**, and **2D SVG/HTML interactive visualization**.

---

## Capabilities & Triggers

Activate this skill when the user:
1. Asks to solve an **Integer Linear Program (ILP)**, **Mixed-Integer Linear Program (MILP)**, or binary `0-1` optimization problem.
2. Requests **Branch and Bound** (`--method branch-and-bound` / `bnb`) with node-by-node tree trace, bounds, pruning reasons (`PRUNED_BOUND`, `INFEASIBLE`, `INTEGER_INCUMBENT`), and search strategy selection (`best-bound`, `dfs`, `bfs`).
3. Requests **Branch and Cut** (`--method branch-and-cut` / `bnc`), combining LP relaxation, cutting plane rounds (Chvátal-Gomory or GMI cuts via Dual Simplex), cut inheritance across subtrees, and branching.
4. Requests **Gomory Fractional Cuts** (`--method cutting-plane`) or **Gomory Mixed-Integer (GMI) Cuts** (`--cut-type gmi`) with step-by-step floor/fractional derivations.
5. Asks for an analysis of **Total Unimodularity (TUM)** on a constraint matrix or asks if a polyhedron is **Ideal** ($P = P_I$, Hoffman-Kruskal theorem).
6. Asks to visualize the continuous feasible region, integer lattice points $\mathbb{Z}^2$, branching hyperplanes, or cutting planes in 2D (`--visualize-2d`).

---

## Core Engine CLI

The standalone deterministic Python engine is located at:
`skills/ilp-solver/scripts/solve_ilp.py`

### Standard Invocations

```bash
# 1. Pure Chvatal-Gomory Cutting Planes with 2D visualization
python3 skills/ilp-solver/scripts/solve_ilp.py \
  --lp skills/ilp-solver/examples/gomory_2d.lp \
  --visualize-2d cuts_visualization.html --verify-scipy

# 2. Branch and Bound with Best-Bound search and 2D tree visualizer
python3 skills/ilp-solver/scripts/solve_ilp.py \
  --example bnb --method branch-and-bound \
  --visualize-2d bnb_visualization.html --verify-scipy

# 3. Branch and Cut (CG/GMI cuts + Dual Simplex + Branching)
python3 skills/ilp-solver/scripts/solve_ilp.py \
  --example bnb --method branch-and-cut --verify-scipy

# 4. Mixed-Integer Linear Program (MILP) solved via Branch and Cut (GMI cuts)
python3 skills/ilp-solver/scripts/solve_ilp.py \
  --example milp --method branch-and-cut --output-md milp_report.md --verify-scipy

# 5. Analyze Total Unimodularity and Ideal Polyhedron on a TUM problem
python3 skills/ilp-solver/scripts/solve_ilp.py --example tum --verify-scipy
```

### CLI Options Reference

| Flag | Description |
|---|---|
| `--lp <text_or_path>` | LP/ILP/MILP formulation string (`;` or newline separated) or path to `.lp` file. Supports variable types (`integer`, `binary`, `continuous`). |
| `--example {gomory2d,tum,bnb,milp}` | Run a built-in benchmark problem. |
| `--method {cutting-plane,cg,branch-and-bound,bnb,branch-and-cut,bnc}` | Solving algorithm (default: `cutting-plane`). |
| `--strategy {best-bound,dfs,bfs}` | Node selection rule for Branch & Bound / Branch & Cut (default: `best-bound`). |
| `--branch-rule {most-fractional,first-fractional}` | Variable selection rule for branching (default: `most-fractional`). |
| `--cut-type {auto,cg,gmi}` | Cut formula (`auto` uses CG when integer-safe, GMI otherwise). |
| `--max-cuts <N>` | Max cuts in `cutting-plane` mode (default: `15`). |
| `--cuts-per-node <N>` | Max cutting-plane rounds per node in `branch-and-cut` mode (default: `2`). |
| `--cut-depth <N>` | Max tree depth at which cuts are generated in `branch-and-cut` (default: `5`). |
| `--max-nodes <N>` | Node exploration limit for tree search (default: `200`). |
| `--visualize-2d <path.html>` | Export interactive standalone 2D HTML/SVG visualizer (`cutting-plane`, `branch-and-bound`, or `branch-and-cut`). |
| `--output-md <path.md>` | Export full Markdown summary report. |
| `--format {text,json}` | Terminal stdout format (default: `text`). |
| `--verify-scipy` | Cross-verify optimal value/status with `scipy.optimize.milp` (HiGHS). |

---

## Step-by-Step Workflow for Agents

### Stage 1: Pre-Solve Total Unimodularity (TUM) Analysis
Always check whether cutting planes or branching are even necessary:
1. **Necessity Check**: Verify that every entry $a_{ij} \in \{-1, 0, 1\}$. If any entry is not, $A$ is not TUM.
2. **Sufficient Criteria**: Check Poincaré / Heller-Tompkins criterion, Consecutive Ones property, and Bareiss submatrix determinants $\in \{-1, 0, 1\}$.
3. **Ideal Polyhedron Determination (Hoffman-Kruskal)**:
   - If $A$ is TUM and $b \in \mathbb{Z}^m$, $P = \{x \ge 0 : Ax \le b\} = P_I$, so the root LP relaxation yields an integer optimal solution directly with zero cuts or branches.

### Stage 2: Root Continuous LP Relaxation (Step 0)
1. Standardize all constraints (`<=`, `>=`, `=`, binary upper bounds `<= 1`) and solve via exact Two-Phase Simplex with Bland's rule.
2. Report the continuous relaxation solution $x^*_{\text{LP}}$ and objective $z_{\text{LP}}$.
3. If all integer/binary decision variables are already integers, stop and report the optimal solution.

### Stage 3: Solve by Requested Exact Algorithm
- **Mode A: Cutting Planes (`--method cutting-plane`)**:
  - Pick fractional integer-constrained basic row $x_B + \sum_{j \in N} \bar{a}_{ij} x_j = \bar{b}_i$ with $f_0 = \bar{b}_i - \lfloor \bar{b}_i \rfloor \in (0, 1)$.
  - Derive CG cut ($\sum f_j x_j \ge f_0$) or GMI cut ($\sum g_j x_j \ge f_0$), re-project into decision variables ($\alpha^T x \le \beta$), append $-s_{\text{cut}} + \sum (-g_j) x_j = -f_0$, and re-optimize via Dual Simplex.
- **Mode B: Branch and Bound (`--method branch-and-bound`)**:
  - Maintain global incumbent $(x^*, z^*)$.
  - Pop next open leaf by `--strategy` (`best-bound`, `dfs`, `bfs`).
  - Report each node's branch domain ($l_j \le x_j \le u_j$), LP relaxation bound, and exact outcome: `INFEASIBLE`, `PRUNED_BOUND`, `INTEGER_INCUMBENT`, `INTEGER_FATHOMED`, or `BRANCHED` ($x_r \le \lfloor \bar{x}_r \rfloor \mid x_r \ge \lceil \bar{x}_r \rceil$).
- **Mode C: Branch and Cut (`--method branch-and-cut`)**:
  - Combine node LP relaxation with up to `--cuts-per-node` CG/GMI cut rounds (with Dual Simplex re-optimization and immediate `INFEASIBLE`/`PRUNED_BOUND` checks) before branching.
  - Propagate valid cuts cleanly down each subtree (root cuts globally, local cuts within subtree).

### Stage 4: Final Solution, Cross-Check & 2D Graphical Visualization
1. Report status, optimal objective $z^*$, solution vector, and integrality gap $|z_{\text{LP}} - z^*|$.
2. If $n = 2$ and visualization is requested or helpful, generate the interactive HTML plot with `--visualize-2d <filename>.html`.
3. Cross-verify against HiGHS using `--verify-scipy` when SciPy is installed.

---

## Directory Reference

- [solve_ilp.py](./scripts/solve_ilp.py): Core deterministic ILP/MILP Python solver & 2D visualizer.
- [gomory_cuts.md](./references/gomory_cuts.md): Chvátal-Gomory & Dual Simplex derivation reference.
- [branch_and_bound_cut.md](./references/branch_and_bound_cut.md): Branch & Bound, Branch & Cut, and GMI cut mathematical reference.
- [total_unimodularity.md](./references/total_unimodularity.md): TUM criteria & Hoffman-Kruskal reference.
- [examples/gomory_2d.lp](./examples/gomory_2d.lp), [examples/bnb_2d.lp](./examples/bnb_2d.lp), [examples/milp_2d.lp](./examples/milp_2d.lp), [examples/tum_network.lp](./examples/tum_network.lp): Benchmark inputs.

---

## Constraints & Style Rules

1. **Deterministic Calculations Only**: Never guess tableau fractions, node bounds, or cuts mentally; always invoke `solve_ilp.py`.
2. **No Emojis**: Do not use any emojis in prose, headers, or lists.
3. **Exact Rational Fractions**: Express fractional quantities using LaTeX (e.g., $\frac{103}{26}$, $\frac{165}{4}$, $\frac{55}{3}$) alongside decimal approximations.
4. **File Links**: Create clickable markdown links for all generated files and artifacts.
