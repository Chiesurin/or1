---
name: ilp-solver
description: Solves Integer Linear Programs (ILP) using the Chvátal-Gomory Cutting Plane Algorithm, performs Total Unimodularity (TUM) and Ideal Polyhedron analysis, and generates interactive 2D graphical visualizations of cuts. Use when the user asks to solve an ILP, analyze Total Unimodularity, derive Gomory or Chvátal cuts, inspect the integer hull, or visualize cutting planes on 2D polyhedra.
---

# Integer Linear Programming (ILP) Solver Skill

This skill provides deterministic mathematical solving and visual explanations for Integer Linear Programming problems. It implements the **Chvátal-Gomory Cutting Plane Algorithm**, pre-solve **Total Unimodularity (TUM)** analysis for **Ideal Polyhedra**, and automatic **2D SVG/HTML graphical rendering**.

---

## Capabilities & Triggers

Activate this skill when the user:
1. Asks to solve an **Integer Linear Program (ILP)**, Mixed-Integer Linear Program, or pure integer optimization problem.
2. Requests **Gomory Fractional Cuts** or **Chvátal-Gomory Cuts** with step-by-step derivations.
3. Asks for an analysis of **Total Unimodularity (TUM)** on a constraint matrix or asks if a polyhedron is **Ideal** ($P = P_I$, Hoffman-Kruskal theorem).
4. Asks to visualize the continuous feasible region, integer lattice points $\mathbb{Z}^2$, and cutting planes in 2D.
5. Asks to trace the **Dual Simplex** re-optimization iterations after cut injection.

---

## Core Engine CLI

The standalone deterministic Python engine is located at:
`skills/ilp-solver/scripts/solve_ilp.py`

### Standard Invocations

```bash
# 1. Solve an ILP formulation directly from prompt text
python3 skills/ilp-solver/scripts/solve_ilp.py --lp "max 2*x1 + 3*x2 s.t. 2*x1 + 5*x2 <= 16; 6*x1 + 2*x2 <= 27; x1, x2 >= 0"

# 2. Solve from an .lp file and generate 2D interactive HTML visualization
python3 skills/ilp-solver/scripts/solve_ilp.py --lp skills/ilp-solver/examples/gomory_2d.lp --visualize-2d cuts_visualization.html

# 3. Analyze Total Unimodularity and Ideal Polyhedron on a TUM problem
python3 skills/ilp-solver/scripts/solve_ilp.py --example tum

# 4. Export Markdown report and cross-verify with SciPy HiGHS MILP solver
python3 skills/ilp-solver/scripts/solve_ilp.py --example gomory2d --output-md report.md --verify-scipy
```

---

## Step-by-Step Workflow for Agents

When answering a user request, execute the solver script and structure your explanation according to the following stages:

### Stage 1: Pre-Solve Total Unimodularity (TUM) Analysis
Always check whether cutting planes are even necessary:
1. **Necessity Check**: Verify that every entry $a_{ij} \in \{-1, 0, 1\}$. If any entry is not, $A$ is not TUM.
2. **Sufficient Criteria**:
   - Check Poincaré / Heller-Tompkins criterion for network flows and bipartite incidence.
   - Check Consecutive Ones (Interval Matrix) property.
   - For small matrices, check that all square submatrix determinants $\in \{-1, 0, 1\}$.
3. **Ideal Polyhedron Determination (Hoffman-Kruskal)**:
   - If $A$ is TUM and $b$ is integer, state clearly:
     $$P = \{x \ge 0 : Ax \le b\} \implies P = P_I = \text{conv}(P \cap \mathbb{Z}^n)$$
   - The continuous LP relaxation is guaranteed to have integer vertices. **No cuts are needed**, and the continuous Simplex solution directly solves the ILP.

### Stage 2: Continuous LP Relaxation (Step 0)
1. Display the initial standardized tableau with slack variables.
2. Report the optimal continuous solution $x^*_{\text{LP}}$ and objective value $z_{\text{LP}}$.
3. If $x^*_{\text{LP}} \in \mathbb{Z}^n$, stop and report the optimal integer solution.
4. If $x^*_{\text{LP}} \notin \mathbb{Z}^n$, identify all basic decision variables with fractional RHS values.

### Stage 3: Chvátal-Gomory Cut Derivations
For each cut iteration:
1. **Generating Row Selection**: State the chosen basic variable $x_B$ and its tableau row:
   $$x_B + \sum_{j \in N} \bar{a}_{ij} x_j = \bar{b}_i$$
2. **Floor Decomposition**: Show the exact integer floor and fractional part for every coefficient:
   $$\lfloor \bar{b}_i \rfloor, \quad f_0 = \bar{b}_i - \lfloor \bar{b}_i \rfloor \in (0, 1)$$
   $$\lfloor \bar{a}_{ij} \rfloor, \quad f_j = \bar{a}_{ij} - \lfloor \bar{a}_{ij} \rfloor \in [0, 1)$$
3. **Gomory Fractional Cut in Tableau Form**:
   $$\sum_{j \in N} f_j x_j \ge f_0 \iff -s_{\text{cut}} + \sum_{j \in N} (-f_j) x_j = -f_0$$
4. **Re-projected Cut in Original Decision Variables**:
   Substitute slack variable definitions to express the cut directly in decision variables:
   $$\alpha_1 x_1 + \alpha_2 x_2 + \dots \le \beta$$
   Explain how this cut strictly separates $x^*_{\text{LP}}$ while preserving all integer lattice points $\mathbb{Z}^n$.
5. **Dual Simplex Re-optimization**:
   - State leaving variable: $s_{\text{cut}}$ with negative RHS $-f_0$.
   - Execute the Dual Minimum Ratio Test across $a_{\text{cut}, j} < 0$:
     $$\theta^* = \min_{j : a_{\text{cut}, j} < 0} \frac{c_j}{|a_{\text{cut}, j}|}$$
   - State entering variable and show the updated tableau.

### Stage 4: Final Solution & 2D Graphical Visualization
1. Report the final integer optimal solution vector $x^*_{\text{ILP}}$ and optimal objective value $z^*_{\text{ILP}}$.
2. If the problem is 2D ($n=2$), generate the interactive HTML visualizer with `--visualize-2d <filename>.html`.
3. Provide a clear summary comparing continuous relaxation vs. integer optimum (integrality gap $z_{\text{LP}} - z^*_{\text{ILP}}$).

---

## Constraints & Style Rules

1. **Deterministic Calculations Only**: Never guess tableau fractions or cuts mentally; always run `solve_ilp.py` to extract exact mathematical values.
2. **No Emojis**: Do not use any emojis in prose, headers, or lists.
3. **Exact Rational Fractions**: Express fractional quantities using LaTeX (e.g. $\frac{103}{26}$, $\frac{21}{13}$, $\frac{166}{13}$) alongside their decimal approximations.
4. **File Links**: Create clickable markdown links for all generated files and artifacts.
