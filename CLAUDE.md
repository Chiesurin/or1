# Claude Code Instructions: OR1 Simplex Tableau Solver

Welcome to `or1`, a pure-Python Linear Programming and Simplex Tableau solver with exact rational arithmetic, Bland's anti-cycling rule, Two-Phase simplex, an interactive Web GUI, and a dedicated AI Agent skill.

---

## Primary Directive for Simplex Tableaus & Linear Programs

> **CRITICAL RULE**: When asked to solve a Simplex tableau, pivot a matrix, perform ratio tests, or solve a Linear Program:
> - **DO NOT perform mental arithmetic or calculate matrix operations manually.** (LLMs frequently make fraction arithmetic errors on matrix pivots).
> - **ALWAYS invoke the deterministic solver engine via Bash**:
>   ```bash
>   python3 skills/tableau-solver/scripts/solve_tableau.py [OPTIONS]
>   ```
> - Use the solver's output to present the complete walkthrough, step-by-step calculations, row operations ($R_i \leftarrow R_i - c \cdot R_{\text{pivot}}$), and optimal solution to the user.

---

## CLI Commands & Usage

### 1. Solve a Tableau from a CSV File
```bash
python3 skills/tableau-solver/scripts/solve_tableau.py --csv csv/tableau.csv
```
*(Auto-resolves paths across `csv/`, project root, and `examples/`).*

### 2. Solve a Tableau Passed Inline as Prompt Text
```bash
python3 skills/tableau-solver/scripts/solve_tableau.py --raw "Basis,RHS,x1,x2\n-z,0,-3,-5\nx3,4,1,0\nx4,12,0,2\nx5,18,3,2"
```

### 3. Solve an Algebraic Linear Program
```bash
python3 skills/tableau-solver/scripts/solve_tableau.py --lp "min -3*x1 - 5*x2 s.t. x1 <= 4; 2*x2 <= 12; 3*x1 + 2*x2 <= 18; x1, x2 >= 0"
```

### 4. Force Two-Phase Simplex
```bash
python3 skills/tableau-solver/scripts/solve_tableau.py --csv csv/tableau.csv --two-phase
```

### 5. Generate a Complete Solution Markdown File (.md)
```bash
python3 skills/tableau-solver/scripts/solve_tableau.py --csv csv/tableau.csv --output-md solution.md
```

### 6. Launch the Interactive Web GUI
```bash
python3 source/gui.py
```
*(Launches zero-dependency HTTP server on `http://localhost:8000`).*

---

## Mathematical Conventions

- **Universal Minimization ($\min z$)**: All problems are formulated as minimization. Maximization problems $\max c^T x$ are converted via $\min -c^T x$.
- **Row 0 Identity**: Row 0 is labeled `-z` (or `-w` during Phase 1). The governing equation is $-z + \sum \bar{c}_j x_j = \text{RHS}_0$, meaning the true objective value is $z = -\text{RHS}_0$.
- **Bland's Smallest-Subscript Anti-Cycling Rule**:
  - *Entering*: Smallest variable index $j$ among candidate columns with $c_j < 0$.
  - *Leaving*: Smallest variable index $i$ among rows achieving minimum ratio $\theta = \text{RHS}_i / a_{ij}$ ($a_{ij} > 0$).
  - Prevents cycling on degenerate tableaus (verified on Beale's 1955 problem).
- **Exact Rational Arithmetic**: Uses Python's standard library `fractions.Fraction` for zero floating-point error.

---

## Repository Map

- `skills/tableau-solver/SKILL.md`: Detailed agent runbook and prompt instructions.
- `skills/tableau-solver/scripts/solve_tableau.py`: Core CLI calculation engine.
- `source/gui.py`: Interactive Web GUI server and HTTP API.
- `csv/`: Tableau datasets (`tableau.csv` = Beale 1955, `example_standard.csv` = 2D textbook).
- `INSTRUCTIONS.md`: Human guide explaining how the project works.
