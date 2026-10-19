# Simplex Method with Bland's Rule in Python (CLI & Graphical Interfaces)

A complete, exact implementation of the **Simplex Algorithm** using **Bland's Smallest-Subscript Rule** to prevent cycling and guarantee convergence to optimality.

Includes both a **Terminal CLI** and modern interactive **Graphical User Interfaces (Web & Native Desktop)**.

---

## 🚀 Key Features

- **Exact Fraction Arithmetic**: Uses Python's standard `fractions.Fraction` — zero float rounding errors (e.g. `1/3`, `-11/2`).
- **Bland's Smallest-Subscript Anti-Cycling Rule**:
  - *Entering variable*: Smallest index among non-basic variables with negative reduced cost in Row 0.
  - *Leaving variable*: Smallest basic variable index among rows tied for the minimum non-negative ratio.
- **Two-Phase Simplex**:
  - Automatically handles general constraints ($\le, \ge, =$).
  - Artificial variables are named **$y_1, y_2, \dots$**.
  - Solves Phase 1 ($\min W = \sum y_i$), eliminates artificial variables, and seamlessly transitions to Phase 2.
- **CSV Tableau Import**:
  - Loads initial tableaus directly from CSV files in the folder (e.g., `tableau.csv`).
- **Zero Dependencies**:
  - Pure Python standard library (`http.server`, `tkinter`, `fractions`, `csv`, `urllib`). No `pip install` required!

---

## ❓ What is RHS?

**RHS** stands for **"Right-Hand Side"**.

### In Linear Programming Equations:
In standard linear programming notation:
$$\max \mathbf{c}^T \mathbf{x} \quad \text{subject to} \quad A\mathbf{x} \le \mathbf{b}, \quad \mathbf{x} \ge 0$$
The vector $\mathbf{b}$ represents the constant numerical values on the right-hand side of each constraint inequality or equality.

### In the Simplex Tableau:
The **RHS column** is the rightmost column of the simplex tableau:
1. **Row 0 (Objective Row)**:
   - Tracks the cost equation: **$-z = \text{RHS}_0$** (or $-w = \text{RHS}_0$ in Phase 1).
   - In the tableau, the basic variable is labeled **`-z`**, and the RHS shows the exact raw value. To find the minimum value of $z$, remember to change the sign: $z = -\text{RHS}_0$.
2. **Rows $1, \dots, m$ (Constraint Rows)**:
   - Each row corresponds to a basic variable. The entry in the **RHS column** indicates the current numerical value taken by that basic variable:
     $$x_{B(i)} = \text{RHS}_i$$
     (All non-basic variables currently equal $0$).
3. **Feasibility**:
   - For a basic solution to be feasible, all constraint RHS values must be non-negative: $\text{RHS}_i \ge 0$.
4. **Minimum Ratio Test**:
   - When choosing which variable leaves the basis, we compute:
     $$\theta_i = \frac{\text{RHS}_i}{a_{ij}} \quad \text{for all } a_{ij} > 0$$
     where $a_{ij}$ is the pivot column coefficient. The row with the minimum non-negative ratio $\theta$ leaves the basis.

---

## 🖥️ Graphical User Interface (`gui.py`)

A modern, responsive, standalone web application served locally using Python's standard `http.server`. Opens automatically in your default browser.

**To start:**
```bash
python3 gui.py
```
*(CLI mode is also supported: `python3 gui.py tableau.csv --cli`)*

### 🏠 Landing Screen: 3 Problem Options
When you first open the application, you can choose from three workflows:
1. **📁 Choose CSV File & Run Simplex**: Select a CSV file from the workspace (e.g. `tableau.csv`), upload your own, or preview & edit the tableau before solving.
2. **💡 Run Predefined Benchmarks**: 1-click execution for Beale's cycling benchmark, Two-Phase problem, and Standard 2D problem.
3. **✍️ Input LP Manually**: Formulate an optimization problem from scratch with custom variables, relations ($\le, \ge, =$), and MIN/MAX goal.

---

### 🔍 Sub-Step Elementary Operations & Focused Visualizer
For each iteration of the Simplex Method, the algorithm displays every discrete elementary row operation as an interactive sub-step:
1. **Pivot Selection (Bland's Rule & Minimum $\theta$ Test)**:
   - Selects candidate entering variable with smallest index among negative reduced costs.
   - Evaluates minimum $\theta$ test formatted with bracketed fractions: `[RHS] / [Coeff]` to select the leaving variable.
   - Tie-breaks resolved using Bland's smallest subscript rule on basic variable indices.
2. **Normalize Pivot Row ($R_p \leftarrow R_p / \text{pivot}$)**:
   - Scales the pivot row so the pivot entry becomes a leading 1.
   - **Selective Color Focus**: Inactive rows are dimmed to eliminate distraction, focusing exclusively on the normalized row $R_p$ and the gold pivot element.
3. **Row-by-Row Gauss-Jordan Elimination ($R_i \leftarrow R_i - \text{factor} \cdot R_p$)**:
   - Performs individual elementary row additions/subtractions on each other row.
   - **Selective Color Focus**: All unaffected rows are dimmed. The target row $R_i$ and reference pivot row $R_p$ are illuminated, with the newly zeroed entry highlighted in bright green (`cell-eliminated`).
4. **Iteration Complete (Updated Tableau & Basis)**:
   - Updates the basis labels and presents the full tableau in the clean 1-color-per-cell semantic color scheme.

---

### 🎨 Clean Neutral Tableau & Focused Operation Highlights
To eliminate visual noise and distraction, all static background colors have been removed from data cells:
- **Variable Names (Basis Status)**:
  - 🔴 **Basic Variables**: Displayed in **Red** with a `Basic` tag in the column header and red badges in the basis column.
  - 🔵 **Non-Basic Variables**: Displayed in **Blue** with a `Non-Basic` tag in the column header.
- **Data Cells (Numbers)**:
  - Default state (initial, completed iteration, optimal): Clean white/neutral cells with crisp numbers and standard borders. No rainbow of background colors.
- **Active Step Highlights ONLY**:
  - 🌟 **Pivot Element (Gold)**: Chosen pivot entry $a_{p,s}$ highlighted with a gold outline.
  - 🟦 **Target / Active Row (Blue / Amber)**: Shows which row is currently being scaled or eliminated.
  - 🟩 **Eliminated Entry (Green)**: The entry in the pivot column driven to $0$ by Gauss-Jordan elimination.
  - 🌫️ **Inactive Rows Dimmed**: Rows not participating in the active row operation are softly dimmed to direct the user's attention.
- **Stepper Bar Controls**: Navigate forward/backward by individual **Sub-steps** or jump directly between full **Iterations** (`⏮ Iter`, `◀ Sub-step`, `Sub-step ▶`, `Iter ⏭`).

---

## 📁 File Structure

| File | Description |
|---|---|
| `gui.py` | Complete standalone Simplex Method Web Application & Engine |
| `tableau.csv` | Default CSV file containing **Beale's Cycling Benchmark (1955)** |
| `example_standard.csv` | Standard 2D textbook linear program |
| `README.md` | Documentation |

---

## 📐 CSV Tableau Format

The standard tableau layout places the **RHS column on the left** (immediately following `Basis`), with **no dedicated $z$ column**:

```csv
Basis,RHS,x1,x2,x3,x4,s1,s2,s3
z,0,-10,57,9,24,0,0,0
s1,0,1/2,-11/2,-5/2,9,1,0,0
s2,0,1/2,-3/2,-1/2,1,0,1,0
s3,1,1,0,0,0,0,0,1
```

- **Column ordering**: `Basis, RHS, x1, x2, ...`
- **Cost row**: Row 0 labeled with the objective variable (`z` or `w`).
- **Constraint rows**: Labeled with current basic variables (`s1, s2, ...` or `y1, y2, ...`).
- **Numbers**: Exact rational fractions (`1/2`, `-11/2`), decimals (`0.5` or European comma `0,5`), and integers.
- **Delimiters**: Comma (`,`), semicolon (`;`), and tab (`\t`).
- *(Backwards compatibility)*: If an older CSV with RHS at the end or an explicit $z$ column is loaded, it is automatically parsed and displayed with RHS on the left and $z$ stripped.
