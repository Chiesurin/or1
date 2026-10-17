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
   - Contains the current numerical value of the objective function ($Z$ or $-Z$ depending on sign formulation).
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

## 🖥️ Graphical User Interfaces

### 1. Web Browser GUI (`gui.py`) — Recommended ⭐
A modern, responsive web application served locally using Python's standard `http.server`. Opens automatically in your browser (Chrome, Safari, Firefox, Edge).

**To start:**
```bash
python3 gui.py
```
*(or press `[g]` from the `simplex.py` interactive menu).*

**Features:**
- 📊 **Step-by-Step Tableau Navigator**:
  - Stepper controls: `[⏮ First]`, `[◀ Prev]`, `[Next ▶]`, `[Last ⏭]`, `[▶ Auto-play]`.
  - **Color Highlights**:
    - **Pivot Cell** highlighted in bright amber/yellow.
    - **Entering Variable** column highlighted in soft blue.
    - **Leaving Variable** row highlighted in soft purple.
  - **Detailed Bland's Rule Explanations**:
    - Lists candidate entering variables with reduced costs and Bland indices.
    - Full Minimum Ratio Test table with division breakdown $\theta = \frac{\text{RHS}_i}{a_{ij}}$ and reasons for invalid rows ($a_{ij} \le 0$).
    - Automatic tie detection and Bland's tie-breaking explanation.
    - Exact elementary row operation formulas for each pivot.
  - **Optimal Solution Card**: Displays optimal objective value $Z^*$, basic variables, and alternative optima check.
- 📁 **CSV Manager & Interactive Editor**:
  - Dropdown list of CSV files found in the folder.
  - Direct file upload from disk.
  - **Live Editable Grid**: Modify numbers and fractions (`1/2`, `-3/2`), add rows, save back to disk, and solve immediately.
- ✍️ **Manual LP Form**:
  - Enter any linear program (Max/Min, arbitrary variables and constraints $\le, \ge, =$).
  - Automatically handles Two-Phase simplex with artificial variables $y_1, y_2, \dots$.
- 💡 **1-Click Benchmark Demos**:
  - Beale's cycling problem (1955), Standard 2D textbook problem, Two-Phase problem.

---

### 2. Native Desktop GUI (`gui_tk.py`)
For users who prefer a native OS desktop window:

```bash
python3 gui_tk.py
```

---

## ⌨️ Command-Line Interface (Terminal)

You can also run directly from the terminal:

```bash
# Direct run on a CSV file:
python3 simplex.py tableau.csv

# Interactive menu:
python3 simplex.py
```

---

## 📁 File Structure

| File | Description |
|---|---|
| `simplex.py` | Core Simplex engine with Bland's rule, Two-Phase method ($y_1, y_2$), fraction arithmetic, and CLI |
| `gui.py` | Web Graphical Interface (opens browser at `http://127.0.0.1:5050`) |
| `gui_tk.py` | Native Desktop Tkinter GUI window |
| `tableau.csv` | Default CSV file containing **Beale's Cycling Benchmark (1955)** |
| `example_standard.csv` | Standard 2D textbook linear program |
| `README.md` | This documentation |

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
