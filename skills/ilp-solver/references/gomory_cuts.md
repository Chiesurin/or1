# Chvátal-Gomory Cutting Plane Algorithm

This document details the mathematical derivation of Gomory Fractional Cuts, Chvátal closures, and the Dual Simplex re-optimization mechanism.

---

## 1. Context and Motivation

Consider the pure Integer Linear Program (ILP):
$$\max \; c^T x \quad \text{subject to} \quad Ax \le b, \; x \ge 0, \; x \in \mathbb{Z}^n$$

The continuous LP relaxation drops the integrality constraint:
$$P = \{ x \in \mathbb{R}^n : Ax \le b, \; x \ge 0 \}$$

If the optimal continuous vertex $x^*_{\text{LP}}$ is integer, $x^*_{\text{LP}}$ is optimal for the ILP. If $x^*_{\text{LP}}$ is fractional, cutting planes are hyperplanes that:
1. Separate $x^*_{\text{LP}}$ from the integer hull ($x^*_{\text{LP}}$ violates the cut).
2. Do not exclude any feasible integer point ($P \cap \mathbb{Z}^n$ satisfies the cut).

---

## 2. Derivation of the Gomory Fractional Cut

Let the optimal LP relaxation tableau be in canonical form. Choose a basic variable $x_B$ whose current RHS value $\bar{b}_i$ is fractional:
$$\bar{b}_i \notin \mathbb{Z}$$

The corresponding tableau row equation is:
$$x_B + \sum_{j \in N} \bar{a}_{ij} x_j = \bar{b}_i$$

where $N$ denotes the set of non-basic variable indices.

### Step 1: Floor Decomposition
Decompose each coefficient into its integer floor and fractional part:
$$\bar{a}_{ij} = \lfloor \bar{a}_{ij} \rfloor + f_j, \quad \text{with } 0 \le f_j < 1$$
$$\bar{b}_i = \lfloor \bar{b}_i \rfloor + f_0, \quad \text{with } 0 < f_0 < 1$$

where $f_0 = \bar{b}_i - \lfloor \bar{b}_i \rfloor$ and $f_j = \bar{a}_{ij} - \lfloor \bar{a}_{ij} \rfloor$.

### Step 2: Algebraic Rearrangement
Substitute the decompositions into the row equation:
$$x_B + \sum_{j \in N} (\lfloor \bar{a}_{ij} \rfloor + f_j) x_j = \lfloor \bar{b}_i \rfloor + f_0$$
$$x_B + \sum_{j \in N} \lfloor \bar{a}_{ij} \rfloor x_j - \lfloor \bar{b}_i \rfloor = f_0 - \sum_{j \in N} f_j x_j$$

### Step 3: Integrality Argument
For any integer feasible solution $x \in \mathbb{Z}^n \ge 0$:
- The left-hand side is an integer because $x_B$, $x_j$, $\lfloor \bar{a}_{ij} \rfloor$, and $\lfloor \bar{b}_i \rfloor$ are all integers.
- Therefore, the right-hand side must also be an integer:
  $$f_0 - \sum_{j \in N} f_j x_j \in \mathbb{Z}$$
- Since $x_j \ge 0$ and $f_j \ge 0$, the sum $\sum_{j \in N} f_j x_j \ge 0$.
- Consequently:
  $$f_0 - \sum_{j \in N} f_j x_j \le f_0 < 1$$
- The largest integer strictly strictly less than 1 is at most 0:
  $$f_0 - \sum_{j \in N} f_j x_j \le 0 \iff \sum_{j \in N} f_j x_j \ge f_0$$

This inequality is the **Gomory Fractional Cut** (also known as a **Chvátal Cut**).

---

## 3. Representation in the Simplex Tableau

To inject the cut directly into the existing Simplex tableau:
1. Introduce a non-negative cut slack variable $s_{\text{cut}} \ge 0$:
   $$\sum_{j \in N} f_j x_j - s_{\text{cut}} = f_0 \iff -s_{\text{cut}} + \sum_{j \in N} (-f_j) x_j = -f_0$$
2. In canonical basic form with $s_{\text{cut}}$ as the new basic variable:
   - Row $R_{\text{cut}}$ has $1$ in the $s_{\text{cut}}$ column and $0$ in other basic columns.
   - For each non-basic variable $j \in N$, the coefficient is $-f_j \le 0$.
   - The RHS value is $-f_0 < 0$.

Because Row 0 costs remain $\ge 0$ (dual feasibility is maintained) while the RHS entry is negative (primal infeasibility), the tableau is ideally positioned for immediate re-optimization via **Dual Simplex**.

---

## 4. Dual Simplex Re-optimization Step

1. **Leaving Variable**: The newly added cut row $R_{\text{cut}}$ has RHS $=-f_0 < 0$, so $s_{\text{cut}}$ leaves the basis.
2. **Dual Minimum Ratio Test**: Across non-basic columns with negative entries in the cut row ($a_{\text{cut}, j} < 0$):
   $$\theta^* = \min_{j : a_{\text{cut}, j} < 0} \frac{c_j}{|a_{\text{cut}, j}|}$$
3. **Entering Variable**: Column $j^*$ that achieves the minimum ratio enters the basis.
4. **Elementary Row Operations**: Pivot on $(R_{\text{cut}}, j^*)$ to restore primal feasibility.

---

## 5. Re-projecting Cuts to Decision Variables

To plot cuts in 2D or express them in the original problem statement:
- Each slack variable $s_i$ was originally defined by an initial constraint $s_i = b_i - \sum_j A_{ij} x_j$.
- Substituting this relationship (and definitions of earlier cuts) into $\sum f_j x_j \ge f_0$ yields the cut directly in decision variables:
  $$\alpha_1 x_1 + \alpha_2 x_2 + \dots + \alpha_n x_n \le \beta$$
- This explicit boundary line is used by the 2D visualizer to display the slicing hyperplane.
