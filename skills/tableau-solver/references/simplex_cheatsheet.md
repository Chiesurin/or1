# Simplex Tableau Mathematical Cheatsheet

This reference provides the mathematical definitions, identities, and proofs implemented in `solve_tableau.py`.

---

## 1. Problem Formulation (Universal Minimization)

The solver standardizes every linear program to universal minimization:
$$\min z = \sum_{j=1}^n c_j x_j \quad \text{subject to} \quad A x = b, \quad x \ge 0, \quad b \ge 0$$

If a maximization problem is provided ($\max z_0 = \sum d_j x_j$), it is converted to minimization:
$$\min z = -\max z_0 = \sum (-d_j) x_j$$

---

## 2. Row 0 Convention & Identity

The Simplex Tableau objective row (Row 0) represents the equation:
$$-z + \sum_{j \in N} \bar{c}_j x_j = \text{RHS}_0$$

Where $\bar{c}_j$ are the reduced costs of non-basic variables $x_j$, and basic variables have $\bar{c}_{B_i} = 0$.

### Identity:
At any basic solution (where non-basic variables $x_j = 0$):
$$-z = \text{RHS}_0 \implies z = -\text{RHS}_0$$
When optimal:
$$\min z^* = -\text{RHS}_0$$

---

## 3. Bland's Smallest-Subscript Anti-Cycling Rule (1977)

Bland's rule guarantees finite termination and prevents cycling even on degenerate tableaus (e.g., Beale's counterexample).

### Variable Indexing:
Variables are ordered with 1-based indexing:
1. Decision variables: $x_1, x_2, \dots, x_n$
2. Slack variables: $x_{n+1}, x_{n+2}, \dots$ (or $s_1, s_2, \dots$)
3. Surplus variables: $e_1, e_2, \dots$
4. Artificial variables: $y_1, y_2, \dots$

### Pivot Selection:
1. **Entering Variable ($q$)**:
   $$\text{Candidate set } C = \{ j \mid \bar{c}_j < 0 \}$$
   If $C = \emptyset$, current tableau is **optimal**.
   Otherwise, select:
   $$q = \min \{ \text{index}(j) \mid j \in C \}$$
2. **Leaving Variable ($p$)**:
   Compute minimum ratio test for all rows $i \ge 1$ where $a_{iq} > 0$:
   $$\theta_i = \frac{\text{RHS}_i}{a_{iq}}, \quad \theta^* = \min_{i: a_{iq} > 0} \theta_i$$
   If no $a_{iq} > 0$, problem is **unbounded**.
   If multiple rows tie for $\theta^*$, apply Bland's rule:
   $$p = \operatorname{arg\,min} \{ \text{index}(x_{B_i}) \mid i \text{ such that } \theta_i = \theta^* \}$$

---

## 4. Elementary Row Operations (Gauss-Jordan)

Let $(p, q)$ be the chosen pivot element (row $p$, column $q$):

1. **Normalize Pivot Row**:
   $$R_p \leftarrow \frac{1}{a_{pq}} R_p$$
   (The pivot element becomes 1).
2. **Eliminate All Other Rows ($r \ne p$, including Row 0)**:
   $$R_r \leftarrow R_r - a_{rq} R_p$$
   (The entry in column $q$ of row $r$ becomes 0).
3. **Update Basis**:
   $$x_{B_p} \leftarrow x_q$$

---

## 5. Two-Phase Simplex Method

When the origin $x = 0$ is not feasible (e.g. equality constraints, $\ge$ constraints, or negative RHS values):

### Phase 1:
1. Normalize rows so that $\text{RHS}_i \ge 0$ (multiply rows with $\text{RHS}_i < 0$ by $-1$).
2. Add artificial variables $y_1, \dots, y_m \ge 0$ to form an initial identity basis.
3. Formulate Phase 1 objective:
   $$\min W = \sum_{i=1}^m y_i \quad \iff \quad -w = -\sum_{i=1}^m y_i$$
4. Substitute basic artificial variables to canonicalize Row 0:
   $$-w - \sum_{j} \left( \sum_{i=1}^m a_{ij} \right) x_j = -\sum_{i=1}^m \text{RHS}_i$$
5. Solve Phase 1 using Bland's rule:
   - If $\min W^* > 0$: Problem is **Infeasible**.
   - If $\min W^* = 0$: A feasible basic solution to the original problem has been found.

### Transition & Phase 2:
1. If artificial variables remain basic at value 0, pivot them out with non-artificial variables or drop redundant rows.
2. Remove artificial columns $y_1, \dots, y_m$.
3. Restore the original objective function row (Row 0: $-z$).
4. Canonicalize Row 0 against the current basis ($R_0 \leftarrow R_0 - c_{B_i} R_i$).
5. Solve Phase 2 using Bland's rule to find the optimal solution $z^*$.
