# Branch and Bound & Branch and Cut Algorithms

This reference details the exact mathematical formulations, search strategies, pruning rules, and cutting-plane integration used by `solve_ilp.py` for **Branch and Bound** and **Branch and Cut**.

---

## 1. Branch and Bound (B&B)

Consider a Mixed-Integer or Pure Integer Linear Program:

$$
\max \; c^T x \quad \text{s.t.} \quad Ax \le b, \; x \ge 0, \; x_j \in \mathbb{Z} \;\; \forall j \in I
$$

(Minimization problems follow symmetrically by reversing the bound inequalities.)

### Node Subproblem & Partitioning
Starting from root node $N_0$ (the continuous LP relaxation):
- Maintain a global primal **incumbent** $(x^*, z^*)$ (the best integer-feasible solution discovered so far).
- At any active node $N_k$, solve its continuous LP relaxation (using Two-Phase Simplex with Bland's rule) subject to original constraints and path branching bounds $[l_j, u_j]$.

### Fathoming / Pruning Rules
Let $z_{\text{LP}}(N_k)$ be the optimal LP objective value at node $N_k$ with solution $\bar{x}$:

1. **Pruned by Infeasibility (`INFEASIBLE`)**: Node LP relaxation has no feasible point. The entire subtree contains no feasible integer solution.
2. **Pruned by Bound / Dominance (`PRUNED_BOUND`)**:
   - For maximization: $z_{\text{LP}}(N_k) \le z^*$
   - For minimization: $z_{\text{LP}}(N_k) \ge z^*$
   No integer solution in the subtree can strictly improve upon the current incumbent $z^*$.
3. **Pruned by Integrality (`INTEGER_INCUMBENT` / `INTEGER_FATHOMED`)**:
   - Every constrained variable satisfies $\bar{x}_j \in \mathbb{Z}$ for all $j \in I$.
   - If strictly better than $z^*$, update incumbent $z^* \leftarrow z_{\text{LP}}(N_k)$ and $x^* \leftarrow \bar{x}$ (`INTEGER_INCUMBENT`). Otherwise fathom (`INTEGER_FATHOMED`).

### Branching Rule
If $\bar{x}_r \notin \mathbb{Z}$ for some $r \in I$ and $N_k$ is not pruned:
- Partition the feasible region into two disjoint child subproblems:
  $$N_{\text{left}}: x_r \le \lfloor \bar{x}_r \rfloor \qquad \text{and} \qquad N_{\text{right}}: x_r \ge \lceil \bar{x}_r \rceil$$
- Any integer point satisfies either $x_r \le \lfloor \bar{x}_r \rfloor$ or $x_r \ge \lceil \bar{x}_r \rceil$, and the fractional relaxation point $\bar{x}$ violates both.
- Variable selection modes:
  - `most-fractional`: maximizes $\min(f_r, 1 - f_r)$ where $f_r = \bar{x}_r - \lfloor \bar{x}_r \rfloor$.
  - `first-fractional`: chooses the smallest index $r \in I$ with fractional value.

### Node Selection Strategies (`--strategy`)
- **`best-bound`**: Pops the open leaf with the strongest LP relaxation bound ($\max z_{\text{LP}}$ for max, $\min z_{\text{LP}}$ for min). Minimizes dual gap rapidly.
- **`dfs` (Depth-First Search)**: Pops the newest open leaf first, diving quickly to reach an early integer incumbent and prune siblings cleanly.
- **`bfs` (Breadth-First Search)**: Explores the tree level by level.

---

## 2. Branch and Cut (B&C)

**Branch and Cut** strengthens Branch and Bound by tightening node LP relaxations with valid cutting planes before branching:

1. At node $N_k$ (when $\text{depth}(N_k) \le \text{cut\_depth}$), solve its LP relaxation.
2. If fractional and not dominated by $z^*$, apply up to `cuts-per-node` cutting-plane rounds:
   - Derive either a **Chvátal-Gomory (CG)** fractional cut (when all non-basic variables in the generating row are guaranteed integer) or a **Gomory Mixed-Integer (GMI)** cut (when continuous variables or non-integer slacks are present):
     $$
     g_j =
     \begin{cases}
     f_j & j \in N_I, \; f_j \le f_0 \\
     \frac{f_0 (1 - f_j)}{1 - f_0} & j \in N_I, \; f_j > f_0 \\
     \bar{a}_{ij} & j \in N_C, \; \bar{a}_{ij} \ge 0 \\
     \frac{f_0 |\bar{a}_{ij}|}{1 - f_0} & j \in N_C, \; \bar{a}_{ij} < 0
     \end{cases}
     \qquad \sum_{j \in N} g_j x_j \ge f_0
     $$
   - Re-project the cut into original decision variables ($\alpha^T x \le \beta$).
   - Inject $-s_{\text{cut}} + \sum (-g_j) x_j = -f_0$ and re-optimize via **Dual Simplex**.
   - Check immediately whether Dual Simplex proved the node `INFEASIBLE` or tightened the bound below/above $z^*$ (`PRUNED_BOUND`).
3. **Cut Validity Inheritance**:
   - Cuts generated at Root $N_0$ are **globally valid** for the entire tree.
   - Cuts generated at node $N_k$ use $N_k$'s branching bounds and are **locally valid** for all descendants in $N_k$'s subtree.
4. Once cuts at $N_k$ either reach `cuts-per-node` or integrality holds, proceed with B&B integrality check or branching.
