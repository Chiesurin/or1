# Total Unimodularity (TUM) and Ideal Polyhedra

This document outlines the theoretical foundations, necessary conditions, sufficient criteria, and algorithmic tests for Total Unimodularity and Ideal Polyhedra in Integer Linear Programming (ILP).

---

## 1. Definition

A matrix $A \in \mathbb{Z}^{m \times n}$ is **Totally Unimodular (TUM)** if every square submatrix (minor) of $A$ has determinant in $\{ -1, 0, 1 \}$.

> [!NOTE]
> This applies to square submatrices of **all** orders $k \times k$ where $1 \le k \le \min(m, n)$, not only full-rank submatrices.

---

## 2. The Hoffman-Kruskal Theorem (1956)

The fundamental theorem connecting Total Unimodularity to Integer Programming states:

Let $A \in \mathbb{Z}^{m \times n}$ be an integer matrix. Then:
1. $A$ is Totally Unimodular (TUM) if and only if for every integer vector $b \in \mathbb{Z}^m$, the polyhedron:
   $$P = \{ x \in \mathbb{R}^n : Ax \le b, \; x \ge 0 \}$$
   has only integer vertices (extreme points).
2. Consequently, the convex hull of integer points $P_I = \text{conv}(P \cap \mathbb{Z}^n)$ equals $P$ itself ($P = P_I$).
   Such a polyhedron is called an **Ideal Polyhedron**.

### Practical Implication for Solvers
When $A$ is TUM and $b$ is integer:
- The continuous LP relaxation solved by the standard Simplex algorithm **guarantees** an integer optimal solution.
- Cutting planes (Chvátal-Gomory) and Branch-and-Bound are completely unnecessary.

---

## 3. Necessary Conditions

Before computing submatrix determinants, test these fast necessity criteria:
1. **Entry Range**: Every entry $a_{ij}$ must satisfy $a_{ij} \in \{ -1, 0, 1 \}$.
   - If any entry $|a_{ij}| \ge 2$ or is fractional, $A$ cannot be TUM.
2. **Column/Row Scalability**: Any scalar multiplication by non-unit values destroys TUM.

---

## 4. Sufficient Conditions for TUM

Several structural properties guarantee that $A$ is TUM without evaluating minors:

### A. Poincaré / Heller-Tompkins Bipartite Criterion (1910, 1956)
A matrix $A \in \{ -1, 0, 1 \}^{m \times n}$ is TUM if:
1. Each column contains at most two non-zero entries.
2. The rows can be partitioned into two disjoint subsets $R_1$ and $R_2$ ($R_1 \cup R_2 = \{1, \dots, m\}$) such that for each column $j$:
   - If the two non-zero entries have the **same sign**, one entry is in $R_1$ and the other in $R_2$.
   - If the two non-zero entries have **opposite signs**, both entries are in the same subset (either both in $R_1$ or both in $R_2$).

**Special Cases**:
- **Node-Arc Incidence Matrix of Directed Graphs**: Every column has exactly one $+1$ and one $-1$. Letting $R_1 = \{1, \dots, m\}$ and $R_2 = \emptyset$ satisfies the criterion. Hence, network flow constraint matrices are always TUM.
- **Bipartite Graph Incidence Matrix**: Undirected graphs are bipartite if and only if their node-edge incidence matrix is TUM.

### B. Consecutive Ones Property (Interval Matrix)
If in every column (or every row) of a $\{0, 1\}$-matrix, all the $1$s appear in consecutive positions, then $A$ is Totally Unimodular.

### C. Ghouila-Houri Characterization (1962)
A matrix $A$ is TUM if and only if for every subset of rows $R \subseteq \{1, \dots, m\}$, there exists a partition of $R$ into $R_1$ and $R_2$ such that for every column $j$:
$$\sum_{i \in R_1} a_{ij} - \sum_{i \in R_2} a_{ij} \in \{ -1, 0, 1 \}$$

---

## 5. Algorithmic Workflow in the Skill

The `solve_ilp.py` engine performs pre-solve TUM analysis in order:
1. **Necessity Pass**: Verifies $a_{ij} \in \{-1, 0, 1\}$.
2. **Network Flow / Node-Arc Pass**: Tests column degree $\le 2$ and signed row partitions.
3. **Consecutive Ones Pass**: Tests contiguous ones block per column.
4. **Exhaustive Determinant Pass**: For matrices with $\min(m, n) \le 6$, evaluates determinants of all square submatrices via exact Bareiss integer elimination.
5. **Hoffman-Kruskal Check**: Cross-references TUM status against vector $b \in \mathbb{Z}^m$ to classify whether $P$ is an Ideal Polyhedron.
