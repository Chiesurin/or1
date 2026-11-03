# Simplex Tableau Solution: sample_tableau.csv

> [!NOTE]
> Calculations performed using exact rational arithmetic (`fractions.Fraction`) with Bland's Smallest-Subscript Anti-Cycling Rule.
> **Objective Formulation**: Universal Minimization ($\min z$). Row 0 Identity: $-z = \text{RHS}_0 \implies z^* = -\text{RHS}_0$.

## Simplex Walkthrough

### Iteration 0 (Simplex)

#### Tableau State:
| Basis | RHS | x1 | x2 | x3 | x4 | x5 |
| --- | --- | --- | --- | --- | --- | --- |
| -z | 0 | -3 | -5 | 0 | 0 | 0 |
| x3 | 4 | 1 | 0 | 1 | 0 | 0 |
| x4 | 12 | 0 | 2 | 0 | 1 | 0 |
| x5 | 18 | 3 | 2 | 0 | 0 | 1 |

**Current Basic Solution**: x3 = 4, x4 = 12, x5 = 18
**Objective Equation**: $-z = 0 \implies z = 0$

#### 1. Entering Variable Selection (Bland's Rule)
- **Candidates with negative reduced costs ($c_j < 0$)**: `x1` (cost: $-3$, Bland index: 1), `x2` (cost: $-5$, Bland index: 2)
- **Selected Entering Variable**: `x1` (smallest Bland index).

#### 2. Minimum Ratio Test (Leaving Variable for Column `x1`)
| Row | Basis | RHS | Pivot Col Entry | Ratio $\theta = \text{RHS} / a_{ij}$ | Valid? |
| --- | --- | --- | --- | --- | --- |
| $R_1$ | `x3` | $4$ | $1$ | $\frac{4}{1} = 4$ | Yes |
| $R_2$ | `x4` | $12$ | $0$ | Skipped (coeff = 0) | No |
| $R_3$ | `x5` | $18$ | $3$ | $\frac{18}{3} = 6$ | Yes |

- **Minimum Ratio**: $\theta^* = 4$
- **Selected Leaving Variable**: `x3` ($R_1$)
- **Pivot Element**: $(R_1, \text{x1}) = 1$

#### 3. Elementary Row Operations (Pivot)
- **Normalize Pivot Row**: $R_1 \text{ already normalized (pivot = 1)}$
- **Row Eliminations**:
  - $-z \leftarrow -z + (3) \times R_1$
  - $R_3 \leftarrow R_3 - (3) \times R_1$
- **Basis Updated**: Row $R_1$ basis changes from `x3` $\rightarrow$ `x1`.

---

### Iteration 1 (Simplex)

#### Tableau State:
| Basis | RHS | x1 | x2 | x3 | x4 | x5 |
| --- | --- | --- | --- | --- | --- | --- |
| -z | 12 | 0 | -5 | 3 | 0 | 0 |
| x1 | 4 | 1 | 0 | 1 | 0 | 0 |
| x4 | 12 | 0 | 2 | 0 | 1 | 0 |
| x5 | 6 | 0 | 2 | -3 | 0 | 1 |

**Current Basic Solution**: x1 = 4, x4 = 12, x5 = 6
**Objective Equation**: $-z = 12 \implies z = -12$

#### 1. Entering Variable Selection (Bland's Rule)
- **Candidates with negative reduced costs ($c_j < 0$)**: `x2` (cost: $-5$, Bland index: 2)
- **Selected Entering Variable**: `x2` (smallest Bland index).

#### 2. Minimum Ratio Test (Leaving Variable for Column `x2`)
| Row | Basis | RHS | Pivot Col Entry | Ratio $\theta = \text{RHS} / a_{ij}$ | Valid? |
| --- | --- | --- | --- | --- | --- |
| $R_1$ | `x1` | $4$ | $0$ | Skipped (coeff = 0) | No |
| $R_2$ | `x4` | $12$ | $2$ | $\frac{12}{2} = 6$ | Yes |
| $R_3$ | `x5` | $6$ | $2$ | $\frac{6}{2} = 3$ | Yes |

- **Minimum Ratio**: $\theta^* = 3$
- **Selected Leaving Variable**: `x5` ($R_3$)
- **Pivot Element**: $(R_3, \text{x2}) = 2$

#### 3. Elementary Row Operations (Pivot)
- **Normalize Pivot Row**: $R_3 \leftarrow \frac{1}{2} \times R_3$
- **Row Eliminations**:
  - $-z \leftarrow -z + (5) \times R_3$
  - $R_2 \leftarrow R_2 - (2) \times R_3$
- **Basis Updated**: Row $R_3$ basis changes from `x5` $\rightarrow$ `x2`.

---

### Iteration 2 (Simplex)

#### Tableau State:
| Basis | RHS | x1 | x2 | x3 | x4 | x5 |
| --- | --- | --- | --- | --- | --- | --- |
| -z | 27 | 0 | 0 | -9/2 | 0 | 5/2 |
| x1 | 4 | 1 | 0 | 1 | 0 | 0 |
| x4 | 6 | 0 | 0 | 3 | 1 | -1 |
| x2 | 3 | 0 | 1 | -3/2 | 0 | 1/2 |

**Current Basic Solution**: x1 = 4, x4 = 6, x2 = 3
**Objective Equation**: $-z = 27 \implies z = -27$

#### 1. Entering Variable Selection (Bland's Rule)
- **Candidates with negative reduced costs ($c_j < 0$)**: `x3` (cost: $-9/2$, Bland index: 3)
- **Selected Entering Variable**: `x3` (smallest Bland index).

#### 2. Minimum Ratio Test (Leaving Variable for Column `x3`)
| Row | Basis | RHS | Pivot Col Entry | Ratio $\theta = \text{RHS} / a_{ij}$ | Valid? |
| --- | --- | --- | --- | --- | --- |
| $R_1$ | `x1` | $4$ | $1$ | $\frac{4}{1} = 4$ | Yes |
| $R_2$ | `x4` | $6$ | $3$ | $\frac{6}{3} = 2$ | Yes |
| $R_3$ | `x2` | $3$ | $-3/2$ | Skipped (coeff < 0) | No |

- **Minimum Ratio**: $\theta^* = 2$
- **Selected Leaving Variable**: `x4` ($R_2$)
- **Pivot Element**: $(R_2, \text{x3}) = 3$

#### 3. Elementary Row Operations (Pivot)
- **Normalize Pivot Row**: $R_2 \leftarrow \frac{1}{3} \times R_2$
- **Row Eliminations**:
  - $-z \leftarrow -z + (9/2) \times R_2$
  - $R_1 \leftarrow R_1 - R_2$
  - $R_3 \leftarrow R_3 + (3/2) \times R_2$
- **Basis Updated**: Row $R_2$ basis changes from `x4` $\rightarrow$ `x3`.

---

### Iteration 3 (Simplex)

#### Tableau State:
| Basis | RHS | x1 | x2 | x3 | x4 | x5 |
| --- | --- | --- | --- | --- | --- | --- |
| -z | 36 | 0 | 0 | 0 | 3/2 | 1 |
| x1 | 2 | 1 | 0 | 0 | -1/3 | 1/3 |
| x3 | 2 | 0 | 0 | 1 | 1/3 | -1/3 |
| x2 | 6 | 0 | 1 | 0 | 1/2 | 0 |

**Current Basic Solution**: x1 = 2, x3 = 2, x2 = 6
**Objective Equation**: $-z = 36 \implies z = -36$

> [!TIP]
> **Optimality Reached**: All reduced costs in Row 0 are $\ge 0$. The current tableau is optimal for Simplex.

## Final Solution Summary

- **Termination Status**: `OPTIMAL`
- **Optimal Value**: $\min z^* = -36$
- **Row 0 RHS**: $-z = 36$

### Optimal Variable Assignment:
| Variable | Type | Value |
| --- | --- | --- |
| `x1` | Basic | $2$ |
| `x2` | Basic | $6$ |
| `x3` | Basic | $2$ |
| `x4` | Non-Basic | $0$ |
| `x5` | Non-Basic | $0$ |


## Python Library Integration (SciPy `scipy.optimize.linprog`)

> [!NOTE]
> **SciPy Library Status**: SciPy is not installed in the active Python environment.
> **How to Install**:
> ```bash
> pip3 install --user scipy numpy
> ```
