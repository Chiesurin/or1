# Simplex Tableau Solution: tableau.csv

> [!NOTE]
> Calculations performed using exact rational arithmetic (`fractions.Fraction`) with Bland's Smallest-Subscript Anti-Cycling Rule.
> **Objective Formulation**: Universal Minimization ($\min z$). Row 0 Identity: $-z = \text{RHS}_0 \implies z^* = -\text{RHS}_0$.

## Simplex Walkthrough

### Iteration 0 (Simplex)

#### Tableau State:
| Basis | RHS | x1 | x2 | x3 | x4 | x5 | x6 | x7 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -z | 0 | -10 | 57 | 9 | 24 | 0 | 0 | 0 |
| x5 | 0 | 1/2 | -11/2 | -5/2 | 9 | 1 | 0 | 0 |
| x6 | 0 | 1/2 | -3/2 | -1/2 | 1 | 0 | 1 | 0 |
| x7 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 |

**Current Basic Solution**: x5 = 0, x6 = 0, x7 = 1
**Objective Equation**: $-z = 0 \implies z = 0$

#### 1. Entering Variable Selection (Bland's Rule)
- **Candidates with negative reduced costs ($c_j < 0$)**: `x1` (cost: $-10$, Bland index: 1)
- **Selected Entering Variable**: `x1` (smallest Bland index).

#### 2. Minimum Ratio Test (Leaving Variable for Column `x1`)
| Row | Basis | RHS | Pivot Col Entry | Ratio $\theta = \text{RHS} / a_{ij}$ | Valid? |
| --- | --- | --- | --- | --- | --- |
| $R_1$ | `x5` | $0$ | $1/2$ | $\frac{0}{1/2} = 0$ | Yes |
| $R_2$ | `x6` | $0$ | $1/2$ | $\frac{0}{1/2} = 0$ | Yes |
| $R_3$ | `x7` | $1$ | $1$ | $\frac{1}{1} = 1$ | Yes |

- **Minimum Ratio**: $\theta^* = 0$
- **Selected Leaving Variable**: `x5` ($R_1$)
- > [!IMPORTANT]
  > Degenerate tie between rows [1, 2] (min θ = 0). Bland's Rule selects basic variable with smallest index: 'x5' (Index 5).
- **Pivot Element**: $(R_1, \text{x1}) = 1/2$

#### 3. Elementary Row Operations (Pivot)
- **Normalize Pivot Row**: $R_1 \leftarrow \frac{1}{1/2} \times R_1$
- **Row Eliminations**:
  - $-z \leftarrow -z + (10) \times R_1$
  - $R_2 \leftarrow R_2 - (1/2) \times R_1$
  - $R_3 \leftarrow R_3 - R_1$
- **Basis Updated**: Row $R_1$ basis changes from `x5` $\rightarrow$ `x1`.

---

### Iteration 1 (Simplex)

#### Tableau State:
| Basis | RHS | x1 | x2 | x3 | x4 | x5 | x6 | x7 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -z | 0 | 0 | -53 | -41 | 204 | 20 | 0 | 0 |
| x1 | 0 | 1 | -11 | -5 | 18 | 2 | 0 | 0 |
| x6 | 0 | 0 | 4 | 2 | -8 | -1 | 1 | 0 |
| x7 | 1 | 0 | 11 | 5 | -18 | -2 | 0 | 1 |

**Current Basic Solution**: x1 = 0, x6 = 0, x7 = 1
**Objective Equation**: $-z = 0 \implies z = 0$

#### 1. Entering Variable Selection (Bland's Rule)
- **Candidates with negative reduced costs ($c_j < 0$)**: `x2` (cost: $-53$, Bland index: 2), `x3` (cost: $-41$, Bland index: 3)
- **Selected Entering Variable**: `x2` (smallest Bland index).

#### 2. Minimum Ratio Test (Leaving Variable for Column `x2`)
| Row | Basis | RHS | Pivot Col Entry | Ratio $\theta = \text{RHS} / a_{ij}$ | Valid? |
| --- | --- | --- | --- | --- | --- |
| $R_1$ | `x1` | $0$ | $-11$ | Skipped (coeff < 0) | No |
| $R_2$ | `x6` | $0$ | $4$ | $\frac{0}{4} = 0$ | Yes |
| $R_3$ | `x7` | $1$ | $11$ | $\frac{1}{11} = 1/11$ | Yes |

- **Minimum Ratio**: $\theta^* = 0$
- **Selected Leaving Variable**: `x6` ($R_2$)
- **Pivot Element**: $(R_2, \text{x2}) = 4$

#### 3. Elementary Row Operations (Pivot)
- **Normalize Pivot Row**: $R_2 \leftarrow \frac{1}{4} \times R_2$
- **Row Eliminations**:
  - $-z \leftarrow -z + (53) \times R_2$
  - $R_1 \leftarrow R_1 + (11) \times R_2$
  - $R_3 \leftarrow R_3 - (11) \times R_2$
- **Basis Updated**: Row $R_2$ basis changes from `x6` $\rightarrow$ `x2`.

---

### Iteration 2 (Simplex)

#### Tableau State:
| Basis | RHS | x1 | x2 | x3 | x4 | x5 | x6 | x7 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -z | 0 | 0 | 0 | -29/2 | 98 | 27/4 | 53/4 | 0 |
| x1 | 0 | 1 | 0 | 1/2 | -4 | -3/4 | 11/4 | 0 |
| x2 | 0 | 0 | 1 | 1/2 | -2 | -1/4 | 1/4 | 0 |
| x7 | 1 | 0 | 0 | -1/2 | 4 | 3/4 | -11/4 | 1 |

**Current Basic Solution**: x1 = 0, x2 = 0, x7 = 1
**Objective Equation**: $-z = 0 \implies z = 0$

#### 1. Entering Variable Selection (Bland's Rule)
- **Candidates with negative reduced costs ($c_j < 0$)**: `x3` (cost: $-29/2$, Bland index: 3)
- **Selected Entering Variable**: `x3` (smallest Bland index).

#### 2. Minimum Ratio Test (Leaving Variable for Column `x3`)
| Row | Basis | RHS | Pivot Col Entry | Ratio $\theta = \text{RHS} / a_{ij}$ | Valid? |
| --- | --- | --- | --- | --- | --- |
| $R_1$ | `x1` | $0$ | $1/2$ | $\frac{0}{1/2} = 0$ | Yes |
| $R_2$ | `x2` | $0$ | $1/2$ | $\frac{0}{1/2} = 0$ | Yes |
| $R_3$ | `x7` | $1$ | $-1/2$ | Skipped (coeff < 0) | No |

- **Minimum Ratio**: $\theta^* = 0$
- **Selected Leaving Variable**: `x1` ($R_1$)
- > [!IMPORTANT]
  > Degenerate tie between rows [1, 2] (min θ = 0). Bland's Rule selects basic variable with smallest index: 'x1' (Index 1).
- **Pivot Element**: $(R_1, \text{x3}) = 1/2$

#### 3. Elementary Row Operations (Pivot)
- **Normalize Pivot Row**: $R_1 \leftarrow \frac{1}{1/2} \times R_1$
- **Row Eliminations**:
  - $-z \leftarrow -z + (29/2) \times R_1$
  - $R_2 \leftarrow R_2 - (1/2) \times R_1$
  - $R_3 \leftarrow R_3 + (1/2) \times R_1$
- **Basis Updated**: Row $R_1$ basis changes from `x1` $\rightarrow$ `x3`.

---

### Iteration 3 (Simplex)

#### Tableau State:
| Basis | RHS | x1 | x2 | x3 | x4 | x5 | x6 | x7 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -z | 0 | 29 | 0 | 0 | -18 | -15 | 93 | 0 |
| x3 | 0 | 2 | 0 | 1 | -8 | -3/2 | 11/2 | 0 |
| x2 | 0 | -1 | 1 | 0 | 2 | 1/2 | -5/2 | 0 |
| x7 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 |

**Current Basic Solution**: x3 = 0, x2 = 0, x7 = 1
**Objective Equation**: $-z = 0 \implies z = 0$

#### 1. Entering Variable Selection (Bland's Rule)
- **Candidates with negative reduced costs ($c_j < 0$)**: `x4` (cost: $-18$, Bland index: 4), `x5` (cost: $-15$, Bland index: 5)
- **Selected Entering Variable**: `x4` (smallest Bland index).

#### 2. Minimum Ratio Test (Leaving Variable for Column `x4`)
| Row | Basis | RHS | Pivot Col Entry | Ratio $\theta = \text{RHS} / a_{ij}$ | Valid? |
| --- | --- | --- | --- | --- | --- |
| $R_1$ | `x3` | $0$ | $-8$ | Skipped (coeff < 0) | No |
| $R_2$ | `x2` | $0$ | $2$ | $\frac{0}{2} = 0$ | Yes |
| $R_3$ | `x7` | $1$ | $0$ | Skipped (coeff = 0) | No |

- **Minimum Ratio**: $\theta^* = 0$
- **Selected Leaving Variable**: `x2` ($R_2$)
- **Pivot Element**: $(R_2, \text{x4}) = 2$

#### 3. Elementary Row Operations (Pivot)
- **Normalize Pivot Row**: $R_2 \leftarrow \frac{1}{2} \times R_2$
- **Row Eliminations**:
  - $-z \leftarrow -z + (18) \times R_2$
  - $R_1 \leftarrow R_1 + (8) \times R_2$
- **Basis Updated**: Row $R_2$ basis changes from `x2` $\rightarrow$ `x4`.

---

### Iteration 4 (Simplex)

#### Tableau State:
| Basis | RHS | x1 | x2 | x3 | x4 | x5 | x6 | x7 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -z | 0 | 20 | 9 | 0 | 0 | -21/2 | 141/2 | 0 |
| x3 | 0 | -2 | 4 | 1 | 0 | 1/2 | -9/2 | 0 |
| x4 | 0 | -1/2 | 1/2 | 0 | 1 | 1/4 | -5/4 | 0 |
| x7 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 |

**Current Basic Solution**: x3 = 0, x4 = 0, x7 = 1
**Objective Equation**: $-z = 0 \implies z = 0$

#### 1. Entering Variable Selection (Bland's Rule)
- **Candidates with negative reduced costs ($c_j < 0$)**: `x5` (cost: $-21/2$, Bland index: 5)
- **Selected Entering Variable**: `x5` (smallest Bland index).

#### 2. Minimum Ratio Test (Leaving Variable for Column `x5`)
| Row | Basis | RHS | Pivot Col Entry | Ratio $\theta = \text{RHS} / a_{ij}$ | Valid? |
| --- | --- | --- | --- | --- | --- |
| $R_1$ | `x3` | $0$ | $1/2$ | $\frac{0}{1/2} = 0$ | Yes |
| $R_2$ | `x4` | $0$ | $1/4$ | $\frac{0}{1/4} = 0$ | Yes |
| $R_3$ | `x7` | $1$ | $0$ | Skipped (coeff = 0) | No |

- **Minimum Ratio**: $\theta^* = 0$
- **Selected Leaving Variable**: `x3` ($R_1$)
- > [!IMPORTANT]
  > Degenerate tie between rows [1, 2] (min θ = 0). Bland's Rule selects basic variable with smallest index: 'x3' (Index 3).
- **Pivot Element**: $(R_1, \text{x5}) = 1/2$

#### 3. Elementary Row Operations (Pivot)
- **Normalize Pivot Row**: $R_1 \leftarrow \frac{1}{1/2} \times R_1$
- **Row Eliminations**:
  - $-z \leftarrow -z + (21/2) \times R_1$
  - $R_2 \leftarrow R_2 - (1/4) \times R_1$
- **Basis Updated**: Row $R_1$ basis changes from `x3` $\rightarrow$ `x5`.

---

### Iteration 5 (Simplex)

#### Tableau State:
| Basis | RHS | x1 | x2 | x3 | x4 | x5 | x6 | x7 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -z | 0 | -22 | 93 | 21 | 0 | 0 | -24 | 0 |
| x5 | 0 | -4 | 8 | 2 | 0 | 1 | -9 | 0 |
| x4 | 0 | 1/2 | -3/2 | -1/2 | 1 | 0 | 1 | 0 |
| x7 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 |

**Current Basic Solution**: x5 = 0, x4 = 0, x7 = 1
**Objective Equation**: $-z = 0 \implies z = 0$

#### 1. Entering Variable Selection (Bland's Rule)
- **Candidates with negative reduced costs ($c_j < 0$)**: `x1` (cost: $-22$, Bland index: 1), `x6` (cost: $-24$, Bland index: 6)
- **Selected Entering Variable**: `x1` (smallest Bland index).

#### 2. Minimum Ratio Test (Leaving Variable for Column `x1`)
| Row | Basis | RHS | Pivot Col Entry | Ratio $\theta = \text{RHS} / a_{ij}$ | Valid? |
| --- | --- | --- | --- | --- | --- |
| $R_1$ | `x5` | $0$ | $-4$ | Skipped (coeff < 0) | No |
| $R_2$ | `x4` | $0$ | $1/2$ | $\frac{0}{1/2} = 0$ | Yes |
| $R_3$ | `x7` | $1$ | $1$ | $\frac{1}{1} = 1$ | Yes |

- **Minimum Ratio**: $\theta^* = 0$
- **Selected Leaving Variable**: `x4` ($R_2$)
- **Pivot Element**: $(R_2, \text{x1}) = 1/2$

#### 3. Elementary Row Operations (Pivot)
- **Normalize Pivot Row**: $R_2 \leftarrow \frac{1}{1/2} \times R_2$
- **Row Eliminations**:
  - $-z \leftarrow -z + (22) \times R_2$
  - $R_1 \leftarrow R_1 + (4) \times R_2$
  - $R_3 \leftarrow R_3 - R_2$
- **Basis Updated**: Row $R_2$ basis changes from `x4` $\rightarrow$ `x1`.

---

### Iteration 6 (Simplex)

#### Tableau State:
| Basis | RHS | x1 | x2 | x3 | x4 | x5 | x6 | x7 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -z | 0 | 0 | 27 | -1 | 44 | 0 | 20 | 0 |
| x5 | 0 | 0 | -4 | -2 | 8 | 1 | -1 | 0 |
| x1 | 0 | 1 | -3 | -1 | 2 | 0 | 2 | 0 |
| x7 | 1 | 0 | 3 | 1 | -2 | 0 | -2 | 1 |

**Current Basic Solution**: x5 = 0, x1 = 0, x7 = 1
**Objective Equation**: $-z = 0 \implies z = 0$

#### 1. Entering Variable Selection (Bland's Rule)
- **Candidates with negative reduced costs ($c_j < 0$)**: `x3` (cost: $-1$, Bland index: 3)
- **Selected Entering Variable**: `x3` (smallest Bland index).

#### 2. Minimum Ratio Test (Leaving Variable for Column `x3`)
| Row | Basis | RHS | Pivot Col Entry | Ratio $\theta = \text{RHS} / a_{ij}$ | Valid? |
| --- | --- | --- | --- | --- | --- |
| $R_1$ | `x5` | $0$ | $-2$ | Skipped (coeff < 0) | No |
| $R_2$ | `x1` | $0$ | $-1$ | Skipped (coeff < 0) | No |
| $R_3$ | `x7` | $1$ | $1$ | $\frac{1}{1} = 1$ | Yes |

- **Minimum Ratio**: $\theta^* = 1$
- **Selected Leaving Variable**: `x7` ($R_3$)
- **Pivot Element**: $(R_3, \text{x3}) = 1$

#### 3. Elementary Row Operations (Pivot)
- **Normalize Pivot Row**: $R_3 \text{ already normalized (pivot = 1)}$
- **Row Eliminations**:
  - $-z \leftarrow -z + R_3$
  - $R_1 \leftarrow R_1 + (2) \times R_3$
  - $R_2 \leftarrow R_2 + R_3$
- **Basis Updated**: Row $R_3$ basis changes from `x7` $\rightarrow$ `x3`.

---

### Iteration 7 (Simplex)

#### Tableau State:
| Basis | RHS | x1 | x2 | x3 | x4 | x5 | x6 | x7 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| -z | 1 | 0 | 30 | 0 | 42 | 0 | 18 | 1 |
| x5 | 2 | 0 | 2 | 0 | 4 | 1 | -5 | 2 |
| x1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 1 |
| x3 | 1 | 0 | 3 | 1 | -2 | 0 | -2 | 1 |

**Current Basic Solution**: x5 = 2, x1 = 1, x3 = 1
**Objective Equation**: $-z = 1 \implies z = -1$

> [!TIP]
> **Optimality Reached**: All reduced costs in Row 0 are $\ge 0$. The current tableau is optimal for Simplex.

## Final Solution Summary

- **Termination Status**: `OPTIMAL`
- **Optimal Value**: $\min z^* = -1$
- **Row 0 RHS**: $-z = 1$

### Optimal Variable Assignment:
| Variable | Type | Value |
| --- | --- | --- |
| `x1` | Basic | $1$ |
| `x2` | Non-Basic | $0$ |
| `x3` | Basic | $1$ |
| `x4` | Non-Basic | $0$ |
| `x5` | Basic | $2$ |
| `x6` | Non-Basic | $0$ |
| `x7` | Non-Basic | $0$ |


## 🔬 Python Library Integration (SciPy `scipy.optimize.linprog`)

> [!NOTE]
> **SciPy Library Status**: SciPy is not installed in the active Python environment.
> **How to Install**:
> ```bash
> pip3 install --user scipy numpy
> ```
