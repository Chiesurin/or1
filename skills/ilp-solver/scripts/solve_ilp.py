#!/usr/bin/env python3
"""
Integer Linear Programming (ILP) & Mixed-Integer Linear Programming (MILP) Solver:
- Pre-solve Total Unimodularity (TUM) and Ideal Polyhedron analysis (Hoffman-Kruskal).
- Chvatal-Gomory and Gomory Mixed-Integer (GMI) Cutting Plane Algorithm.
- Branch and Bound (Best-Bound, DFS, BFS; Most-Fractional / First-Fractional branching).
- Branch and Cut (LP relaxation + CG/GMI cuts + Dual Simplex + Branching tree search).
- Interactive 2D SVG/HTML graphical visualization for Cutting Planes, Branch & Bound, and Branch & Cut.
- Optional cross-verification with scipy.optimize.milp (HiGHS).
- Zero mandatory external dependencies (fractions.Fraction pure standard library).
"""

import sys
import os
import re
import math
import json
import itertools
import argparse
from fractions import Fraction
from typing import List, Tuple, Dict, Optional, Set, Any


# ---------------------------------------------------------------------
# FRACTION & ARITHMETIC UTILITIES
# ---------------------------------------------------------------------

def parse_fraction(val_str: str) -> Fraction:
    """Parses a string into an exact Fraction. Supports integers, a/b, and decimals."""
    s = str(val_str).strip()
    if not s:
        raise ValueError("Cannot parse empty value.")
    if re.match(r'^-?\d+,\d+$', s):
        s = s.replace(',', '.')
    if '/' in s:
        parts = s.split('/')
        if len(parts) != 2:
            raise ValueError(f"Invalid fraction: '{s}'")
        return Fraction(int(parts[0].strip()), int(parts[1].strip()))
    return Fraction(s)


def format_fraction(f: Fraction) -> str:
    """Formats a Fraction as an integer if denominator is 1, else 'num/den'."""
    if f.denominator == 1:
        return str(f.numerator)
    return f"{f.numerator}/{f.denominator}"


def format_approx(f: Fraction) -> str:
    """Formats a Fraction with decimal approximation if fractional."""
    s = format_fraction(f)
    if f.denominator != 1:
        s += f" (~{float(f):.4g})"
    return s


def floor_frac(f: Fraction) -> int:
    """Exact floor function for Fraction (largest integer <= f)."""
    return f.numerator // f.denominator


def frac_part(f: Fraction) -> Fraction:
    """Exact fractional part of f: f - floor(f) in [0, 1)."""
    return f - Fraction(floor_frac(f), 1)


# ---------------------------------------------------------------------
# TOTAL UNIMODULARITY (TUM) & IDEAL POLYHEDRON ANALYZER
# ---------------------------------------------------------------------

def compute_determinant(matrix: List[List[int]]) -> int:
    """Computes exact determinant of an integer square matrix via Bareiss algorithm."""
    n = len(matrix)
    if n == 0:
        return 1
    if n == 1:
        return matrix[0][0]
    if n == 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]

    # Bareiss integer division algorithm
    mat = [row[:] for row in matrix]
    sign = 1
    prev = 1

    for k in range(n - 1):
        # Pivot selection
        if mat[k][k] == 0:
            swap_idx = -1
            for i in range(k + 1, n):
                if mat[i][k] != 0:
                    swap_idx = i
                    break
            if swap_idx == -1:
                return 0
            mat[k], mat[swap_idx] = mat[swap_idx], mat[k]
            sign = -sign

        pivot = mat[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = mat[i][j] * pivot - mat[i][k] * mat[k][j]
                mat[i][j] = num // prev
        prev = pivot

    return sign * mat[n - 1][n - 1]


def analyze_tum(A: List[List[Fraction]], b: List[Fraction], var_names: List[str]) -> Dict[str, Any]:
    """
    Performs rigorous Total Unimodularity (TUM) analysis on constraint matrix A:
    1. Tests necessity conditions (all entries in {-1, 0, 1}).
    2. Tests sufficient conditions (Poincare / Heller-Tompkins bipartite criterion, consecutive ones).
    3. For small-to-medium matrices (min(m, n) <= 8), exhaustively verifies all submatrix determinants.
    4. Evaluates Hoffman-Kruskal theorem to determine if polyhedron P = {x >= 0 : Ax <= b} is Ideal (P = P_I).
    """
    m = len(A)
    n = len(A[0]) if m > 0 else 0

    report = {
        "m": m,
        "n": n,
        "is_tum": False,
        "is_ideal": False,
        "necessity_passed": True,
        "non_tum_reason": None,
        "sufficient_method": None,
        "all_submatrices_checked": False,
        "details": []
    }

    if m == 0 or n == 0:
        report["is_tum"] = True
        report["is_ideal"] = True
        return report

    # 1. Necessity check: every entry must be in {-1, 0, 1}
    int_A: List[List[int]] = []
    for i in range(m):
        int_row = []
        for j in range(n):
            val = A[i][j]
            if val.denominator != 1 or val.numerator not in (-1, 0, 1):
                report["necessity_passed"] = False
                report["is_tum"] = False
                report["non_tum_reason"] = (
                    f"Entry A[{i+1},{j+1}] = {format_fraction(val)} is not in {{-1, 0, 1}}. "
                    "A matrix can only be Totally Unimodular if all its entries belong to {-1, 0, 1}."
                )
                report["details"].append(report["non_tum_reason"])
                # Even if not TUM, check if b is integer for reporting
                b_int = all(x.denominator == 1 for x in b)
                report["details"].append(f"Vector b is {'integral' if b_int else 'non-integral'}.")
                report["details"].append("Conclusion: Polyhedron is NOT guaranteed to be ideal. Cutting planes may be required.")
                return report
            int_row.append(val.numerator)
        int_A.append(int_row)

    report["details"].append("[Necessity Check PASSED] All entries of A belong to {-1, 0, 1}.")

    # 2. Check Heller-Tompkins / Poincare sufficient condition (bipartite network flow structure)
    # Condition: each column has at most two non-zeros.
    max_nonzeros_col = max(sum(1 for i in range(m) if int_A[i][j] != 0) for j in range(n))
    if max_nonzeros_col <= 2:
        # Check if rows can be partitioned into R1 and R2
        # For standard node-arc incidence: each column has one +1 and one -1 (or one +-1)
        is_node_arc = True
        for j in range(n):
            col_vals = [int_A[i][j] for i in range(m) if int_A[i][j] != 0]
            if len(col_vals) == 2 and sum(col_vals) != 0:
                is_node_arc = False
                break
        if is_node_arc:
            report["is_tum"] = True
            report["sufficient_method"] = "Network Flow / Node-Arc Incidence Criterion (Poincaré 1910)"
            report["details"].append(
                "[Sufficient Condition MET] Matrix A is the node-arc incidence matrix of a directed graph "
                "(each column has at most two non-zero entries, summing to 0). Hence, A is Totally Unimodular (TUM)."
            )

    # 3. Check Interval / Consecutive Ones Property
    if not report["is_tum"]:
        consecutive_ones = True
        for j in range(n):
            ones_rows = [i for i in range(m) if int_A[i][j] == 1]
            non_ones = [i for i in range(m) if int_A[i][j] != 0 and int_A[i][j] != 1]
            if non_ones:
                consecutive_ones = False
                break
            if len(ones_rows) > 1:
                if max(ones_rows) - min(ones_rows) + 1 != len(ones_rows):
                    consecutive_ones = False
                    break
        if consecutive_ones:
            report["is_tum"] = True
            report["sufficient_method"] = "Consecutive Ones Property (Interval Matrix)"
            report["details"].append(
                "[Sufficient Condition MET] Matrix A possesses the Consecutive Ones property in columns. "
                "Hence, A is Totally Unimodular (TUM)."
            )

    # 4. Exhaustive Determinant Check (for matrices with min(m, n) <= 6 to guarantee instant runtime)
    max_k = min(m, n)
    if not report["is_tum"] and max_k <= 6:
        all_tum = True
        bad_submatrix = None

        for k in range(1, max_k + 1):
            for row_comb in itertools.combinations(range(m), k):
                for col_comb in itertools.combinations(range(n), k):
                    sub = [[int_A[r][c] for c in col_comb] for r in row_comb]
                    det = compute_determinant(sub)
                    if det not in (-1, 0, 1):
                        all_tum = False
                        bad_submatrix = (row_comb, col_comb, det)
                        break
                if not all_tum:
                    break
            if not all_tum:
                break

        report["all_submatrices_checked"] = True
        if all_tum:
            report["is_tum"] = True
            report["sufficient_method"] = f"Exhaustive Submatrix Determinant Check (all square minors 1x1 to {max_k}x{max_k})"
            report["details"].append(
                f"[Exhaustive Check PASSED] Every square submatrix of A (sizes 1x1 to {max_k}x{max_k}) "
                "has determinant in {-1, 0, 1}. Hence, A is Totally Unimodular (TUM)."
            )
        else:
            r_idx, c_idx, det_val = bad_submatrix
            report["is_tum"] = False
            report["non_tum_reason"] = (
                f"Square submatrix with rows {[r+1 for r in r_idx]} and columns {[c+1 for c in c_idx]} "
                f"has determinant = {det_val} not in {{-1, 0, 1}}."
            )
            report["details"].append(f"[Determinant Check FAILED] {report['non_tum_reason']}")

    elif not report["is_tum"]:
        report["details"].append(
            f"Matrix size ({m}x{n}) exceeds fast exhaustive threshold; necessary conditions met but TUM not certified."
        )

    # 5. Ideal Polyhedron Verification (Hoffman-Kruskal Theorem)
    b_is_integer = all(val.denominator == 1 for val in b)
    if report["is_tum"] and b_is_integer:
        report["is_ideal"] = True
        report["details"].append(
            "[Hoffman-Kruskal Theorem] Matrix A is Totally Unimodular (TUM) AND vector b is integral. "
            "THEREFORE, the polyhedron P = {x >= 0 : Ax <= b} is an IDEAL POLYHEDRON (P = P_I). "
            "All extreme points (basic feasible solutions) are integers! "
            "The continuous LP relaxation directly solves the Integer Linear Program without any cutting planes."
        )
    elif report["is_tum"] and not b_is_integer:
        report["is_ideal"] = False
        report["details"].append(
            "[Polyhedron Status] Matrix A is TUM, but vector b is not strictly integral. "
            "Integrality of all extreme points is not guaranteed."
        )
    else:
        report["is_ideal"] = False
        report["details"].append(
            "[Polyhedron Status] Matrix A is NOT certified Totally Unimodular. "
            "The LP relaxation may produce fractional extreme points requiring Chvátal-Gomory cuts."
        )

    return report


# ---------------------------------------------------------------------
# PROBLEM FORMULATION DATA MODEL
# ---------------------------------------------------------------------

class LinearProgram:
    """Represents a Linear Programming problem for continuous and integer solving."""

    def __init__(self, obj_type: str, obj_coeffs: List[Fraction],
                 constraints: List[Tuple[List[Fraction], str, Fraction]],
                 var_names: Optional[List[str]] = None,
                 var_types: Optional[List[str]] = None):
        self.obj_type = obj_type.lower()
        if self.obj_type not in ("min", "max"):
            raise ValueError("Objective type must be 'min' or 'max'")
        self.obj_coeffs = [Fraction(c) for c in obj_coeffs]
        self.num_vars = len(self.obj_coeffs)
        self.constraints = []
        for coeffs, rel, rhs in constraints:
            if rel not in ("<=", ">=", "="):
                raise ValueError(f"Unsupported relation: '{rel}'")
            self.constraints.append(([Fraction(c) for c in coeffs], rel, Fraction(rhs)))

        if var_names:
            self.var_names = list(var_names)
        else:
            self.var_names = [f"x{j+1}" for j in range(self.num_vars)]

        if var_types:
            norm_types = []
            for vt in var_types:
                vtl = vt.lower()
                if vtl in ("int", "integer"):
                    norm_types.append("integer")
                elif vtl in ("bin", "binary"):
                    norm_types.append("binary")
                elif vtl in ("cont", "continuous", "real"):
                    norm_types.append("continuous")
                else:
                    raise ValueError(f"Unsupported variable type: '{vt}'")
            self.var_types = norm_types
        else:
            self.var_types = ["integer"] * self.num_vars

    def is_2d(self) -> bool:
        return self.num_vars == 2

    def integer_var_names(self) -> List[str]:
        return [
            name for name, vt in zip(self.var_names, self.var_types)
            if vt in ("integer", "binary")
        ]

    def effective_constraints(self) -> List[Tuple[List[Fraction], str, Fraction]]:
        """Returns original constraints plus implicit x_j <= 1 upper bounds for binary variables."""
        constrs = list(self.constraints)
        for j, vt in enumerate(self.var_types):
            if vt == "binary":
                row = [Fraction(0)] * self.num_vars
                row[j] = Fraction(1)
                constrs.append((row, "<=", Fraction(1)))
        return constrs


# ---------------------------------------------------------------------
# TABLEAU DATA STRUCTURE FOR SIMPLEX & DUAL SIMPLEX
# ---------------------------------------------------------------------

class SimplexTableau:
    """
    Standard Simplex Tableau with exact Fraction arithmetic.
    Row 0: -z equation for minimization: -z + sum(c_j * x_j) = RHS_0.
    Constraints: sum(a_ij * x_j) = b_i for i = 1..m.
    """

    def __init__(self, var_names: List[str], tableau: List[List[Fraction]],
                 basic_vars: List[str], obj_name: str = "-z",
                 orig_decision_vars: Optional[List[str]] = None):
        self.var_names = list(var_names)
        self.tableau = [[Fraction(v) for v in row] for row in tableau]
        self.basic_vars = list(basic_vars)
        self.obj_name = obj_name
        self.orig_decision_vars = list(orig_decision_vars) if orig_decision_vars else list(var_names)

    def num_rows(self) -> int:
        return len(self.tableau)

    def num_cols(self) -> int:
        return len(self.var_names)

    def copy(self) -> 'SimplexTableau':
        return SimplexTableau(
            self.var_names,
            self.tableau,
            self.basic_vars,
            self.obj_name,
            self.orig_decision_vars
        )

    def format_ascii(self, title: str = "") -> str:
        """Formats the tableau as a clean ASCII grid with RHS on the left."""
        col_headers = ["Basis", "RHS"] + self.var_names
        table_cells = [col_headers]

        # Row 0
        r0 = [self.obj_name, format_fraction(self.tableau[0][-1])]
        for j in range(self.num_cols()):
            r0.append(format_fraction(self.tableau[0][j]))
        table_cells.append(r0)

        # Rows 1..m
        for i in range(1, self.num_rows()):
            b_var = self.basic_vars[i - 1]
            r_cells = [b_var, format_fraction(self.tableau[i][-1])]
            for j in range(self.num_cols()):
                r_cells.append(format_fraction(self.tableau[i][j]))
            table_cells.append(r_cells)

        col_widths = []
        for c in range(len(col_headers)):
            w = max(len(table_cells[r][c]) for r in range(len(table_cells)))
            col_widths.append(max(w, 5))

        sep = "+" + "+".join("-" * (w + 2) for w in col_widths) + "+"
        lines = []
        if title:
            lines.append(f"--- {title} ---")
        lines.append(sep)
        # Header
        h_str = "| " + " | ".join(f"{table_cells[0][c]:<{col_widths[c]}}" for c in range(len(col_headers))) + " |"
        lines.append(h_str)
        lines.append(sep)
        # Row 0
        r0_str = "| " + " | ".join(f"{table_cells[1][c]:>{col_widths[c]}}" if c > 0 else f"{table_cells[1][c]:<{col_widths[c]}}" for c in range(len(col_headers))) + " |"
        lines.append(r0_str)
        lines.append(sep)
        # Constraint rows
        for r in range(2, len(table_cells)):
            row_str = "| " + " | ".join(f"{table_cells[r][c]:>{col_widths[c]}}" if c > 0 else f"{table_cells[r][c]:<{col_widths[c]}}" for c in range(len(col_headers))) + " |"
            lines.append(row_str)
        lines.append(sep)
        return "\n".join(lines)


# ---------------------------------------------------------------------
# PRIMAL & DUAL SIMPLEX SOLVER ENGINES
# ---------------------------------------------------------------------

def solve_primal_simplex(tab: SimplexTableau, max_iterations: int = 500) -> Tuple[str, List[str]]:
    """
    Solves continuous LP relaxation using Primal Simplex with Bland's anti-cycling rule.
    Returns (status, logs).
    """
    logs = []
    iteration = 0
    var_order = {name: idx for idx, name in enumerate(tab.var_names)}

    while iteration < max_iterations:
        # Check optimality: all reduced costs in Row 0 >= 0
        candidates = []
        for j in range(tab.num_cols()):
            if tab.tableau[0][j] < 0:
                candidates.append(j)

        if not candidates:
            logs.append(f"[Primal Simplex] Optimality reached at Iteration {iteration}.")
            return "OPTIMAL", logs

        # Bland's smallest-subscript rule for entering variable
        pivot_col = min(candidates, key=lambda j: var_order.get(tab.var_names[j], j))
        entering_var = tab.var_names[pivot_col]

        # Minimum Ratio Test (theta) with Bland's smallest-index tie breaking
        tied_rows: List[Tuple[Fraction, int, int]] = []
        min_theta = None

        for i in range(1, tab.num_rows()):
            coeff = tab.tableau[i][pivot_col]
            if coeff > 0:
                theta = tab.tableau[i][-1] / coeff
                b_idx = var_order.get(tab.basic_vars[i - 1], 10**9 + i)
                if min_theta is None or theta < min_theta:
                    min_theta = theta
                    tied_rows = [(theta, b_idx, i)]
                elif theta == min_theta:
                    tied_rows.append((theta, b_idx, i))

        if not tied_rows:
            logs.append(f"[Primal Simplex] Unbounded problem detected for entering variable '{entering_var}'.")
            return "UNBOUNDED", logs

        tied_rows.sort(key=lambda item: item[1])
        leaving_row = tied_rows[0][2]
        pivot_elem = tab.tableau[leaving_row][pivot_col]

        # Elementary Row Operations
        for c in range(len(tab.tableau[leaving_row])):
            tab.tableau[leaving_row][c] /= pivot_elem

        for i in range(tab.num_rows()):
            if i != leaving_row:
                factor = tab.tableau[i][pivot_col]
                if factor != 0:
                    for c in range(len(tab.tableau[i])):
                        tab.tableau[i][c] -= factor * tab.tableau[leaving_row][c]

        tab.basic_vars[leaving_row - 1] = entering_var
        iteration += 1

    return "ITERATION_LIMIT", logs


def solve_dual_simplex(tab: SimplexTableau, max_iterations: int = 200) -> Tuple[str, List[str]]:
    """
    Executes Dual Simplex steps until all RHS values >= 0 (primal feasibility restored).
    Used after injecting a Chvatal-Gomory / GMI cut row.
    """
    logs = []
    iteration = 0
    var_order = {name: idx for idx, name in enumerate(tab.var_names)}

    while iteration < max_iterations:
        # Check primal feasibility: all RHS_i >= 0 for i >= 1
        negative_rows = [i for i in range(1, tab.num_rows()) if tab.tableau[i][-1] < 0]

        if not negative_rows:
            logs.append("[Dual Simplex] Primal feasibility restored.")
            return "OPTIMAL", logs

        # Most negative RHS for leaving row, tie-broken by Bland index
        most_neg = min(tab.tableau[r][-1] for r in negative_rows)
        tied_leaving = [r for r in negative_rows if tab.tableau[r][-1] == most_neg]
        leaving_row = min(tied_leaving, key=lambda r: var_order.get(tab.basic_vars[r - 1], 10**9 + r))
        leaving_var = tab.basic_vars[leaving_row - 1]

        # Dual Minimum Ratio Test across negative coefficients in the leaving row: min(c_j / |a_rj|)
        ratios = []
        for j in range(tab.num_cols()):
            coeff = tab.tableau[leaving_row][j]
            if coeff < 0:
                cost = tab.tableau[0][j]
                ratio = cost / abs(coeff)
                ratios.append((ratio, j))

        if not ratios:
            logs.append(f"[Dual Simplex] Infeasibility detected: leaving row '{leaving_var}' has no negative coefficients.")
            return "INFEASIBLE", logs

        # Bland's smallest-index tie-breaking on min ratio
        min_ratio = min(r[0] for r in ratios)
        candidates = [j for ratio, j in ratios if ratio == min_ratio]
        pivot_col = min(candidates, key=lambda j: var_order.get(tab.var_names[j], j))
        entering_var = tab.var_names[pivot_col]

        pivot_elem = tab.tableau[leaving_row][pivot_col]

        for c in range(len(tab.tableau[leaving_row])):
            tab.tableau[leaving_row][c] /= pivot_elem

        for i in range(tab.num_rows()):
            if i != leaving_row:
                factor = tab.tableau[i][pivot_col]
                if factor != 0:
                    for c in range(len(tab.tableau[i])):
                        tab.tableau[i][c] -= factor * tab.tableau[leaving_row][c]

        tab.basic_vars[leaving_row - 1] = entering_var
        iteration += 1

    return "ITERATION_LIMIT", logs


def _alloc_unique_var(prefix: str, counter_start: int, used: Set[str]) -> Tuple[str, int]:
    idx = counter_start
    while f"{prefix}{idx}" in used:
        idx += 1
    name = f"{prefix}{idx}"
    used.add(name)
    return name, idx + 1


def solve_lp_relaxation_twophase(
    lp: LinearProgram,
    extra_constraints: Optional[List[Tuple[List[Fraction], str, Fraction]]] = None
) -> Dict[str, Any]:
    """
    Constructs and solves the continuous LP relaxation of `lp` plus any `extra_constraints`
    (such as branching bounds or inherited cuts) using exact Two-Phase Simplex with Bland's rule.

    Returns a dictionary containing:
    - status: 'OPTIMAL', 'INFEASIBLE', 'UNBOUNDED', or 'ITERATION_LIMIT'
    - tableau: solved SimplexTableau ready for inspection or Dual Simplex cut injection (if OPTIMAL)
    - initial_tableau_ascii: ASCII rendering before Simplex pivots
    - solution: Dict[var_name, Fraction] (if OPTIMAL)
    - opt_val: Fraction (if OPTIMAL)
    - var_definitions: affine definitions of all variables in terms of original decision variables
    - var_integrality: Dict[var_name, bool] indicating whether each variable must be integer
    - used_two_phase: bool
    """
    all_constraints = lp.effective_constraints()
    if extra_constraints:
        all_constraints.extend(extra_constraints)

    norm_constraints: List[Tuple[List[Fraction], str, Fraction]] = []
    for coeffs, rel, rhs in all_constraints:
        r_coeffs = [Fraction(c) for c in coeffs]
        r_rel = rel
        r_rhs = Fraction(rhs)
        if r_rhs < 0:
            r_coeffs = [-c for c in r_coeffs]
            r_rhs = -r_rhs
            if r_rel == "<=":
                r_rel = ">="
            elif r_rel == ">=":
                r_rel = "<="
        norm_constraints.append((r_coeffs, r_rel, r_rhs))

    needs_phase1 = any(rel in (">=", "=") for _, rel, _ in norm_constraints)

    var_definitions: Dict[str, Tuple[Dict[str, Fraction], Fraction]] = {}
    var_integrality: Dict[str, bool] = {}
    used_names: Set[str] = set(lp.var_names)

    for j, var in enumerate(lp.var_names):
        var_definitions[var] = ({var: Fraction(1)}, Fraction(0))
        var_integrality[var] = (lp.var_types[j] in ("integer", "binary"))

    slacks: List[str] = []
    surpluses: List[str] = []
    artificials: List[str] = []
    row_aux: List[Dict[str, Optional[str]]] = []

    s_counter = 1
    e_counter = 1
    y_counter = 1

    for r_coeffs, r_rel, r_rhs in norm_constraints:
        aux: Dict[str, Optional[str]] = {"slack": None, "surplus": None, "artificial": None}
        all_int_terms = (r_rhs.denominator == 1) and all(
            c.denominator == 1 and var_integrality[lp.var_names[j]]
            for j, c in enumerate(r_coeffs) if c != 0
        )
        if r_rel == "<=":
            s_name, s_counter = _alloc_unique_var("s", s_counter, used_names)
            slacks.append(s_name)
            aux["slack"] = s_name
            gammas = {lp.var_names[j]: -r_coeffs[j] for j in range(lp.num_vars)}
            var_definitions[s_name] = (gammas, r_rhs)
            var_integrality[s_name] = all_int_terms
        elif r_rel == ">=":
            e_name, e_counter = _alloc_unique_var("e", e_counter, used_names)
            y_name, y_counter = _alloc_unique_var("y", y_counter, used_names)
            surpluses.append(e_name)
            artificials.append(y_name)
            aux["surplus"] = e_name
            aux["artificial"] = y_name
            gammas = {lp.var_names[j]: r_coeffs[j] for j in range(lp.num_vars)}
            var_definitions[e_name] = (gammas, -r_rhs)
            var_integrality[e_name] = all_int_terms
        elif r_rel == "=":
            y_name, y_counter = _alloc_unique_var("y", y_counter, used_names)
            artificials.append(y_name)
            aux["artificial"] = y_name
        row_aux.append(aux)

    full_var_names = list(lp.var_names) + slacks + surpluses + artificials
    n_orig = lp.num_vars
    n_full = len(full_var_names)
    m = len(norm_constraints)

    matrix_A: List[List[Fraction]] = []
    basic_vars: List[str] = []
    for i, (r_coeffs, r_rel, r_rhs) in enumerate(norm_constraints):
        row = [Fraction(0)] * (n_full + 1)
        for j in range(n_orig):
            row[j] = r_coeffs[j]
        aux = row_aux[i]
        if aux["slack"]:
            row[full_var_names.index(aux["slack"])] = Fraction(1)
            basic_vars.append(aux["slack"])
        if aux["surplus"]:
            row[full_var_names.index(aux["surplus"])] = Fraction(-1)
        if aux["artificial"]:
            row[full_var_names.index(aux["artificial"])] = Fraction(1)
            basic_vars.append(aux["artificial"])
        row[-1] = r_rhs
        matrix_A.append(row)

    # Original objective row (-z + sum c'_j x_j = 0)
    orig_row0 = [Fraction(0)] * (n_full + 1)
    for j in range(n_orig):
        c = lp.obj_coeffs[j]
        orig_row0[j] = -c if lp.obj_type == "max" else c

    if not needs_phase1:
        tab = SimplexTableau(
            full_var_names,
            [list(orig_row0)] + [list(r) for r in matrix_A],
            basic_vars,
            obj_name="-z",
            orig_decision_vars=lp.var_names
        )
        init_ascii = tab.copy().format_ascii("Initial Tableau (LP Relaxation)")
        status, logs = solve_primal_simplex(tab)
    else:
        # Phase 1: minimize W = sum(artificials)
        phase1_row0 = [Fraction(0)] * (n_full + 1)
        for a_name in artificials:
            phase1_row0[full_var_names.index(a_name)] = Fraction(1)
        for i, aux in enumerate(row_aux):
            if aux["artificial"]:
                for col in range(n_full + 1):
                    phase1_row0[col] -= matrix_A[i][col]

        p1_tab = SimplexTableau(
            full_var_names,
            [phase1_row0] + [list(r) for r in matrix_A],
            basic_vars,
            obj_name="-w",
            orig_decision_vars=lp.var_names
        )
        init_ascii = p1_tab.copy().format_ascii("Initial Phase 1 Tableau (LP Relaxation)")
        p1_status, logs = solve_primal_simplex(p1_tab)

        if p1_status != "OPTIMAL":
            return {
                "status": p1_status,
                "tableau": p1_tab,
                "initial_tableau_ascii": init_ascii,
                "solution": None,
                "opt_val": None,
                "var_definitions": var_definitions,
                "var_integrality": var_integrality,
                "used_two_phase": True,
                "logs": logs
            }

        if p1_tab.tableau[0][-1] != 0:
            return {
                "status": "INFEASIBLE",
                "tableau": p1_tab,
                "initial_tableau_ascii": init_ascii,
                "solution": None,
                "opt_val": None,
                "var_definitions": var_definitions,
                "var_integrality": var_integrality,
                "used_two_phase": True,
                "logs": logs
            }

        # Drive any degenerate artificial variables out of the basis
        art_set = set(artificials)
        redundant_rows: Set[int] = set()
        for r_idx in range(1, p1_tab.num_rows()):
            if p1_tab.basic_vars[r_idx - 1] in art_set:
                pivoted = False
                for c_idx, v_name in enumerate(p1_tab.var_names):
                    if v_name not in art_set and p1_tab.tableau[r_idx][c_idx] != 0:
                        p_val = p1_tab.tableau[r_idx][c_idx]
                        p1_tab.tableau[r_idx] = [val / p_val for val in p1_tab.tableau[r_idx]]
                        for r_other in range(p1_tab.num_rows()):
                            if r_other != r_idx and p1_tab.tableau[r_other][c_idx] != 0:
                                factor = p1_tab.tableau[r_other][c_idx]
                                p1_tab.tableau[r_other] = [
                                    p1_tab.tableau[r_other][c] - factor * p1_tab.tableau[r_idx][c]
                                    for c in range(len(p1_tab.tableau[r_other]))
                                ]
                        p1_tab.basic_vars[r_idx - 1] = v_name
                        pivoted = True
                        break
                if not pivoted:
                    redundant_rows.add(r_idx)

        # Transition to Phase 2: strip artificial columns and redundant zero rows
        keep_cols = [c for c, name in enumerate(p1_tab.var_names) if name not in art_set]
        p2_var_names = [p1_tab.var_names[c] for c in keep_cols]
        p2_rows: List[List[Fraction]] = []
        p2_basic: List[str] = []
        for r_idx in range(1, p1_tab.num_rows()):
            if r_idx in redundant_rows:
                continue
            p2_rows.append([p1_tab.tableau[r_idx][c] for c in keep_cols] + [p1_tab.tableau[r_idx][-1]])
            p2_basic.append(p1_tab.basic_vars[r_idx - 1])

        p2_row0 = [Fraction(0)] * (len(p2_var_names) + 1)
        for j, var_name in enumerate(p2_var_names):
            orig_idx = full_var_names.index(var_name)
            p2_row0[j] = orig_row0[orig_idx]
        p2_row0[-1] = orig_row0[-1]

        for r_idx, b_var in enumerate(p2_basic):
            c_idx = p2_var_names.index(b_var)
            factor = p2_row0[c_idx]
            if factor != 0:
                for col in range(len(p2_row0)):
                    p2_row0[col] -= factor * p2_rows[r_idx][col]

        tab = SimplexTableau(
            p2_var_names,
            [p2_row0] + p2_rows,
            p2_basic,
            obj_name="-z",
            orig_decision_vars=lp.var_names
        )
        status, p2_logs = solve_primal_simplex(tab)
        logs.extend(p2_logs)

    if status != "OPTIMAL":
        return {
            "status": status,
            "tableau": tab,
            "initial_tableau_ascii": init_ascii,
            "solution": None,
            "opt_val": None,
            "var_definitions": var_definitions,
            "var_integrality": var_integrality,
            "used_two_phase": needs_phase1,
            "logs": logs
        }

    solution: Dict[str, Fraction] = {}
    for name in lp.var_names:
        if name in tab.basic_vars:
            r_idx = tab.basic_vars.index(name) + 1
            solution[name] = tab.tableau[r_idx][-1]
        else:
            solution[name] = Fraction(0)

    raw_rhs0 = tab.tableau[0][-1]
    opt_val = -raw_rhs0 if lp.obj_type == "min" else raw_rhs0

    return {
        "status": "OPTIMAL",
        "tableau": tab,
        "initial_tableau_ascii": init_ascii,
        "solution": solution,
        "opt_val": opt_val,
        "var_definitions": var_definitions,
        "var_integrality": var_integrality,
        "used_two_phase": needs_phase1,
        "logs": logs
    }


# ---------------------------------------------------------------------
# CUTTING PLANE ALGORITHM & CHVÁTAL-GOMORY GENERATOR
# ---------------------------------------------------------------------

class CutDetail:
    """Detailed mathematical trace of an injected Chvátal-Gomory or GMI cut."""

    def __init__(self, cut_index: int, generating_var: str, generating_row: int,
                 rhs_fractional: Fraction, floor_rhs: int, frac_rhs: Fraction,
                 slack_coefficients: Dict[str, Fraction],
                 tableau_inequality_str: str,
                 orig_vars_inequality_str: str,
                 orig_vars_coefficients: Dict[str, Fraction],
                 orig_vars_rhs: Fraction,
                 tableau_snapshot_before: str,
                 tableau_snapshot_after: str,
                 cut_kind: str = "CG",
                 node_id: int = 0,
                 opt_val_before: Optional[Fraction] = None,
                 opt_val_after: Optional[Fraction] = None):
        self.cut_index = cut_index
        self.generating_var = generating_var
        self.generating_row = generating_row
        self.rhs_fractional = rhs_fractional
        self.floor_rhs = floor_rhs
        self.frac_rhs = frac_rhs
        self.slack_coefficients = dict(slack_coefficients)
        self.tableau_inequality_str = tableau_inequality_str
        self.orig_vars_inequality_str = orig_vars_inequality_str
        self.orig_vars_coefficients = dict(orig_vars_coefficients)
        self.orig_vars_rhs = orig_vars_rhs
        self.tableau_snapshot_before = tableau_snapshot_before
        self.tableau_snapshot_after = tableau_snapshot_after
        self.cut_kind = cut_kind
        self.node_id = node_id
        self.opt_val_before = opt_val_before
        self.opt_val_after = opt_val_after


def reconstruct_cut_in_original_vars(
    lp: LinearProgram,
    f_coeffs: Dict[str, Fraction],
    f0: Fraction,
    var_definitions: Dict[str, Tuple[Dict[str, Fraction], Fraction]]
) -> Tuple[Dict[str, Fraction], Fraction, str]:
    """
    Substitutes the affine definitions v = const_v + sum(gamma_{v,j} * x_j) of all variables v
    into the tableau cut sum(f_v * v) >= f0, obtaining the cut directly in decision variables:
    sum(alpha_j * x_j) <= beta
    """
    alpha = {var: Fraction(0) for var in lp.var_names}
    sum_f_const = Fraction(0)

    for v, f_val in f_coeffs.items():
        if f_val == 0:
            continue
        if v in var_definitions:
            gammas, const_v = var_definitions[v]
            sum_f_const += f_val * const_v
            for var in lp.var_names:
                if var in gammas:
                    alpha[var] -= f_val * gammas[var]
        else:
            if v in lp.var_names:
                alpha[v] -= f_val

    beta = sum_f_const - f0

    # Format algebraic string
    terms = []
    for var in lp.var_names:
        c = alpha[var]
        if c != 0:
            if c == 1:
                terms.append(var)
            elif c == -1:
                terms.append(f"-{var}")
            else:
                terms.append(f"{format_fraction(c)}*{var}")
    lhs_str = " + ".join(terms) if terms else "0"
    lhs_str = lhs_str.replace("+ -", "- ")
    ineq_str = f"{lhs_str} <= {format_fraction(beta)}"

    return alpha, beta, ineq_str


def _generate_and_apply_cut_round(
    lp: LinearProgram,
    tab: SimplexTableau,
    var_definitions: Dict[str, Tuple[Dict[str, Fraction], Fraction]],
    var_integrality: Dict[str, bool],
    cut_index: int,
    cut_type: str = "auto",
    node_id: int = 0
) -> Optional[Tuple[str, CutDetail]]:
    """
    Finds a fractional integer-constrained basic variable in `tab`, derives either a
    Chvátal-Gomory (CG) or Gomory Mixed-Integer (GMI) cut, injects it, and runs Dual Simplex.
    Returns None if all integer-constrained basic variables are integral; otherwise (dual_status, CutDetail).
    """
    int_var_set = set(lp.integer_var_names())
    fractional_vars = []
    for i in range(1, tab.num_rows()):
        b_var = tab.basic_vars[i - 1]
        if b_var in int_var_set:
            rhs_val = tab.tableau[i][-1]
            f0 = frac_part(rhs_val)
            if f0 != 0:
                fractional_vars.append((i, b_var, rhs_val, f0))

    if not fractional_vars:
        return None

    # Prefer row with largest f0 (matching canonical CG order) when cut_type is cg/auto for pure ILP,
    # and deterministic tie-breaker on row index.
    def frac_score(item):
        return (item[3], -item[0])

    chosen_row, chosen_var, chosen_rhs, f0 = max(fractional_vars, key=frac_score)

    raw_rhs0_before = tab.tableau[0][-1]
    opt_before = -raw_rhs0_before if lp.obj_type == "min" else raw_rhs0_before

    # Decide CG vs GMI: pure CG is valid iff all nonzero non-basic coefficients in chosen_row are integer-guaranteed vars
    basic_set = set(tab.basic_vars)
    can_use_cg = True
    for j in range(tab.num_cols()):
        v_name = tab.var_names[j]
        if v_name in basic_set:
            continue
        if tab.tableau[chosen_row][j] != 0 and not var_integrality.get(v_name, False):
            can_use_cg = False
            break

    use_kind = "CG" if (cut_type == "cg" or (cut_type == "auto" and can_use_cg)) else "GMI"

    cut_name = f"s_cut{cut_index}"
    while cut_name in var_definitions or cut_name in tab.var_names:
        cut_name = f"{cut_name}_n{node_id}"
    tab_before_cut = tab.copy().format_ascii(f"Tableau before Cut {cut_index}")

    f_coeffs: Dict[str, Fraction] = {}
    cut_row: List[Fraction] = []
    ineq_terms: List[str] = []

    for j in range(tab.num_cols()):
        var = tab.var_names[j]
        coeff = tab.tableau[chosen_row][j]
        if var in basic_set:
            gj = Fraction(0)
        elif use_kind == "CG":
            gj = frac_part(coeff)
        else:
            # Gomory Mixed-Integer (GMI) cut formula
            if var_integrality.get(var, False):
                fj = frac_part(coeff)
                if fj <= f0:
                    gj = fj
                else:
                    gj = (f0 * (Fraction(1) - fj)) / (Fraction(1) - f0)
            else:
                if coeff >= 0:
                    gj = coeff
                else:
                    gj = (f0 * abs(coeff)) / (Fraction(1) - f0)

        f_coeffs[var] = gj
        cut_row.append(-gj)
        if gj != 0:
            gj_str = format_fraction(gj)
            ineq_terms.append(f"{gj_str}*{var}" if gj_str != "1" else var)

    cut_row.append(Fraction(1))
    cut_row.append(-f0)

    for i in range(tab.num_rows()):
        tab.tableau[i].insert(-1, Fraction(0))

    tab.var_names.append(cut_name)
    tab.tableau.append(cut_row)
    tab.basic_vars.append(cut_name)

    tableau_ineq_str = (" + ".join(ineq_terms) if ineq_terms else "0") + f" >= {format_fraction(f0)}"

    orig_coeffs, orig_rhs, orig_ineq_str = reconstruct_cut_in_original_vars(lp, f_coeffs, f0, var_definitions)
    cut_gammas = {var: -orig_coeffs[var] for var in lp.var_names}
    var_definitions[cut_name] = (cut_gammas, orig_rhs)
    var_integrality[cut_name] = (
        use_kind == "CG"
        and orig_rhs.denominator == 1
        and all(
            orig_coeffs[v].denominator == 1 and var_integrality.get(v, False)
            for v in lp.var_names if orig_coeffs[v] != 0
        )
    )

    dual_status, _ = solve_dual_simplex(tab)
    tab_after_cut = tab.copy().format_ascii(f"Tableau after Cut {cut_index} (Dual Simplex)")

    raw_rhs0_after = tab.tableau[0][-1]
    opt_after = -raw_rhs0_after if lp.obj_type == "min" else raw_rhs0_after

    cut_detail = CutDetail(
        cut_index=cut_index,
        generating_var=chosen_var,
        generating_row=chosen_row,
        rhs_fractional=chosen_rhs,
        floor_rhs=floor_frac(chosen_rhs),
        frac_rhs=f0,
        slack_coefficients=f_coeffs,
        tableau_inequality_str=tableau_ineq_str,
        orig_vars_inequality_str=orig_ineq_str,
        orig_vars_coefficients=orig_coeffs,
        orig_vars_rhs=orig_rhs,
        tableau_snapshot_before=tab_before_cut,
        tableau_snapshot_after=tab_after_cut,
        cut_kind=use_kind,
        node_id=node_id,
        opt_val_before=opt_before,
        opt_val_after=opt_after
    )
    return dual_status, cut_detail


def solve_chvatal_gomory(
    lp: LinearProgram,
    max_cuts: int = 15,
    cut_type: str = "auto"
) -> Dict[str, Any]:
    """
    Solves an ILP/MILP using pure cutting planes (Chvatal-Gomory / GMI Cuts) at the root relaxation.
    Supports <=, >=, and = constraints via Two-Phase Simplex + Dual Simplex re-optimization.
    """
    lp_res = solve_lp_relaxation_twophase(lp)
    current_status = lp_res["status"]
    tab: SimplexTableau = lp_res["tableau"]
    var_definitions = lp_res["var_definitions"]
    var_integrality = lp_res["var_integrality"]

    relax_tab_before = lp_res["initial_tableau_ascii"]
    relax_tab_after = tab.copy().format_ascii("Optimal LP Relaxation Tableau")
    relax_solution = lp_res["solution"] if lp_res["solution"] is not None else {v: Fraction(0) for v in lp.var_names}
    relax_opt_val = lp_res["opt_val"]

    cuts_history: List[CutDetail] = []
    cut_count = 0

    while current_status == "OPTIMAL" and cut_count < max_cuts:
        step_res = _generate_and_apply_cut_round(
            lp, tab, var_definitions, var_integrality,
            cut_index=cut_count + 1, cut_type=cut_type, node_id=0
        )
        if step_res is None:
            break
        cut_count += 1
        current_status, cut_detail = step_res
        cuts_history.append(cut_detail)
        if current_status != "OPTIMAL":
            break

    int_var_set = set(lp.integer_var_names())
    final_solution = {}
    is_integer = (current_status == "OPTIMAL")
    if current_status == "OPTIMAL":
        for name in lp.var_names:
            if name in tab.basic_vars:
                r_idx = tab.basic_vars.index(name) + 1
                val = tab.tableau[r_idx][-1]
            else:
                val = Fraction(0)
            final_solution[name] = val
            if name in int_var_set and val.denominator != 1:
                is_integer = False
        raw_rhs0_final = tab.tableau[0][-1]
        final_opt_val = -raw_rhs0_final if lp.obj_type == "min" else raw_rhs0_final
    else:
        final_solution = {name: Fraction(0) for name in lp.var_names}
        final_opt_val = None

    return {
        "method": "cutting-plane",
        "status": current_status,
        "is_integer": is_integer,
        "relax_solution": relax_solution,
        "relax_opt_val": relax_opt_val,
        "relax_tab_before": relax_tab_before,
        "relax_tab_after": relax_tab_after,
        "cuts_history": cuts_history,
        "total_cuts": len(cuts_history),
        "final_solution": final_solution,
        "final_opt_val": final_opt_val,
        "final_tableau_ascii": tab.format_ascii("Final Tableau")
    }


# ---------------------------------------------------------------------
# BRANCH AND BOUND & BRANCH AND CUT ENGINES
# ---------------------------------------------------------------------

class BBNode:
    """Represents a single node in a Branch-and-Bound or Branch-and-Cut search tree."""

    def __init__(
        self,
        node_id: int,
        parent_id: Optional[int],
        depth: int,
        branch_desc: str,
        var_bounds: Dict[str, Tuple[int, Optional[int]]],
        inherited_cuts: List[CutDetail],
        branch_line: Optional[Tuple[str, str, int]] = None
    ):
        self.node_id = node_id
        self.parent_id = parent_id
        self.depth = depth
        self.branch_desc = branch_desc
        self.var_bounds = dict(var_bounds)
        self.inherited_cuts = list(inherited_cuts)
        self.branch_line = branch_line  # (var_name, "<=" | ">=", int_bound) if created by branching

        self.lp_status: str = "PENDING"
        self.lp_val_before_cuts: Optional[Fraction] = None
        self.lp_sol_before_cuts: Optional[Dict[str, Fraction]] = None
        self.cuts_generated: List[CutDetail] = []
        self.active_cuts: List[CutDetail] = list(inherited_cuts)
        self.lp_val: Optional[Fraction] = None
        self.lp_sol: Optional[Dict[str, Fraction]] = None
        self.outcome: str = "PENDING"
        self.outcome_reason: str = ""
        self.branch_var: Optional[str] = None
        self.branch_val: Optional[Fraction] = None
        self.children_ids: List[int] = []
        self.incumbent_at_step: Optional[Fraction] = None
        self.incumbent_sol_at_step: Optional[Dict[str, Fraction]] = None

    def bounds_summary(self, var_names: List[str]) -> str:
        parts = []
        for v in var_names:
            lb, ub = self.var_bounds[v]
            if ub is not None:
                parts.append(f"{lb} <= {v} <= {ub}")
            elif lb > 0:
                parts.append(f"{v} >= {lb}")
        return ", ".join(parts) if parts else "Root domain (x >= 0)"

    def to_extra_constraints(self, lp: LinearProgram) -> List[Tuple[List[Fraction], str, Fraction]]:
        """Builds explicit branching bound constraints (excluding default x >= 0 and binary <= 1 handled by lp) plus inherited cuts."""
        extra: List[Tuple[List[Fraction], str, Fraction]] = []
        for j, name in enumerate(lp.var_names):
            lb, ub = self.var_bounds[name]
            default_lb = 0
            default_ub = 1 if lp.var_types[j] == "binary" else None
            if lb > default_lb:
                coeffs = [Fraction(0)] * lp.num_vars
                coeffs[j] = Fraction(1)
                extra.append((coeffs, ">=", Fraction(lb)))
            if ub is not None and (default_ub is None or ub < default_ub):
                coeffs = [Fraction(0)] * lp.num_vars
                coeffs[j] = Fraction(1)
                extra.append((coeffs, "<=", Fraction(ub)))
        for cut in self.inherited_cuts:
            coeffs = [cut.orig_vars_coefficients[v] for v in lp.var_names]
            extra.append((coeffs, "<=", cut.orig_vars_rhs))
        return extra


def _is_better_obj(candidate: Fraction, incumbent: Optional[Fraction], obj_type: str) -> bool:
    if incumbent is None:
        return True
    return (candidate > incumbent) if obj_type == "max" else (candidate < incumbent)


def _is_bound_dominated(bound_val: Fraction, incumbent: Optional[Fraction], obj_type: str) -> bool:
    if incumbent is None:
        return False
    return (bound_val <= incumbent) if obj_type == "max" else (bound_val >= incumbent)


def _extract_sol_from_tableau(lp: LinearProgram, tab: SimplexTableau) -> Tuple[Dict[str, Fraction], Fraction]:
    sol: Dict[str, Fraction] = {}
    for name in lp.var_names:
        if name in tab.basic_vars:
            r_idx = tab.basic_vars.index(name) + 1
            sol[name] = tab.tableau[r_idx][-1]
        else:
            sol[name] = Fraction(0)
    raw_rhs0 = tab.tableau[0][-1]
    opt_val = -raw_rhs0 if lp.obj_type == "min" else raw_rhs0
    return sol, opt_val


def _check_integer_feasibility(lp: LinearProgram, sol: Dict[str, Fraction]) -> List[Tuple[int, str, Fraction]]:
    """Returns [(var_idx, var_name, val)] for all fractional integer/binary decision variables."""
    int_set = set(lp.integer_var_names())
    frac_list: List[Tuple[int, str, Fraction]] = []
    for j, name in enumerate(lp.var_names):
        if name in int_set and sol[name].denominator != 1:
            frac_list.append((j, name, sol[name]))
    return frac_list


def _select_branch_var(frac_vars: List[Tuple[int, str, Fraction]], branch_rule: str) -> Tuple[int, str, Fraction]:
    if branch_rule == "first-fractional":
        return frac_vars[0]
    # most-fractional: maximize distance from nearest integer = min(f, 1 - f), tie-break smallest var index
    def score(item):
        f = frac_part(item[2])
        dist = min(f, Fraction(1) - f)
        return (dist, -item[0])
    return max(frac_vars, key=score)


def _solve_tree_search(
    lp: LinearProgram,
    enable_cuts: bool = False,
    strategy: str = "best-bound",
    branch_rule: str = "most-fractional",
    cut_type: str = "auto",
    cuts_per_node: int = 2,
    cut_depth: int = 5,
    max_nodes: int = 200
) -> Dict[str, Any]:
    """
    Unified exact Branch-and-Bound / Branch-and-Cut tree search engine.
    - Every created node is solved once with Two-Phase Simplex (and Dual Simplex cuts when enable_cuts=True).
    - Open nodes store their exact node LP bound for Best-Bound, DFS, or BFS selection.
    """
    init_bounds: Dict[str, Tuple[int, Optional[int]]] = {
        name: (0, 1 if lp.var_types[j] == "binary" else None)
        for j, name in enumerate(lp.var_names)
    }

    nodes: List[BBNode] = []
    global_cuts_history: List[CutDetail] = []
    incumbent_val: Optional[Fraction] = None
    incumbent_sol: Optional[Dict[str, Fraction]] = None
    incumbent_node_id: Optional[int] = None
    node_counter = 0
    global_cut_counter = 0

    relax_solution: Dict[str, Fraction] = {v: Fraction(0) for v in lp.var_names}
    relax_opt_val: Optional[Fraction] = None
    relax_tab_before = ""
    relax_tab_after = ""
    root_status = "OPTIMAL"

    # Open queue items: list of BBNode objects awaiting branching
    open_queue: List[BBNode] = []

    def evaluate_node(node: BBNode) -> None:
        nonlocal incumbent_val, incumbent_sol, incumbent_node_id, global_cut_counter
        nonlocal relax_solution, relax_opt_val, relax_tab_before, relax_tab_after, root_status

        extra = node.to_extra_constraints(lp)
        lp_res = solve_lp_relaxation_twophase(lp, extra_constraints=extra)
        node.lp_status = lp_res["status"]

        if node.node_id == 0:
            root_status = lp_res["status"]
            relax_tab_before = lp_res["initial_tableau_ascii"]
            relax_tab_after = lp_res["tableau"].copy().format_ascii("Root LP Relaxation Tableau")
            if lp_res["solution"] is not None:
                relax_solution = dict(lp_res["solution"])
            relax_opt_val = lp_res["opt_val"]

        if lp_res["status"] == "INFEASIBLE":
            node.outcome = "INFEASIBLE"
            node.outcome_reason = "Node LP relaxation is infeasible under branching bounds."
            node.incumbent_at_step = incumbent_val
            node.incumbent_sol_at_step = dict(incumbent_sol) if incumbent_sol else None
            return

        if lp_res["status"] == "UNBOUNDED":
            node.outcome = "UNBOUNDED"
            node.outcome_reason = "Node LP relaxation is unbounded."
            node.incumbent_at_step = incumbent_val
            node.incumbent_sol_at_step = dict(incumbent_sol) if incumbent_sol else None
            return

        if lp_res["status"] != "OPTIMAL":
            node.outcome = lp_res["status"]
            node.outcome_reason = f"Node LP solver terminated with status {lp_res['status']}."
            node.incumbent_at_step = incumbent_val
            node.incumbent_sol_at_step = dict(incumbent_sol) if incumbent_sol else None
            return

        tab: SimplexTableau = lp_res["tableau"]
        var_definitions = lp_res["var_definitions"]
        var_integrality = lp_res["var_integrality"]

        cur_sol, cur_val = _extract_sol_from_tableau(lp, tab)
        node.lp_val_before_cuts = cur_val
        node.lp_sol_before_cuts = dict(cur_sol)
        node.lp_val = cur_val
        node.lp_sol = dict(cur_sol)

        # Bound pruning BEFORE spending time generating cuts
        if _is_bound_dominated(cur_val, incumbent_val, lp.obj_type):
            node.outcome = "PRUNED_BOUND"
            comp = "<=" if lp.obj_type == "max" else ">="
            node.outcome_reason = (
                f"LP bound {format_fraction(cur_val)} {comp} Incumbent {format_fraction(incumbent_val)}."
            )
            node.incumbent_at_step = incumbent_val
            node.incumbent_sol_at_step = dict(incumbent_sol) if incumbent_sol else None
            return

        # Optional cutting-plane strengthening in Branch and Cut mode
        if enable_cuts and node.depth <= cut_depth and cuts_per_node > 0:
            node_cut_rounds = 0
            while node_cut_rounds < cuts_per_node:
                frac_vars_now = _check_integer_feasibility(lp, cur_sol)
                if not frac_vars_now:
                    break
                step_res = _generate_and_apply_cut_round(
                    lp, tab, var_definitions, var_integrality,
                    cut_index=global_cut_counter + 1,
                    cut_type=cut_type,
                    node_id=node.node_id
                )
                if step_res is None:
                    break
                global_cut_counter += 1
                node_cut_rounds += 1
                dual_status, cut_detail = step_res
                node.cuts_generated.append(cut_detail)
                node.active_cuts.append(cut_detail)
                global_cuts_history.append(cut_detail)

                if dual_status == "INFEASIBLE":
                    node.lp_status = "INFEASIBLE"
                    node.outcome = "INFEASIBLE"
                    node.outcome_reason = (
                        f"Cut {cut_detail.cut_index} ({cut_detail.orig_vars_inequality_str}) "
                        "eliminated all remaining LP feasibility in this subtree (Dual Simplex INFEASIBLE)."
                    )
                    node.incumbent_at_step = incumbent_val
                    node.incumbent_sol_at_step = dict(incumbent_sol) if incumbent_sol else None
                    return
                if dual_status != "OPTIMAL":
                    node.lp_status = dual_status
                    node.outcome = dual_status
                    node.outcome_reason = f"Dual Simplex ended with status {dual_status} after Cut {cut_detail.cut_index}."
                    node.incumbent_at_step = incumbent_val
                    node.incumbent_sol_at_step = dict(incumbent_sol) if incumbent_sol else None
                    return

                cur_sol, cur_val = _extract_sol_from_tableau(lp, tab)
                node.lp_val = cur_val
                node.lp_sol = dict(cur_sol)

                if _is_bound_dominated(cur_val, incumbent_val, lp.obj_type):
                    node.outcome = "PRUNED_BOUND"
                    comp = "<=" if lp.obj_type == "max" else ">="
                    node.outcome_reason = (
                        f"After Cut {cut_detail.cut_index}, bound {format_fraction(cur_val)} {comp} "
                        f"Incumbent {format_fraction(incumbent_val)}."
                    )
                    node.incumbent_at_step = incumbent_val
                    node.incumbent_sol_at_step = dict(incumbent_sol) if incumbent_sol else None
                    return

        # Integrality check on node optimum (after any cuts)
        frac_vars = _check_integer_feasibility(lp, cur_sol)
        if not frac_vars:
            if _is_better_obj(cur_val, incumbent_val, lp.obj_type):
                incumbent_val = cur_val
                incumbent_sol = dict(cur_sol)
                incumbent_node_id = node.node_id
                node.outcome = "INTEGER_INCUMBENT"
                node.outcome_reason = f"New incumbent integer solution found: z = {format_fraction(cur_val)}."
            else:
                node.outcome = "INTEGER_FATHOMED"
                node.outcome_reason = (
                    f"Integer feasible (z = {format_fraction(cur_val)}), dominated by incumbent "
                    f"{format_fraction(incumbent_val)}."
                )
            node.incumbent_at_step = incumbent_val
            node.incumbent_sol_at_step = dict(incumbent_sol) if incumbent_sol else None
            return

        # Fractional candidate eligible for branching
        _, b_var, b_val = _select_branch_var(frac_vars, branch_rule)
        node.branch_var = b_var
        node.branch_val = b_val
        node.outcome = "BRANCH_CANDIDATE"
        node.incumbent_at_step = incumbent_val
        node.incumbent_sol_at_step = dict(incumbent_sol) if incumbent_sol else None
        open_queue.append(node)

    def pop_next_node() -> BBNode:
        if strategy == "dfs":
            return open_queue.pop()
        if strategy == "bfs":
            return open_queue.pop(0)
        # best-bound: max lp_val for max, min lp_val for min (tie-break lowest node_id)
        if lp.obj_type == "max":
            best_idx = max(range(len(open_queue)), key=lambda i: (open_queue[i].lp_val, -open_queue[i].node_id))
        else:
            best_idx = min(range(len(open_queue)), key=lambda i: (open_queue[i].lp_val, open_queue[i].node_id))
        return open_queue.pop(best_idx)

    # Evaluate Root Node 0
    root_node = BBNode(
        node_id=0,
        parent_id=None,
        depth=0,
        branch_desc="Root Relaxation",
        var_bounds=init_bounds,
        inherited_cuts=[]
    )
    node_counter = 1
    nodes.append(root_node)
    evaluate_node(root_node)

    hit_node_limit = False

    while open_queue:
        if len(nodes) >= max_nodes:
            hit_node_limit = True
            for rem in open_queue:
                rem.outcome = "NODE_LIMIT"
                rem.outcome_reason = f"Reached max node limit ({max_nodes})."
            break

        parent = pop_next_node()

        # Re-check bound dominance before expanding children (incumbent may have improved since parent was queued)
        if _is_bound_dominated(parent.lp_val, incumbent_val, lp.obj_type):
            parent.outcome = "PRUNED_BOUND"
            comp = "<=" if lp.obj_type == "max" else ">="
            parent.outcome_reason = (
                f"Pruned before branching: node bound {format_fraction(parent.lp_val)} {comp} "
                f"updated incumbent {format_fraction(incumbent_val)}."
            )
            continue

        b_var = parent.branch_var
        b_val = parent.branch_val
        assert b_var is not None and b_val is not None
        fl = floor_frac(b_val)
        cl = fl + 1

        parent.outcome = "BRANCHED"
        parent.outcome_reason = f"Branched on {b_var} = {format_fraction(b_val)} ({b_var} <= {fl} | {b_var} >= {cl})."

        # Left child: b_var <= fl
        left_bounds = dict(parent.var_bounds)
        old_lb, old_ub = left_bounds[b_var]
        new_ub = fl if old_ub is None else min(old_ub, fl)
        left_bounds[b_var] = (old_lb, new_ub)
        left_node = BBNode(
            node_id=node_counter,
            parent_id=parent.node_id,
            depth=parent.depth + 1,
            branch_desc=f"{b_var} <= {fl}",
            var_bounds=left_bounds,
            inherited_cuts=parent.active_cuts,
            branch_line=(b_var, "<=", fl)
        )
        node_counter += 1

        # Right child: b_var >= cl
        right_bounds = dict(parent.var_bounds)
        r_lb, r_ub = right_bounds[b_var]
        new_lb = max(r_lb, cl)
        right_bounds[b_var] = (new_lb, r_ub)
        right_node = BBNode(
            node_id=node_counter,
            parent_id=parent.node_id,
            depth=parent.depth + 1,
            branch_desc=f"{b_var} >= {cl}",
            var_bounds=right_bounds,
            inherited_cuts=parent.active_cuts,
            branch_line=(b_var, ">=", cl)
        )
        node_counter += 1

        parent.children_ids = [left_node.node_id, right_node.node_id]

        # For DFS stack order: evaluate left then right; pushing [right, left] if both branch pops left next.
        # Record left & right sequentially in node list.
        nodes.append(left_node)
        evaluate_node(left_node)
        nodes.append(right_node)
        evaluate_node(right_node)

        if strategy == "dfs" and len(open_queue) >= 2 and open_queue[-2] is left_node and open_queue[-1] is right_node:
            # Ensure left child is popped before right child under LIFO stack
            open_queue[-2], open_queue[-1] = open_queue[-1], open_queue[-2]

    if root_status == "UNBOUNDED":
        final_status = "UNBOUNDED"
    elif hit_node_limit and incumbent_sol is None:
        final_status = "NODE_LIMIT"
    elif incumbent_sol is None:
        final_status = "INFEASIBLE"
    else:
        final_status = "OPTIMAL"

    pruned_bound_count = sum(1 for n in nodes if n.outcome == "PRUNED_BOUND")
    infeasible_count = sum(1 for n in nodes if n.outcome == "INFEASIBLE")
    int_leaf_count = sum(1 for n in nodes if n.outcome in ("INTEGER_INCUMBENT", "INTEGER_FATHOMED"))
    branched_count = sum(1 for n in nodes if n.outcome == "BRANCHED")

    return {
        "method": "branch-and-cut" if enable_cuts else "branch-and-bound",
        "strategy": strategy,
        "branch_rule": branch_rule,
        "cut_type": cut_type if enable_cuts else None,
        "status": final_status,
        "is_integer": (incumbent_sol is not None),
        "relax_solution": relax_solution,
        "relax_opt_val": relax_opt_val,
        "relax_tab_before": relax_tab_before,
        "relax_tab_after": relax_tab_after,
        "cuts_history": global_cuts_history,
        "total_cuts": len(global_cuts_history),
        "nodes": nodes,
        "total_nodes": len(nodes),
        "incumbent_node_id": incumbent_node_id,
        "stats": {
            "branched": branched_count,
            "integer_leaves": int_leaf_count,
            "pruned_by_bound": pruned_bound_count,
            "infeasible": infeasible_count,
        },
        "final_solution": incumbent_sol if incumbent_sol is not None else {v: Fraction(0) for v in lp.var_names},
        "final_opt_val": incumbent_val
    }


def solve_branch_and_bound(
    lp: LinearProgram,
    strategy: str = "best-bound",
    branch_rule: str = "most-fractional",
    max_nodes: int = 200
) -> Dict[str, Any]:
    """Solves an ILP/MILP using exact Branch-and-Bound."""
    return _solve_tree_search(
        lp,
        enable_cuts=False,
        strategy=strategy,
        branch_rule=branch_rule,
        max_nodes=max_nodes
    )


def solve_branch_and_cut(
    lp: LinearProgram,
    strategy: str = "best-bound",
    branch_rule: str = "most-fractional",
    cut_type: str = "auto",
    cuts_per_node: int = 2,
    cut_depth: int = 5,
    max_nodes: int = 200
) -> Dict[str, Any]:
    """Solves an ILP/MILP using exact Branch-and-Cut (cutting planes + Dual Simplex + Branch-and-Bound)."""
    return _solve_tree_search(
        lp,
        enable_cuts=True,
        strategy=strategy,
        branch_rule=branch_rule,
        cut_type=cut_type,
        cuts_per_node=cuts_per_node,
        cut_depth=cut_depth,
        max_nodes=max_nodes
    )


# ---------------------------------------------------------------------
# 2D GRAPHICAL VISUALIZER (SVG & INTERACTIVE HTML)
# ---------------------------------------------------------------------

def compute_polygon_2d(
    constraints: List[Tuple[List[Fraction], str, Fraction]],
    bbox: Tuple[float, float, float, float]
) -> List[Tuple[float, float]]:
    """
    Computes vertices of 2D convex polygon (or segment/point) defined by linear constraints and bounding box.
    bbox: (xmin, xmax, ymin, ymax)
    """
    xmin, xmax, ymin, ymax = bbox
    all_halfplanes = []

    all_halfplanes.append(([-1.0, 0.0], -xmin))  # x >= xmin
    all_halfplanes.append(([1.0, 0.0], xmax))    # x <= xmax
    all_halfplanes.append(([0.0, -1.0], -ymin))  # y >= ymin
    all_halfplanes.append(([0.0, 1.0], ymax))    # y <= ymax

    for coeffs, rel, rhs in constraints:
        c1, c2 = float(coeffs[0]), float(coeffs[1])
        r = float(rhs)
        if rel == "<=":
            all_halfplanes.append(([c1, c2], r))
        elif rel == ">=":
            all_halfplanes.append(([-c1, -c2], -r))
        elif rel == "=":
            all_halfplanes.append(([c1, c2], r))
            all_halfplanes.append(([-c1, -c2], -r))

    points = []
    n_hp = len(all_halfplanes)
    for i in range(n_hp):
        (a1, a2), b1 = all_halfplanes[i]
        for j in range(i + 1, n_hp):
            (c1, c2), b2 = all_halfplanes[j]
            det = a1 * c2 - a2 * c1
            if abs(det) > 1e-9:
                px = (b1 * c2 - b2 * a2) / det
                py = (a1 * b2 - c1 * b1) / det

                feasible = True
                for (h1, h2), hb in all_halfplanes:
                    if h1 * px + h2 * py > hb + 1e-6:
                        feasible = False
                        break
                if feasible:
                    points.append((round(px, 5), round(py, 5)))

    unique_pts: List[Tuple[float, float]] = []
    for p in points:
        if not any(math.hypot(p[0] - u[0], p[1] - u[1]) < 1e-4 for u in unique_pts):
            unique_pts.append(p)

    if len(unique_pts) < 3:
        return unique_pts

    cx = sum(p[0] for p in unique_pts) / len(unique_pts)
    cy = sum(p[1] for p in unique_pts) / len(unique_pts)
    unique_pts.sort(key=lambda p: math.atan2(p[1] - cy, p[0] - cx))
    return unique_pts


def _compute_2d_bbox_and_grid(lp: LinearProgram) -> Tuple[List[Tuple[List[Fraction], str, Fraction]], Tuple[float, float, float, float], float, List[Dict[str, Any]], List[Dict[str, Any]]]:
    raw_constraints = lp.effective_constraints()
    probe_poly = compute_polygon_2d(raw_constraints, (0.0, 50.0, 0.0, 50.0))
    if probe_poly:
        max_x = max(max(p[0] for p in probe_poly) + 1.5, 5.0)
        max_y = max(max(p[1] for p in probe_poly) + 1.5, 5.0)
    else:
        max_x, max_y = 10.0, 10.0
    max_coord = max(max_x, max_y)
    bbox = (0.0, max_coord, 0.0, max_coord)

    def _point_feasible(px: float, py: float) -> bool:
        if px < -1e-6 or py < -1e-6:
            return False
        for coeffs, rel, rhs in raw_constraints:
            val = float(coeffs[0]) * px + float(coeffs[1]) * py
            r = float(rhs)
            if rel == "<=" and val > r + 1e-6:
                return False
            elif rel == ">=" and val < r - 1e-6:
                return False
            elif rel == "=" and abs(val - r) > 1e-6:
                return False
        return True

    t0, t1 = lp.var_types[0], lp.var_types[1]
    int_points: List[Dict[str, Any]] = []
    mixed_segments: List[Dict[str, Any]] = []

    if t0 in ("integer", "binary") and t1 in ("integer", "binary"):
        for ix in range(int(max_coord) + 1):
            for iy in range(int(max_coord) + 1):
                int_points.append({"x": ix, "y": iy, "in_initial_poly": _point_feasible(float(ix), float(iy))})
    else:
        if t0 in ("integer", "binary"):
            for ix in range(int(max_coord) + 1):
                seg_constrs = list(raw_constraints) + [([Fraction(1), Fraction(0)], "=", Fraction(ix))]
                pts = compute_polygon_2d(seg_constrs, bbox)
                if len(pts) >= 2:
                    pts.sort(key=lambda p: p[1])
                    mixed_segments.append({"x1": pts[0][0], "y1": pts[0][1], "x2": pts[-1][0], "y2": pts[-1][1]})
                elif len(pts) == 1:
                    mixed_segments.append({"x1": pts[0][0], "y1": pts[0][1], "x2": pts[0][0], "y2": pts[0][1]})
        elif t1 in ("integer", "binary"):
            for iy in range(int(max_coord) + 1):
                seg_constrs = list(raw_constraints) + [([Fraction(0), Fraction(1)], "=", Fraction(iy))]
                pts = compute_polygon_2d(seg_constrs, bbox)
                if len(pts) >= 2:
                    pts.sort(key=lambda p: p[0])
                    mixed_segments.append({"x1": pts[0][0], "y1": pts[0][1], "x2": pts[-1][0], "y2": pts[-1][1]})
                elif len(pts) == 1:
                    mixed_segments.append({"x1": pts[0][0], "y1": pts[0][1], "x2": pts[0][0], "y2": pts[0][1]})

    return raw_constraints, bbox, max_coord, int_points, mixed_segments


def _generate_2d_tree_visualizer(
    lp: LinearProgram,
    result: Dict[str, Any],
    tum_report: Dict[str, Any],
    output_html_path: str
) -> str:
    """Generates standalone interactive 2D HTML5+SVG visualizer for Branch-and-Bound and Branch-and-Cut."""
    raw_constraints, bbox, max_coord, int_points, mixed_segments = _compute_2d_bbox_and_grid(lp)
    init_poly = compute_polygon_2d(raw_constraints, bbox)

    v0, v1 = lp.var_names[0], lp.var_names[1]
    nodes_json = []
    for node in result["nodes"]:
        node_constrs = list(raw_constraints) + node.to_extra_constraints(lp)
        for c in node.cuts_generated:
            node_constrs.append(([c.orig_vars_coefficients[v0], c.orig_vars_coefficients[v1]], "<=", c.orig_vars_rhs))
        poly_post = compute_polygon_2d(node_constrs, bbox)

        branch_lines = []
        for idx_v, name in enumerate((v0, v1)):
            lb, ub = node.var_bounds[name]
            if lb > 0:
                branch_lines.append({"var_idx": idx_v, "var": name, "rel": ">=", "val": lb, "label": f"{name} >= {lb}"})
            def_ub = 1 if lp.var_types[idx_v] == "binary" else None
            if ub is not None and (def_ub is None or ub < def_ub):
                branch_lines.append({"var_idx": idx_v, "var": name, "rel": "<=", "val": ub, "label": f"{name} <= {ub}"})

        active_cuts_json = [
            {
                "index": c.cut_index,
                "kind": c.cut_kind,
                "node_id": c.node_id,
                "formula_slack": c.tableau_inequality_str,
                "formula_orig": c.orig_vars_inequality_str,
                "c1": float(c.orig_vars_coefficients[v0]),
                "c2": float(c.orig_vars_coefficients[v1]),
                "rhs": float(c.orig_vars_rhs)
            }
            for c in node.active_cuts
        ]

        nodes_json.append({
            "node_id": node.node_id,
            "parent_id": node.parent_id,
            "depth": node.depth,
            "branch_desc": node.branch_desc,
            "bounds_summary": node.bounds_summary(lp.var_names),
            "lp_status": node.lp_status,
            "lp_val_before": float(node.lp_val_before_cuts) if node.lp_val_before_cuts is not None else None,
            "lp_val_before_str": format_approx(node.lp_val_before_cuts) if node.lp_val_before_cuts is not None else None,
            "lp_val": float(node.lp_val) if node.lp_val is not None else None,
            "lp_val_str": format_approx(node.lp_val) if node.lp_val is not None else "N/A",
            "lp_sol": (
                {"x": float(node.lp_sol[v0]), "y": float(node.lp_sol[v1]), "x_str": format_fraction(node.lp_sol[v0]), "y_str": format_fraction(node.lp_sol[v1])}
                if node.lp_sol is not None else None
            ),
            "outcome": node.outcome,
            "outcome_reason": node.outcome_reason,
            "children_ids": node.children_ids,
            "poly": poly_post,
            "branch_lines": branch_lines,
            "cuts_generated_count": len(node.cuts_generated),
            "active_cuts": active_cuts_json,
            "incumbent_val_str": format_approx(node.incumbent_at_step) if node.incumbent_at_step is not None else "None",
            "incumbent_pt": (
                {"x": float(node.incumbent_sol_at_step[v0]), "y": float(node.incumbent_sol_at_step[v1]),
                 "x_str": format_fraction(node.incumbent_sol_at_step[v0]), "y_str": format_fraction(node.incumbent_sol_at_step[v1])}
                if node.incumbent_sol_at_step is not None else None
            )
        })

    relax_x = float(result["relax_solution"][v0])
    relax_y = float(result["relax_solution"][v1])
    relax_opt = float(result["relax_opt_val"]) if result["relax_opt_val"] is not None else 0.0
    has_final = result["final_opt_val"] is not None
    final_x = float(result["final_solution"][v0]) if has_final else 0.0
    final_y = float(result["final_solution"][v1]) if has_final else 0.0
    final_opt = float(result["final_opt_val"]) if has_final else 0.0

    method_title = "Branch & Cut" if result["method"] == "branch-and-cut" else "Branch & Bound"
    obj_badge = f"{lp.obj_type.upper()} {format_fraction(lp.obj_coeffs[0])}*{v0} + {format_fraction(lp.obj_coeffs[1])}*{v1}"

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>2D {method_title} Visualizer</title>
  <style>
    :root {{
      --primary: #2563eb;
      --bg: #0f172a;
      --card: #1e293b;
      --border: #334155;
      --text: #f8fafc;
      --text-muted: #94a3b8;
      --green: #22c55e;
      --red: #ef4444;
      --amber: #f59e0b;
      --cyan: #38bdf8;
      --purple: #a855f7;
    }}
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      background: var(--bg);
      color: var(--text);
      display: flex;
      flex-direction: column;
      height: 100vh;
      overflow: hidden;
    }}
    header {{
      background: #111827;
      padding: 12px 24px;
      border-bottom: 1px solid var(--border);
      display: flex;
      justify-content: space-between;
      align-items: center;
    }}
    h1 {{ font-size: 17px; font-weight: 700; color: #fff; }}
    .subtitle {{ font-size: 12px; color: var(--text-muted); }}
    .main-layout {{
      display: flex;
      flex: 1;
      height: calc(100vh - 54px);
      overflow: hidden;
    }}
    .canvas-panel {{
      flex: 1;
      background: #020617;
      display: flex;
      align-items: center;
      justify-content: center;
      padding: 16px;
    }}
    .sidebar {{
      width: 480px;
      background: var(--card);
      border-left: 1px solid var(--border);
      display: flex;
      flex-direction: column;
      overflow-y: auto;
    }}
    .section {{
      padding: 14px 18px;
      border-bottom: 1px solid var(--border);
    }}
    .section-title {{
      font-size: 12px;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.5px;
      color: var(--text-muted);
      margin-bottom: 8px;
    }}
    .badge {{
      display: inline-block;
      padding: 2px 7px;
      border-radius: 4px;
      font-size: 11px;
      font-weight: 700;
    }}
    .badge-green {{ background: rgba(34, 197, 94, 0.15); color: #4ade80; border: 1px solid #16a34a; }}
    .badge-red {{ background: rgba(239, 68, 68, 0.15); color: #f87171; border: 1px solid #dc2626; }}
    .badge-amber {{ background: rgba(245, 158, 11, 0.15); color: #fcd34d; border: 1px solid #d97706; }}
    .badge-cyan {{ background: rgba(56, 189, 248, 0.15); color: #38bdf8; border: 1px solid #0284c7; }}
    .badge-purple {{ background: rgba(168, 85, 247, 0.15); color: #d8b4fe; border: 1px solid #9333ea; }}
    .btn {{
      background: var(--primary);
      color: white;
      border: none;
      padding: 7px 13px;
      border-radius: 6px;
      font-size: 13px;
      font-weight: 600;
      cursor: pointer;
    }}
    .btn:disabled {{ opacity: 0.4; cursor: not-allowed; }}
    .btn-secondary {{ background: #334155; color: #e2e8f0; }}
    .btn-group {{ display: flex; gap: 8px; margin-top: 8px; }}
    .card-box {{
      background: rgba(15, 23, 42, 0.65);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 10px 12px;
      font-size: 13px;
      line-height: 1.45;
    }}
    .code-box {{
      font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
      background: #0f172a;
      padding: 5px 8px;
      border-radius: 5px;
      color: #38bdf8;
      font-size: 12px;
      margin-top: 4px;
      border: 1px solid #1e293b;
    }}
    .tree-list {{
      max-height: 220px;
      overflow-y: auto;
      display: flex;
      flex-direction: column;
      gap: 4px;
    }}
    .tree-node {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      padding: 6px 10px;
      border-radius: 6px;
      border: 1px solid var(--border);
      background: #0f172a;
      cursor: pointer;
      font-size: 12px;
    }}
    .tree-node.active {{
      border-color: var(--cyan);
      background: rgba(56, 189, 248, 0.12);
    }}
    svg {{ width: 100%; height: 100%; max-width: 800px; max-height: 800px; }}
  </style>
</head>
<body>
  <header>
    <div>
      <h1>Integer / Mixed-Integer LP: {method_title} 2D Visualizer</h1>
      <div class="subtitle">Strategy: {result['strategy']} | Branch Rule: {result['branch_rule']} | Nodes: {result['total_nodes']} | Cuts: {result['total_cuts']}</div>
    </div>
    <div>
      <span class="badge badge-purple">{obj_badge}</span>
    </div>
  </header>

  <div class="main-layout">
    <div class="canvas-panel">
      <svg id="plot-svg" viewBox="-50 -50 700 700"></svg>
    </div>

    <div class="sidebar">
      <div class="section">
        <div class="section-title">Pre-Solve: Total Unimodularity (TUM)</div>
        <div class="card-box">
          <strong>Status:</strong>
          <span class="badge {'badge-green' if tum_report['is_tum'] else 'badge-amber'}">
            {'IDEAL POLYHEDRON (TUM)' if tum_report['is_ideal'] else ('TUM' if tum_report['is_tum'] else 'NOT CERTIFIED TUM')}
          </span>
        </div>
      </div>

      <div class="section">
        <div class="section-title">Tree Step Navigator</div>
        <div style="display: flex; justify-content: space-between; align-items: center;">
          <span id="step-label" style="font-weight: 700; font-size: 14px; color: #38bdf8;">Node 0</span>
          <span id="step-counter" class="badge badge-amber">Node 0 of {max(len(result['nodes']) - 1, 0)}</span>
        </div>
        <div class="btn-group">
          <button id="btn-prev" class="btn btn-secondary" onclick="prevNode()">Previous</button>
          <button id="btn-next" class="btn" onclick="nextNode()">Next Node</button>
          <button class="btn btn-secondary" onclick="selectNode(0)">Reset</button>
        </div>
      </div>

      <div class="section">
        <div class="section-title">Search Tree Nodes (Click to Inspect)</div>
        <div id="tree-list" class="tree-list"></div>
      </div>

      <div class="section">
        <div class="section-title">Selected Node Details</div>
        <div id="node-detail" class="card-box"></div>
      </div>

      <div class="section">
        <div class="section-title">Final Solution Summary</div>
        <div class="card-box">
          <div style="color: #f87171;"><strong>Root LP ({lp.obj_type.upper()}):</strong> <code>{format_fraction(result['relax_opt_val']) if result['relax_opt_val'] is not None else 'N/A'}</code> at ({v0}={format_fraction(result['relax_solution'][v0])}, {v1}={format_fraction(result['relax_solution'][v1])})</div>
          <div style="color: #4ade80; margin-top: 4px;"><strong>Optimal {('MILP' if len(lp.integer_var_names()) < lp.num_vars else 'ILP')} ({result['status']}):</strong> <code>{format_fraction(result['final_opt_val']) if has_final else 'N/A'}</code> at ({v0}={format_fraction(result['final_solution'][v0])}, {v1}={format_fraction(result['final_solution'][v1])})</div>
        </div>
      </div>
    </div>
  </div>

  <script>
    const maxCoord = {max_coord};
    const svgSize = 600;
    const originX = 50;
    const originY = 550;
    const scale = (svgSize - 100) / maxCoord;
    const v0Name = {json.dumps(v0)};
    const v1Name = {json.dumps(v1)};

    function toSvgX(x) {{ return originX + x * scale; }}
    function toSvgY(y) {{ return originY - y * scale; }}

    const initPoly = {json.dumps(init_poly)};
    const intPoints = {json.dumps(int_points)};
    const mixedSegments = {json.dumps(mixed_segments)};
    const nodesData = {json.dumps(nodes_json)};
    const lpOpt = {{ x: {relax_x}, y: {relax_y}, val: {relax_opt} }};
    const hasFinal = {json.dumps(has_final)};
    const finalOpt = {{ x: {final_x}, y: {final_y}, val: {final_opt} }};

    let currentIdx = 0;

    function outcomeBadgeClass(out) {{
      if (out === 'INTEGER_INCUMBENT') return 'badge-green';
      if (out === 'INFEASIBLE' || out === 'UNBOUNDED') return 'badge-red';
      if (out === 'PRUNED_BOUND' || out === 'INTEGER_FATHOMED') return 'badge-amber';
      if (out === 'BRANCHED') return 'badge-cyan';
      return 'badge-purple';
    }}

    function renderTreeList() {{
      const el = document.getElementById('tree-list');
      el.innerHTML = nodesData.map((n, idx) => {{
        const indent = Math.min(n.depth, 6) * 14;
        const cls = outcomeBadgeClass(n.outcome);
        return `
          <div class="tree-node ${{idx === currentIdx ? 'active' : ''}}" onclick="selectNode(${{idx}})" style="padding-left: ${{10 + indent}}px;">
            <div>
              <strong>N${{n.node_id}}</strong>: ${{n.branch_desc}}
              <span style="color: #94a3b8; margin-left: 6px;">${{n.lp_status === 'OPTIMAL' ? 'z=' + n.lp_val_str : n.lp_status}}</span>
            </div>
            <span class="badge ${{cls}}">${{n.outcome}}</span>
          </div>
        `;
      }}).join('');
    }}

    function renderPlot() {{
      const svg = document.getElementById('plot-svg');
      const node = nodesData[currentIdx];
      let html = '';

      for (let i = 0; i <= maxCoord; i++) {{
        const sx = toSvgX(i);
        const sy = toSvgY(i);
        html += `<line x1="${{sx}}" y1="${{toSvgY(0)}}" x2="${{sx}}" y2="${{toSvgY(maxCoord)}}" stroke="#1e293b" stroke-width="1" />`;
        html += `<line x1="${{toSvgX(0)}}" y1="${{sy}}" x2="${{toSvgX(maxCoord)}}" y2="${{sy}}" stroke="#1e293b" stroke-width="1" />`;
        html += `<text x="${{sx}}" y="${{originY + 18}}" fill="#64748b" font-size="11" text-anchor="middle">${{i}}</text>`;
        html += `<text x="${{originX - 12}}" y="${{sy + 4}}" fill="#64748b" font-size="11" text-anchor="end">${{i}}</text>`;
      }}

      html += `<line x1="${{toSvgX(0)}}" y1="${{originY}}" x2="${{toSvgX(maxCoord)}}" y2="${{originY}}" stroke="#94a3b8" stroke-width="2" />`;
      html += `<line x1="${{originX}}" y1="${{toSvgY(0)}}" x2="${{originX}}" y2="${{toSvgY(maxCoord)}}" stroke="#94a3b8" stroke-width="2" />`;
      html += `<text x="${{toSvgX(maxCoord) + 12}}" y="${{originY + 4}}" fill="#94a3b8" font-size="13" font-weight="bold">${{v0Name}}</text>`;
      html += `<text x="${{originX}}" y="${{toSvgY(maxCoord) - 12}}" fill="#94a3b8" font-size="13" font-weight="bold" text-anchor="middle">${{v1Name}}</text>`;

      if (initPoly.length >= 3) {{
        const initPts = initPoly.map(p => `${{toSvgX(p[0])}},${{toSvgY(p[1])}}`).join(' ');
        html += `<polygon points="${{initPts}}" fill="rgba(37, 99, 235, 0.08)" stroke="#3b82f6" stroke-dasharray="4,4" stroke-width="1.5" />`;
      }}

      if (node.poly.length >= 3) {{
        const curPts = node.poly.map(p => `${{toSvgX(p[0])}},${{toSvgY(p[1])}}`).join(' ');
        html += `<polygon points="${{curPts}}" fill="rgba(56, 189, 248, 0.24)" stroke="#38bdf8" stroke-width="2.5" />`;
      }} else if (node.poly.length === 2) {{
        html += `<line x1="${{toSvgX(node.poly[0][0])}}" y1="${{toSvgY(node.poly[0][1])}}" x2="${{toSvgX(node.poly[1][0])}}" y2="${{toSvgY(node.poly[1][1])}}" stroke="#38bdf8" stroke-width="4" />`;
      }} else if (node.poly.length === 1) {{
        html += `<circle cx="${{toSvgX(node.poly[0][0])}}" cy="${{toSvgY(node.poly[0][1])}}" r="5" fill="#38bdf8" />`;
      }}

      mixedSegments.forEach(seg => {{
        if (seg.x1 === seg.x2 && seg.y1 === seg.y2) {{
          html += `<circle cx="${{toSvgX(seg.x1)}}" cy="${{toSvgY(seg.y1)}}" r="3" fill="#38bdf8" />`;
        }} else {{
          html += `<line x1="${{toSvgX(seg.x1)}}" y1="${{toSvgY(seg.y1)}}" x2="${{toSvgX(seg.x2)}}" y2="${{toSvgY(seg.y2)}}" stroke="#38bdf8" stroke-width="2" opacity="0.65" />`;
        }}
      }});

      intPoints.forEach(pt => {{
        const isFeas = pt.in_initial_poly;
        const fill = isFeas ? '#38bdf8' : '#475569';
        const r = isFeas ? 3.5 : 2;
        html += `<circle cx="${{toSvgX(pt.x)}}" cy="${{toSvgY(pt.y)}}" r="${{r}}" fill="${{fill}}" opacity="${{isFeas ? 1 : 0.4}}" />`;
      }});

      node.branch_lines.forEach(bl => {{
        if (bl.var_idx === 0) {{
          const sx = toSvgX(bl.val);
          html += `<line x1="${{sx}}" y1="${{toSvgY(0)}}" x2="${{sx}}" y2="${{toSvgY(maxCoord)}}" stroke="#a855f7" stroke-dasharray="6,4" stroke-width="2" />`;
          html += `<text x="${{sx + 4}}" y="${{toSvgY(maxCoord) + 14}}" fill="#d8b4fe" font-size="11" font-weight="bold">${{bl.label}}</text>`;
        }} else {{
          const sy = toSvgY(bl.val);
          html += `<line x1="${{toSvgX(0)}}" y1="${{sy}}" x2="${{toSvgX(maxCoord)}}" y2="${{sy}}" stroke="#a855f7" stroke-dasharray="6,4" stroke-width="2" />`;
          html += `<text x="${{toSvgX(maxCoord) - 55}}" y="${{sy - 6}}" fill="#d8b4fe" font-size="11" font-weight="bold">${{bl.label}}</text>`;
        }}
      }});

      node.active_cuts.forEach(cut => {{
        let p1, p2;
        if (Math.abs(cut.c2) > 1e-6) {{
          p1 = {{ x: 0, y: cut.rhs / cut.c2 }};
          p2 = {{ x: maxCoord, y: (cut.rhs - cut.c1 * maxCoord) / cut.c2 }};
        }} else if (Math.abs(cut.c1) > 1e-6) {{
          p1 = {{ x: cut.rhs / cut.c1, y: 0 }};
          p2 = {{ x: cut.rhs / cut.c1, y: maxCoord }};
        }}
        if (p1 && p2) {{
          html += `<line x1="${{toSvgX(p1.x)}}" y1="${{toSvgY(p1.y)}}" x2="${{toSvgX(p2.x)}}" y2="${{toSvgY(p2.y)}}" stroke="#f59e0b" stroke-width="2.2" />`;
          const midX = (p1.x + p2.x) / 2;
          const midY = (p1.y + p2.y) / 2;
          if (midX >= 0 && midX <= maxCoord && midY >= 0 && midY <= maxCoord) {{
            html += `<text x="${{toSvgX(midX)}}" y="${{toSvgY(midY) - 6}}" fill="#f59e0b" font-size="11" font-weight="bold">Cut ${{cut.index}}</text>`;
          }}
        }}
      }});

      html += `<circle cx="${{toSvgX(lpOpt.x)}}" cy="${{toSvgY(lpOpt.y)}}" r="5" fill="#ef4444" stroke="#ffffff" stroke-width="1.5" />`;
      html += `<text x="${{toSvgX(lpOpt.x) + 8}}" y="${{toSvgY(lpOpt.y) - 8}}" fill="#f87171" font-size="11" font-weight="bold">x*_LP (${{lpOpt.x.toFixed(2)}}, ${{lpOpt.y.toFixed(2)}})</text>`;

      if (node.lp_sol) {{
        html += `<circle cx="${{toSvgX(node.lp_sol.x)}}" cy="${{toSvgY(node.lp_sol.y)}}" r="6.5" fill="none" stroke="#38bdf8" stroke-width="2.5" />`;
      }}

      if (node.incumbent_pt) {{
        html += `<circle cx="${{toSvgX(node.incumbent_pt.x)}}" cy="${{toSvgY(node.incumbent_pt.y)}}" r="7" fill="#22c55e" stroke="#ffffff" stroke-width="2" />`;
        html += `<text x="${{toSvgX(node.incumbent_pt.x) + 8}}" y="${{toSvgY(node.incumbent_pt.y) + 16}}" fill="#4ade80" font-size="11" font-weight="bold">Incumbent (${{node.incumbent_pt.x_str}}, ${{node.incumbent_pt.y_str}})</text>`;
      }}

      svg.innerHTML = html;

      document.getElementById('step-label').textContent = `Node ${{node.node_id}} (${{node.branch_desc}})`;
      document.getElementById('step-counter').textContent = `Step ${{currentIdx}} of ${{Math.max(nodesData.length - 1, 0)}}`;
      document.getElementById('btn-prev').disabled = (currentIdx === 0);
      document.getElementById('btn-next').disabled = (currentIdx === nodesData.length - 1);

      let detailHtml = `
        <div><strong>Domain Bounds:</strong> <code>${{node.bounds_summary}}</code></div>
        <div style="margin-top: 6px;"><strong>Node LP Relaxation:</strong> ${{node.lp_status}}${{node.lp_val_before_str ? ' | Before cuts: ' + node.lp_val_before_str : ''}}</div>
        ${{node.lp_sol ? `<div><strong>Node Optimum:</strong> <code>${{v0Name}}=${{node.lp_sol.x_str}}, ${{v1Name}}=${{node.lp_sol.y_str}}</code> (z = ${{node.lp_val_str}})</div>` : ''}}
        <div style="margin-top: 6px;"><strong>Outcome:</strong> <span class="badge ${{outcomeBadgeClass(node.outcome)}}">${{node.outcome}}</span></div>
        <div style="color: #94a3b8; font-size: 12px; margin-top: 4px;">${{node.outcome_reason}}</div>
        <div style="margin-top: 6px;"><strong>Best Incumbent at Step:</strong> <code>${{node.incumbent_val_str}}</code></div>
      `;
      if (node.active_cuts.length > 0) {{
        detailHtml += `<div style="margin-top: 8px;"><strong>Active Cuts (${{node.active_cuts.length}}):</strong></div>`;
        node.active_cuts.forEach(c => {{
          detailHtml += `<div class="code-box">[${{c.kind}} Cut ${{c.index}} @ N${{c.node_id}}] ${{c.formula_orig}}</div>`;
        }});
      }}
      document.getElementById('node-detail').innerHTML = detailHtml;
      renderTreeList();
    }}

    function selectNode(idx) {{ currentIdx = idx; renderPlot(); }}
    function nextNode() {{ if (currentIdx < nodesData.length - 1) {{ currentIdx++; renderPlot(); }} }}
    function prevNode() {{ if (currentIdx > 0) {{ currentIdx--; renderPlot(); }} }}

    renderPlot();
  </script>
</body>
</html>
"""
    with open(output_html_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    return output_html_path


def generate_2d_html_visualizer(
    lp: LinearProgram,
    result: Dict[str, Any],
    tum_report: Dict[str, Any],
    output_html_path: str
) -> str:
    """
    Generates a standalone, interactive HTML5 + SVG visualizer for 2-variable (2D) problems across:
    - Chvatal-Gomory / GMI Cutting Planes
    - Branch and Bound
    - Branch and Cut
    """
    if not lp.is_2d():
        raise ValueError("Visualization is only supported for 2-variable (2D) problems.")

    if result.get("method") in ("branch-and-bound", "branch-and-cut"):
        return _generate_2d_tree_visualizer(lp, result, tum_report, output_html_path)

    raw_constraints, bbox, max_coord, int_points, mixed_segments = _compute_2d_bbox_and_grid(lp)

    polygons_by_step = []
    polygons_by_step.append({
        "step_name": "Continuous Relaxation (No Cuts)",
        "poly": compute_polygon_2d(raw_constraints, bbox)
    })

    current_constraints = list(raw_constraints)
    for cut in result["cuts_history"]:
        coeffs = [cut.orig_vars_coefficients[lp.var_names[0]], cut.orig_vars_coefficients[lp.var_names[1]]]
        current_constraints.append((coeffs, "<=", cut.orig_vars_rhs))
        polygons_by_step.append({
            "step_name": f"Cut {cut.cut_index}: {cut.orig_vars_inequality_str}",
            "poly": compute_polygon_2d(current_constraints, bbox)
        })

    relax_x = float(result["relax_solution"][lp.var_names[0]])
    relax_y = float(result["relax_solution"][lp.var_names[1]])
    relax_opt = float(result["relax_opt_val"]) if result["relax_opt_val"] is not None else 0.0

    has_final = result["final_opt_val"] is not None
    final_x = float(result["final_solution"][lp.var_names[0]]) if has_final else 0.0
    final_y = float(result["final_solution"][lp.var_names[1]]) if has_final else 0.0
    final_opt = float(result["final_opt_val"]) if has_final else 0.0

    cuts_json = []
    for c in result["cuts_history"]:
        c1 = float(c.orig_vars_coefficients[lp.var_names[0]])
        c2 = float(c.orig_vars_coefficients[lp.var_names[1]])
        cr = float(c.orig_vars_rhs)
        cuts_json.append({
            "index": c.cut_index,
            "kind": c.cut_kind,
            "generating_var": c.generating_var,
            "formula_slack": c.tableau_inequality_str,
            "formula_orig": c.orig_vars_inequality_str,
            "c1": c1,
            "c2": c2,
            "rhs": cr
        })

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>2D Cutting Plane Visualizer: Chvátal-Gomory Cuts</title>
  <style>
    :root {{
      --primary: #2563eb;
      --bg: #0f172a;
      --card: #1e293b;
      --border: #334155;
      --text: #f8fafc;
      --text-muted: #94a3b8;
      --green: #22c55e;
      --red: #ef4444;
      --amber: #f59e0b;
      --purple: #a855f7;
    }}
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      background: var(--bg);
      color: var(--text);
      display: flex;
      flex-direction: column;
      height: 100vh;
      overflow: hidden;
    }}
    header {{
      background: #111827;
      padding: 12px 24px;
      border-bottom: 1px solid var(--border);
      display: flex;
      justify-content: space-between;
      align-items: center;
    }}
    h1 {{ font-size: 17px; font-weight: 700; color: #fff; }}
    .subtitle {{ font-size: 12px; color: var(--text-muted); }}
    .main-layout {{
      display: flex;
      flex: 1;
      height: calc(100vh - 54px);
      overflow: hidden;
    }}
    .canvas-panel {{
      flex: 1;
      position: relative;
      background: #020617;
      display: flex;
      align-items: center;
      justify-content: center;
      padding: 16px;
    }}
    .sidebar {{
      width: 440px;
      background: var(--card);
      border-left: 1px solid var(--border);
      display: flex;
      flex-direction: column;
      overflow-y: auto;
    }}
    .section {{
      padding: 18px 20px;
      border-bottom: 1px solid var(--border);
    }}
    .section-title {{
      font-size: 13px;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.5px;
      color: var(--text-muted);
      margin-bottom: 10px;
    }}
    .badge {{
      display: inline-block;
      padding: 3px 8px;
      border-radius: 4px;
      font-size: 11px;
      font-weight: 700;
    }}
    .badge-green {{ background: rgba(34, 197, 94, 0.15); color: #4ade80; border: 1px solid #16a34a; }}
    .badge-amber {{ background: rgba(245, 158, 11, 0.15); color: #fcd34d; border: 1px solid #d97706; }}
    .badge-purple {{ background: rgba(168, 85, 247, 0.15); color: #d8b4fe; border: 1px solid #9333ea; }}
    .btn {{
      background: var(--primary);
      color: white;
      border: none;
      padding: 8px 14px;
      border-radius: 6px;
      font-size: 13px;
      font-weight: 600;
      cursor: pointer;
      transition: opacity 0.15s;
    }}
    .btn:hover {{ opacity: 0.9; }}
    .btn-secondary {{
      background: #334155;
      color: #e2e8f0;
    }}
    .btn-group {{
      display: flex;
      gap: 8px;
      margin-top: 10px;
    }}
    .card-box {{
      background: rgba(15, 23, 42, 0.6);
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 12px 14px;
      font-size: 13px;
      line-height: 1.5;
    }}
    .code-box {{
      font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
      background: #0f172a;
      padding: 6px 10px;
      border-radius: 6px;
      color: #38bdf8;
      font-size: 12px;
      margin-top: 4px;
      border: 1px solid #1e293b;
    }}
    svg {{
      width: 100%;
      height: 100%;
      max-width: 800px;
      max-height: 800px;
    }}
  </style>
</head>
<body>
  <header>
    <div>
      <h1>Integer Linear Programming & Chvátal-Gomory 2D Visualizer</h1>
      <div class="subtitle">Cutting Plane Method Step-by-Step with Total Unimodularity (TUM) Analysis</div>
    </div>
    <div>
      <span class="badge badge-purple">{lp.obj_type.upper()} {lp.obj_coeffs[0]}*{lp.var_names[0]} + {lp.obj_coeffs[1]}*{lp.var_names[1]}</span>
    </div>
  </header>

  <div class="main-layout">
    <div class="canvas-panel">
      <svg id="plot-svg" viewBox="-50 -50 700 700"></svg>
    </div>

    <div class="sidebar">
      <div class="section">
        <div class="section-title">Pre-Solve: Total Unimodularity (TUM) Analysis</div>
        <div class="card-box">
          <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 6px;">
            <strong>Matrix A Property:</strong>
            <span class="badge {'badge-green' if tum_report['is_tum'] else 'badge-amber'}">
              {'TOTALLY UNIMODULAR (TUM)' if tum_report['is_tum'] else 'NOT CERTIFIED TUM'}
            </span>
          </div>
          <div><strong>Polyhedron Status:</strong> {'IDEAL POLYHEDRON (P = P_I)' if tum_report['is_ideal'] else 'Non-Ideal Polyhedron'}</div>
          <div style="color: var(--text-muted); font-size: 12px; margin-top: 6px;">
            {tum_report['details'][-1]}
          </div>
        </div>
      </div>

      <div class="section">
        <div class="section-title">Cutting Plane Stepper</div>
        <div style="display: flex; justify-content: space-between; align-items: center;">
          <span id="step-label" style="font-weight: 700; font-size: 14px; color: #38bdf8;">Continuous LP Relaxation</span>
          <span id="step-counter" class="badge badge-amber">Step 0 of {len(result['cuts_history'])}</span>
        </div>
        <div class="btn-group">
          <button id="btn-prev" class="btn btn-secondary" onclick="prevStep()">Previous</button>
          <button id="btn-next" class="btn" onclick="nextStep()">Next Cut</button>
          <button class="btn btn-secondary" onclick="resetStep()">Reset</button>
        </div>
      </div>

      <div class="section" id="cut-info-section">
        <div class="section-title">Current Cut Algebraic Derivation</div>
        <div id="cut-info-body" class="card-box">
          Select 'Next Cut' to visualize sequential Chvátal-Gomory planes.
        </div>
      </div>

      <div class="section">
        <div class="section-title">Solution Comparison</div>
        <div class="card-box" style="margin-bottom: 8px;">
          <div style="color: #f87171; font-weight: 700;">Continuous Relaxation Optimum (x*_LP):</div>
          <div>{lp.var_names[0]} = <code>{format_fraction(result['relax_solution'][lp.var_names[0]])}</code>, {lp.var_names[1]} = <code>{format_fraction(result['relax_solution'][lp.var_names[1]])}</code></div>
          <div>Objective ({lp.obj_type.upper()} z): <code>{format_fraction(result['relax_opt_val']) if result['relax_opt_val'] is not None else 'N/A'}</code> (~{relax_opt:.4g})</div>
        </div>
        <div class="card-box">
          <div style="color: #4ade80; font-weight: 700;">Final Integer Optimum (x*_ILP):</div>
          <div>{lp.var_names[0]} = <code>{format_fraction(result['final_solution'][lp.var_names[0]])}</code>, {lp.var_names[1]} = <code>{format_fraction(result['final_solution'][lp.var_names[1]])}</code></div>
          <div>Objective ({lp.obj_type.upper()} z*): <code>{format_fraction(result['final_opt_val']) if has_final else 'N/A'}</code> (~{final_opt:.4g})</div>
        </div>
      </div>
    </div>
  </div>

  <script>
    const maxCoord = {max_coord};
    const svgSize = 600;
    const originX = 50;
    const originY = 550;
    const scale = (svgSize - 100) / maxCoord;

    function toSvgX(x) {{ return originX + x * scale; }}
    function toSvgY(y) {{ return originY - y * scale; }}

    const polygons = {json.dumps(polygons_by_step)};
    const intPoints = {json.dumps(int_points)};
    const mixedSegments = {json.dumps(mixed_segments)};
    const cutsData = {json.dumps(cuts_json)};
    const lpOpt = {{ x: {relax_x}, y: {relax_y}, val: {relax_opt} }};
    const hasFinal = {json.dumps(has_final)};
    const ilpOpt = {{ x: {final_x}, y: {final_y}, val: {final_opt} }};

    let currentStep = 0;

    function renderPlot() {{
      const svg = document.getElementById('plot-svg');
      let html = '';

      for (let i = 0; i <= maxCoord; i++) {{
        const sx = toSvgX(i);
        const sy = toSvgY(i);
        html += `<line x1="${{sx}}" y1="${{toSvgY(0)}}" x2="${{sx}}" y2="${{toSvgY(maxCoord)}}" stroke="#1e293b" stroke-width="1" />`;
        html += `<line x1="${{toSvgX(0)}}" y1="${{sy}}" x2="${{toSvgX(maxCoord)}}" y2="${{sy}}" stroke="#1e293b" stroke-width="1" />`;
        html += `<text x="${{sx}}" y="${{originY + 18}}" fill="#64748b" font-size="11" text-anchor="middle">${{i}}</text>`;
        html += `<text x="${{originX - 12}}" y="${{sy + 4}}" fill="#64748b" font-size="11" text-anchor="end">${{i}}</text>`;
      }}

      html += `<line x1="${{toSvgX(0)}}" y1="${{originY}}" x2="${{toSvgX(maxCoord)}}" y2="${{originY}}" stroke="#94a3b8" stroke-width="2" />`;
      html += `<line x1="${{originX}}" y1="${{toSvgY(0)}}" x2="${{originX}}" y2="${{toSvgY(maxCoord)}}" stroke="#94a3b8" stroke-width="2" />`;
      html += `<text x="${{toSvgX(maxCoord) + 12}}" y="${{originY + 4}}" fill="#94a3b8" font-size="13" font-weight="bold">{lp.var_names[0]}</text>`;
      html += `<text x="${{originX}}" y="${{toSvgY(maxCoord) - 12}}" fill="#94a3b8" font-size="13" font-weight="bold" text-anchor="middle">{lp.var_names[1]}</text>`;

      if (polygons[0].poly.length >= 3) {{
        const initPts = polygons[0].poly.map(p => `${{toSvgX(p[0])}},${{toSvgY(p[1])}}`).join(' ');
        html += `<polygon points="${{initPts}}" fill="rgba(37, 99, 235, 0.08)" stroke="#3b82f6" stroke-dasharray="4,4" stroke-width="1.5" />`;
      }}

      const cur = polygons[currentStep].poly;
      if (cur.length >= 3) {{
        const curPts = cur.map(p => `${{toSvgX(p[0])}},${{toSvgY(p[1])}}`).join(' ');
        html += `<polygon points="${{curPts}}" fill="rgba(56, 189, 248, 0.22)" stroke="#38bdf8" stroke-width="2.5" />`;
      }} else if (cur.length === 2) {{
        html += `<line x1="${{toSvgX(cur[0][0])}}" y1="${{toSvgY(cur[0][1])}}" x2="${{toSvgX(cur[1][0])}}" y2="${{toSvgY(cur[1][1])}}" stroke="#38bdf8" stroke-width="3.5" />`;
      }} else if (cur.length === 1) {{
        html += `<circle cx="${{toSvgX(cur[0][0])}}" cy="${{toSvgY(cur[0][1])}}" r="5" fill="#38bdf8" />`;
      }}

      for (let c = 0; c < currentStep; c++) {{
        const cut = cutsData[c];
        let p1, p2;
        if (Math.abs(cut.c2) > 1e-6) {{
          p1 = {{ x: 0, y: cut.rhs / cut.c2 }};
          p2 = {{ x: maxCoord, y: (cut.rhs - cut.c1 * maxCoord) / cut.c2 }};
        }} else {{
          p1 = {{ x: cut.rhs / cut.c1, y: 0 }};
          p2 = {{ x: cut.rhs / cut.c1, y: maxCoord }};
        }}
        html += `<line x1="${{toSvgX(p1.x)}}" y1="${{toSvgY(p1.y)}}" x2="${{toSvgX(p2.x)}}" y2="${{toSvgY(p2.y)}}" stroke="#f59e0b" stroke-width="2.5" />`;
        const midX = (p1.x + p2.x) / 2;
        const midY = (p1.y + p2.y) / 2;
        if (midX >= 0 && midX <= maxCoord && midY >= 0 && midY <= maxCoord) {{
          html += `<text x="${{toSvgX(midX)}}" y="${{toSvgY(midY) - 8}}" fill="#f59e0b" font-size="11" font-weight="bold">Cut ${{cut.index}}</text>`;
        }}
      }}

      mixedSegments.forEach(seg => {{
        html += `<line x1="${{toSvgX(seg.x1)}}" y1="${{toSvgY(seg.y1)}}" x2="${{toSvgX(seg.x2)}}" y2="${{toSvgY(seg.y2)}}" stroke="#38bdf8" stroke-width="2" opacity="0.65" />`;
      }});

      intPoints.forEach(pt => {{
        const isFeas = pt.in_initial_poly;
        const fill = isFeas ? '#38bdf8' : '#475569';
        const r = isFeas ? 3.5 : 2;
        html += `<circle cx="${{toSvgX(pt.x)}}" cy="${{toSvgY(pt.y)}}" r="${{r}}" fill="${{fill}}" opacity="${{isFeas ? 1 : 0.4}}" />`;
      }});

      html += `<circle cx="${{toSvgX(lpOpt.x)}}" cy="${{toSvgY(lpOpt.y)}}" r="6" fill="#ef4444" stroke="#ffffff" stroke-width="2" />`;
      html += `<text x="${{toSvgX(lpOpt.x) + 8}}" y="${{toSvgY(lpOpt.y) - 8}}" fill="#f87171" font-size="12" font-weight="bold">x*_LP (${{lpOpt.x.toFixed(2)}}, ${{lpOpt.y.toFixed(2)}})</text>`;

      if (hasFinal && currentStep === cutsData.length) {{
        html += `<circle cx="${{toSvgX(ilpOpt.x)}}" cy="${{toSvgY(ilpOpt.y)}}" r="7" fill="#22c55e" stroke="#ffffff" stroke-width="2.5" />`;
        html += `<text x="${{toSvgX(ilpOpt.x) + 8}}" y="${{toSvgY(ilpOpt.y) + 16}}" fill="#4ade80" font-size="12" font-weight="bold">x*_ILP (${{ilpOpt.x}}, ${{ilpOpt.y}})</text>`;
      }}

      svg.innerHTML = html;
      updateSidebar();
    }}

    function updateSidebar() {{
      const total = cutsData.length;
      document.getElementById('step-counter').textContent = `Step ${{currentStep}} of ${{total}}`;
      document.getElementById('btn-prev').disabled = (currentStep === 0);
      document.getElementById('btn-next').disabled = (currentStep === total);

      if (currentStep === 0) {{
        document.getElementById('step-label').textContent = "Continuous LP Relaxation";
        document.getElementById('cut-info-body').innerHTML = `
          <div><strong>Base Problem:</strong> Continuous relaxation solved to optimality via Simplex.</div>
          <div style="margin-top: 6px;">Continuous optimum: <code>{lp.var_names[0]} = ${{lpOpt.x.toFixed(3)}}, {lp.var_names[1]} = ${{lpOpt.y.toFixed(3)}}</code>.</div>
          <div style="color: var(--amber); margin-top: 6px;">Chvátal-Gomory cuts slice off fractional solutions while preserving all feasible integer points.</div>
        `;
      }} else {{
        const cut = cutsData[currentStep - 1];
        document.getElementById('step-label').textContent = `${{cut.kind}} Cut ${{cut.index}} Applied`;
        document.getElementById('cut-info-body').innerHTML = `
          <div><strong>Generating Variable:</strong> <code>${{cut.generating_var}}</code> (${{cut.kind}})</div>
          <div style="margin-top: 6px;"><strong>Tableau Slack Form:</strong></div>
          <div class="code-box">${{cut.formula_slack}}</div>
          <div style="margin-top: 6px;"><strong>Original Decision Variables Form:</strong></div>
          <div class="code-box" style="color: #4ade80;">${{cut.formula_orig}}</div>
          <div style="margin-top: 6px; font-size: 12px; color: var(--text-muted);">
            The yellow hyperplane separates the fractional relaxation point without eliminating any integer feasible points.
          </div>
        `;
      }}
    }}

    function nextStep() {{ if (currentStep < cutsData.length) {{ currentStep++; renderPlot(); }} }}
    function prevStep() {{ if (currentStep > 0) {{ currentStep--; renderPlot(); }} }}
    function resetStep() {{ currentStep = 0; renderPlot(); }}

    renderPlot();
  </script>
</body>
</html>
"""
    with open(output_html_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    return output_html_path


# ---------------------------------------------------------------------
# LP PARSER & CLI UTILITIES
# ---------------------------------------------------------------------

def parse_lp_input(text: str) -> LinearProgram:
    """
    Parses a standard algebraic LP / ILP / MILP formulation string (multiline or single-line with ';'):
    min/max c1*x1 + c2*x2
    s.t.
    a11*x1 + a12*x2 <= b1
    x1, x2 >= 0, integer   (or binary / continuous)
    """
    raw_lines = [line.strip() for line in text.strip().splitlines() if line.strip() and not line.strip().startswith("#")]
    lines: List[str] = []
    for raw_line in raw_lines:
        if re.search(r'(?:s\.?t\.?:?|subject\s+to:?)', raw_line, re.I):
            parts = re.split(r'(?:s\.?t\.?:?|subject\s+to:?)', raw_line, maxsplit=1, flags=re.I)
            before = parts[0].strip()
            after = parts[1].strip() if len(parts) > 1 else ""
            if before:
                lines.extend([p.strip() for p in before.split(";") if p.strip()])
            if after:
                lines.extend([p.strip() for p in after.split(";") if p.strip()])
        elif ";" in raw_line:
            lines.extend([p.strip() for p in raw_line.split(";") if p.strip()])
        else:
            lines.append(raw_line)

    if len(lines) < 2:
        raise ValueError("LP formulation requires at least 1 objective line and 1 constraint line.")

    # 1. Objective
    obj_line = lines[0]
    m_obj = re.match(r'^(min(?:imize)?|max(?:imize)?)\s+(?:[a-zA-Z]\s*=\s*)?(.+)$', obj_line, re.I)
    if not m_obj:
        raise ValueError(f"Could not parse objective line: '{lines[0]}'")
    obj_type = "min" if "min" in m_obj.group(1).lower() else "max"
    obj_expr = m_obj.group(2)

    # 2. Collect constraint lines & optional variable type declarations
    raw_var_tokens = set(re.findall(r'[a-zA-Z]\w*', obj_expr))
    parsed_raw_constraints: List[Tuple[str, str, Fraction]] = []
    declared_types: Dict[str, str] = {}

    bound_type_re = re.compile(
        r'^([a-zA-Z]\w*(?:\s*,\s*[a-zA-Z]\w*)*)\s*>=\s*0(?:\s*,\s*(integer|int|binary|bin|continuous|cont|real))?$',
        re.I
    )
    prefix_type_re = re.compile(
        r'^(integer|int|binary|bin|continuous|cont|real)\s+([a-zA-Z]\w*(?:\s*,\s*[a-zA-Z]\w*)*)$',
        re.I
    )
    suffix_type_re = re.compile(
        r'^([a-zA-Z]\w*(?:\s*,\s*[a-zA-Z]\w*)*)\s+(integer|int|binary|bin|continuous|cont|real)$',
        re.I
    )

    def _normalize_vtype(raw_t: str) -> str:
        t = raw_t.lower()
        if t in ("int", "integer"):
            return "integer"
        if t in ("bin", "binary"):
            return "binary"
        return "continuous"

    for line_clean in lines[1:]:
        if not line_clean:
            continue

        m_bt = bound_type_re.match(line_clean)
        if m_bt:
            v_list = [v.strip() for v in m_bt.group(1).split(",") if v.strip()]
            raw_var_tokens.update(v_list)
            if m_bt.group(2):
                vt = _normalize_vtype(m_bt.group(2))
                for v in v_list:
                    declared_types[v] = vt
            continue

        m_pt = prefix_type_re.match(line_clean)
        if m_pt:
            vt = _normalize_vtype(m_pt.group(1))
            v_list = [v.strip() for v in m_pt.group(2).split(",") if v.strip()]
            raw_var_tokens.update(v_list)
            for v in v_list:
                declared_types[v] = vt
            continue

        m_st = suffix_type_re.match(line_clean)
        if m_st:
            v_list = [v.strip() for v in m_st.group(1).split(",") if v.strip()]
            vt = _normalize_vtype(m_st.group(2))
            raw_var_tokens.update(v_list)
            for v in v_list:
                declared_types[v] = vt
            continue

        m_rel = re.search(r'(<=|>=|=)', line_clean)
        if not m_rel:
            continue
        rel = m_rel.group(1)
        lhs, rhs = line_clean.split(rel, 1)
        raw_var_tokens.update(re.findall(r'[a-zA-Z]\w*', lhs))
        parsed_raw_constraints.append((lhs, rel, parse_fraction(rhs)))

    if not parsed_raw_constraints:
        raise ValueError("LP formulation requires at least 1 valid constraint line.")

    def var_sort_key(name: str):
        m = re.match(r'^([a-zA-Z]+)(\d+)?(.*)$', name)
        if m:
            prefix, num, rest = m.groups()
            return (prefix, int(num) if num else 0, rest)
        return (name, 0, "")

    var_tokens = sorted(raw_var_tokens, key=var_sort_key)
    if not var_tokens:
        raise ValueError("No decision variables found in LP formulation.")

    def parse_expr(expr: str, variables: List[str]) -> List[Fraction]:
        coeffs = {v: Fraction(0) for v in variables}
        clean_expr = expr.replace('-', '+-')
        terms = [t.strip() for t in clean_expr.split('+') if t.strip()]
        for term in terms:
            for v in variables:
                if re.search(r'\b' + re.escape(v) + r'\b', term):
                    c_str = re.sub(r'\s*\*?\s*\b' + re.escape(v) + r'\b', '', term)
                    c_str = re.sub(r'\s+', '', c_str)
                    if c_str == '' or c_str == '+':
                        coeffs[v] += Fraction(1)
                    elif c_str == '-':
                        coeffs[v] += Fraction(-1)
                    else:
                        coeffs[v] += parse_fraction(c_str)
                    break
        return [coeffs[v] for v in variables]

    obj_coeffs = parse_expr(obj_expr, var_tokens)
    constraints = [
        (parse_expr(lhs, var_tokens), rel, rhs_val)
        for lhs, rel, rhs_val in parsed_raw_constraints
    ]
    var_types = [declared_types.get(v, "integer") for v in var_tokens]

    return LinearProgram(obj_type, obj_coeffs, constraints, var_names=var_tokens, var_types=var_types)


def verify_with_scipy_milp(lp: LinearProgram, result: Dict[str, Any]) -> Dict[str, Any]:
    """Cross-verifies the ILP/MILP solution with SciPy HiGHS MILP solver."""
    try:
        from scipy.optimize import milp, LinearConstraint, Bounds
        import numpy as np
    except ImportError:
        return {
            "installed": False,
            "message": "SciPy is not installed in the active Python environment. Cross-check skipped."
        }

    try:
        c = [float(v) if lp.obj_type == "min" else -float(v) for v in lp.obj_coeffs]
        A_rows = []
        lhs_bounds = []
        rhs_bounds = []

        for coeffs, rel, rhs in lp.constraints:
            A_rows.append([float(val) for val in coeffs])
            r = float(rhs)
            if rel == "<=":
                lhs_bounds.append(-np.inf)
                rhs_bounds.append(r)
            elif rel == ">=":
                lhs_bounds.append(r)
                rhs_bounds.append(np.inf)
            elif rel == "=":
                lhs_bounds.append(r)
                rhs_bounds.append(r)

        constraints = LinearConstraint(A_rows, lhs_bounds, rhs_bounds)
        integrality = np.array([1 if vt in ("integer", "binary") else 0 for vt in lp.var_types])
        ub_arr = np.array([1.0 if vt == "binary" else np.inf for vt in lp.var_types])
        bounds = Bounds(lb=np.zeros(lp.num_vars), ub=ub_arr)

        res = milp(c=c, constraints=constraints, integrality=integrality, bounds=bounds)

        our_status = result.get("status")
        if our_status == "INFEASIBLE":
            match = (not res.success and res.status == 2)
            return {
                "installed": True,
                "success": True,
                "scipy_fun": None,
                "scipy_x": None,
                "match": match,
                "message": "MATCH CONFIRMED: Both solver and SciPy HiGHS report INFEASIBLE." if match else f"Discrepancy: SciPy reported status {res.status} ({res.message})."
            }

        if not res.success:
            return {
                "installed": True,
                "success": False,
                "message": res.message
            }

        scipy_x = {lp.var_names[j]: round(float(res.x[j]), 5) for j in range(lp.num_vars)}
        scipy_fun = float(res.fun) if lp.obj_type == "min" else -float(res.fun)

        if result.get("final_opt_val") is None:
            return {
                "installed": True,
                "success": True,
                "scipy_fun": scipy_fun,
                "scipy_x": scipy_x,
                "match": False,
                "message": f"Solver stopped ({our_status}) without optimal value while SciPy found {scipy_fun}."
            }

        our_opt = float(result["final_opt_val"])
        match = abs(our_opt - scipy_fun) < 1e-4

        method_label = {
            "cutting-plane": "Cutting Plane",
            "branch-and-bound": "Branch & Bound",
            "branch-and-cut": "Branch & Cut"
        }.get(result.get("method", ""), "Solver")

        return {
            "installed": True,
            "success": True,
            "scipy_fun": scipy_fun,
            "scipy_x": scipy_x,
            "match": match,
            "message": f"MATCH CONFIRMED: {method_label} ({our_opt}) == SciPy MILP ({scipy_fun})" if match else f"Discrepancy detected ({our_opt} vs {scipy_fun})."
        }
    except Exception as ex:
        return {
            "installed": True,
            "error": str(ex)
        }


# ---------------------------------------------------------------------
# TERMINAL & MARKDOWN REPORT FORMATTERS
# ---------------------------------------------------------------------

def generate_markdown_report(lp: LinearProgram, tum_report: Dict[str, Any], result: Dict[str, Any], scipy_res: Optional[Dict[str, Any]] = None) -> str:
    method = result.get("method", "cutting-plane")
    title_map = {
        "cutting-plane": "Integer Linear Programming Report: Chvátal-Gomory Cutting Planes",
        "branch-and-bound": "Integer Linear Programming Report: Branch and Bound",
        "branch-and-cut": "Integer Linear Programming Report: Branch and Cut"
    }
    lines: List[str] = [f"# {title_map.get(method, 'ILP Solver Report')}", ""]

    lines.append("## Formulation")
    obj_terms = " + ".join(f"{format_fraction(c)}*{v}" for c, v in zip(lp.obj_coeffs, lp.var_names))
    lines.append(f"- **Objective**: `${lp.obj_type.upper()} z = {obj_terms}$`")
    lines.append("- **Constraints**:")
    for coeffs, rel, rhs in lp.constraints:
        lhs = " + ".join(f"{format_fraction(c)}*{v}" for c, v in zip(coeffs, lp.var_names) if c != 0) or "0"
        lines.append(f"  - `{lhs} {rel} {format_fraction(rhs)}`")
    lines.append("- **Variables**: " + ", ".join(f"`{v}` >= 0 ({vt})" for v, vt in zip(lp.var_names, lp.var_types)))
    lines.append("")

    lines.append("## Pre-Solve Total Unimodularity (TUM) Analysis")
    for line in tum_report["details"]:
        lines.append(f"- {line}")
    lines.append("")

    lines.append("## Root LP Relaxation")
    lines.append(f"- **Continuous Objective**: `{format_approx(result['relax_opt_val']) if result['relax_opt_val'] is not None else 'N/A'}`")
    for v, val in result["relax_solution"].items():
        lines.append(f"- `{v}` = `{format_approx(val)}`")
    lines.append("")

    if method in ("branch-and-bound", "branch-and-cut"):
        lines.append(f"## Search Tree Trace ({result['total_nodes']} Nodes Evaluated)")
        lines.append(f"- **Node Selection Strategy**: `{result.get('strategy')}`")
        lines.append(f"- **Branching Rule**: `{result.get('branch_rule')}`")
        lines.append("")
        lines.append("| Node | Parent | Depth | Branch / Domain | LP Relaxation | Outcome | Notes |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- |")
        for n in result["nodes"]:
            p_str = str(n.parent_id) if n.parent_id is not None else "-"
            lp_str = format_approx(n.lp_val) if n.lp_val is not None else n.lp_status
            lines.append(f"| N{n.node_id} | {p_str} | {n.depth} | `{n.branch_desc}` | `{lp_str}` | `{n.outcome}` | {n.outcome_reason} |")
        lines.append("")

    if result["total_cuts"] > 0:
        lines.append(f"## Cutting Planes Applied ({result['total_cuts']})")
        for cut in result["cuts_history"]:
            node_tag = f" @ Node N{cut.node_id}" if method == "branch-and-cut" else ""
            lines.append(f"### Cut {cut.cut_index} [{cut.cut_kind}{node_tag}] (`{cut.generating_var}`)")
            lines.append(f"- **Tableau inequality**: `{cut.tableau_inequality_str}`")
            lines.append(f"- **Original variables**: `{cut.orig_vars_inequality_str}`")
        lines.append("")

    lines.append("## Final Optimal Solution")
    lines.append(f"- **Status**: `{result['status']}`")
    if result["final_opt_val"] is not None:
        lines.append(f"- **Objective Value ({lp.obj_type.upper()} z*)**: `{format_approx(result['final_opt_val'])}`")
        if result["relax_opt_val"] is not None:
            gap = abs(result["relax_opt_val"] - result["final_opt_val"])
            lines.append(f"- **Integrality Gap**: `{format_approx(gap)}`")
        for v, val in result["final_solution"].items():
            lines.append(f"- `{v}` = `{format_approx(val)}`")
    lines.append("")

    if scipy_res:
        lines.append("## SciPy MILP Verification")
        lines.append(f"- {scipy_res.get('message', 'Unavailable')}")
        lines.append("")

    return "\n".join(lines)


def print_terminal_report(lp: LinearProgram, tum_report: Dict[str, Any], result: Dict[str, Any], scipy_res: Optional[Dict[str, Any]] = None):
    """Prints a comprehensive, emoji-free terminal report for Cutting Plane, Branch & Bound, or Branch & Cut."""
    method = result.get("method", "cutting-plane")
    header_map = {
        "cutting-plane": "INTEGER LINEAR PROGRAMMING: CHVÁTAL-GOMORY CUTTING PLANES",
        "branch-and-bound": "INTEGER LINEAR PROGRAMMING: BRANCH AND BOUND",
        "branch-and-cut": "INTEGER LINEAR PROGRAMMING: BRANCH AND CUT"
    }

    print("=" * 82)
    print(header_map.get(method, "INTEGER LINEAR PROGRAMMING SOLVER").center(82))
    print("=" * 82)

    # 1. Formulation
    print("\n[1] PROBLEM FORMULATION:")
    obj_terms = [f"({format_fraction(c)})*{name}" for c, name in zip(lp.obj_coeffs, lp.var_names)]
    print(f"    {lp.obj_type.upper()} z = " + " + ".join(obj_terms))
    print("    Subject to:")
    for coeffs, rel, rhs in lp.constraints:
        lhs_terms = [f"({format_fraction(c)})*{name}" for c, name in zip(coeffs, lp.var_names) if c != 0]
        lhs_str = " + ".join(lhs_terms) if lhs_terms else "0"
        print(f"      {lhs_str} {rel} {format_fraction(rhs)}")
    print("      " + ", ".join(f"{v} >= 0 ({vt})" for v, vt in zip(lp.var_names, lp.var_types)))

    # 2. Total Unimodularity
    print("\n" + "=" * 82)
    print(" [2] TOTAL UNIMODULARITY (TUM) & IDEAL POLYHEDRON ANALYSIS".center(82))
    print("=" * 82)
    for line in tum_report["details"]:
        print(f" • {line}")

    # 3. Continuous LP Relaxation
    print("\n" + "=" * 82)
    print(" [3] ROOT CONTINUOUS LP RELAXATION (STEP 0)".center(82))
    print("=" * 82)
    if result.get("relax_tab_after"):
        print(result["relax_tab_after"])
    if result["relax_opt_val"] is not None:
        print(f"\n Relaxation Status: OPTIMAL")
        print(f" Continuous Objective ({lp.obj_type.upper()} z): {format_approx(result['relax_opt_val'])}")
        print(" Continuous Solution Vector:")
        for var, val in result["relax_solution"].items():
            print(f"   {var:<6} = {format_approx(val)}")
    else:
        print("\n Root LP Relaxation did not reach optimality (INFEASIBLE or UNBOUNDED).")

    # 4. Algorithm Trace (Cuts or Search Tree)
    if method == "cutting-plane":
        if result["total_cuts"] == 0:
            print("\n" + "=" * 82)
            print(" [4] CUTTING PLANE ALGORITHM: NO CUTS NEEDED".center(82))
            print("=" * 82)
            print(" All integer-constrained decision variables are already integers!")
        else:
            print("\n" + "=" * 82)
            print(f" [4] CUTTING PLANES APPLIED (TOTAL CUTS INJECTED: {result['total_cuts']})".center(82))
            print("=" * 82)

            for cut in result["cuts_history"]:
                print(f"\n>>> Cut {cut.cut_index} [{cut.cut_kind}]: Generated from Row {cut.generating_row} ('{cut.generating_var}')")
                print(f"    RHS Value: {format_fraction(cut.rhs_fractional)}  =>  Floor = {cut.floor_rhs}, Fractional Part f0 = {format_fraction(cut.frac_rhs)}")
                print(f"    Tableau Inequality:                 {cut.tableau_inequality_str}")
                print(f"    Re-projected to Decision Variables:  {cut.orig_vars_inequality_str}")
                print(f"    Injected Row in Dual Simplex:")
                print(f"      -s_cut{cut.cut_index} + sum(-g_j * x_j) = -{format_fraction(cut.frac_rhs)}")
                if cut.opt_val_before is not None and cut.opt_val_after is not None:
                    print(f"    Dual Simplex Pivot executed ({format_approx(cut.opt_val_before)} -> {format_approx(cut.opt_val_after)}).")
                else:
                    print(f"    Dual Simplex Pivot executed.")
    else:
        sec_title = " [4] BRANCH AND CUT SEARCH TREE" if method == "branch-and-cut" else " [4] BRANCH AND BOUND SEARCH TREE"
        print("\n" + "=" * 82)
        print(f"{sec_title} (Strategy: {result['strategy']}, Rule: {result['branch_rule']})".center(82))
        print("=" * 82)
        for n in result["nodes"]:
            parent_str = f"N{n.parent_id}" if n.parent_id is not None else "ROOT"
            indent = "  " * min(n.depth, 5)
            print(f"\n{indent}[Node N{n.node_id}] (Parent: {parent_str}, Depth: {n.depth}) | Rule: {n.branch_desc}")
            print(f"{indent}  Domain Bounds: {n.bounds_summary(lp.var_names)}")
            if n.lp_val_before_cuts is not None and n.cuts_generated:
                b_sol_str = ", ".join(f"{k}={format_fraction(v)}" for k, v in n.lp_sol_before_cuts.items())
                print(f"{indent}  LP before cuts: z = {format_approx(n.lp_val_before_cuts)} ({b_sol_str})")
                for cut in n.cuts_generated:
                    print(f"{indent}    -> Cut {cut.cut_index} [{cut.cut_kind}] on '{cut.generating_var}': {cut.orig_vars_inequality_str}")
            if n.lp_sol is not None and n.lp_val is not None:
                sol_str = ", ".join(f"{k}={format_approx(v)}" for k, v in n.lp_sol.items())
                print(f"{indent}  Node LP Optimum: z = {format_approx(n.lp_val)} | ({sol_str})")
            else:
                print(f"{indent}  Node LP Status: {n.lp_status}")
            print(f"{indent}  Outcome: {n.outcome} -- {n.outcome_reason}")
            inc_str = format_approx(n.incumbent_at_step) if n.incumbent_at_step is not None else "None"
            print(f"{indent}  Global Incumbent: {inc_str}")

        stats = result["stats"]
        print("\n" + "-" * 82)
        print(
            f" Tree Summary: Total Nodes={result['total_nodes']} | Branched={stats['branched']} | "
            f"Integer Leaves={stats['integer_leaves']} | Pruned by Bound={stats['pruned_by_bound']} | "
            f"Infeasible={stats['infeasible']} | Total Cuts={result['total_cuts']}"
        )

    # 5. Final Results
    print("\n" + "=" * 82)
    print(" [5] FINAL OPTIMAL RESULTS".center(82))
    print("=" * 82)
    print(f" Status: {result['status']}")
    if result["final_opt_val"] is not None:
        print(f" Optimal Objective ({lp.obj_type.upper()} z*): {format_fraction(result['final_opt_val'])} (~{float(result['final_opt_val']):.4g})")
        if result["relax_opt_val"] is not None:
            gap = abs(result["relax_opt_val"] - result["final_opt_val"])
            print(f" Integrality Gap |z_LP - z*|: {format_fraction(gap)} (~{float(gap):.4g})")
        print(" Optimal Solution Vector:")
        for var, vt in zip(lp.var_names, lp.var_types):
            val = result["final_solution"][var]
            print(f"   {var:<6} = {format_approx(val):<16} [{vt}]")
    else:
        print(" No integer feasible solution found.")

    # 6. SciPy MILP Verification
    if scipy_res:
        print("\n" + "=" * 82)
        print(" [6] SCIPY MILP CROSS-VERIFICATION".center(82))
        print("=" * 82)
        if scipy_res.get("installed") and scipy_res.get("success"):
            print(f" SciPy MILP Solver: scipy.optimize.milp (HiGHS Branch-and-Cut)")
            if scipy_res.get("scipy_fun") is not None:
                print(f" SciPy Optimal Objective: {scipy_res['scipy_fun']}")
            if scipy_res.get("scipy_x") is not None:
                print(f" SciPy Optimal Solution: {scipy_res['scipy_x']}")
            flag = "[PASSED]" if scipy_res.get("match") else "[FAILED]"
            print(f"\n {flag} {scipy_res['message']}")
        else:
            print(f" {scipy_res.get('message', scipy_res.get('error', 'SciPy verification unavailable.'))}")
    print("=" * 82 + "\n")


# ---------------------------------------------------------------------
# CLI ENTRY POINT
# ---------------------------------------------------------------------

def get_builtin_example(name: str) -> str:
    if name == "gomory2d":
        return """max 2*x1 + 3*x2
s.t.
2*x1 + 5*x2 <= 16
6*x1 + 2*x2 <= 27
x1, x2 >= 0"""
    elif name == "tum":
        return """max x1 + x2
s.t.
x1 <= 3
x2 <= 4
x1 + x2 <= 5
x1, x2 >= 0"""
    elif name == "bnb":
        return """max 5*x1 + 8*x2
s.t.
x1 + x2 <= 6
5*x1 + 9*x2 <= 45
x1, x2 >= 0"""
    elif name == "milp":
        return """max 3*x1 + 4*x2
s.t.
2*x1 + 3*x2 <= 13
2*x1 + x2 <= 9
x1 >= 0, integer
x2 >= 0, continuous"""
    else:
        raise ValueError(f"Unknown example '{name}'. Options: 'gomory2d', 'tum', 'bnb', 'milp'")


def main():
    parser = argparse.ArgumentParser(
        description="Integer & Mixed-Integer LP Solver: TUM, Chvatal-Gomory / GMI Cuts, Branch & Bound, Branch & Cut, and 2D Plotting",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--lp", type=str, help="LP formulation string or path to .lp file")
    group.add_argument("--example", choices=["gomory2d", "tum", "bnb", "milp"], help="Run built-in benchmark example")

    parser.add_argument(
        "--method",
        choices=["cutting-plane", "cg", "branch-and-bound", "bnb", "branch-and-cut", "bnc"],
        default="cutting-plane",
        help="Solving method: cutting-plane (default), branch-and-bound (bnb), or branch-and-cut (bnc)"
    )
    parser.add_argument(
        "--strategy",
        choices=["best-bound", "dfs", "bfs"],
        default="best-bound",
        help="Node selection strategy for Branch & Bound / Branch & Cut (default: best-bound)"
    )
    parser.add_argument(
        "--branch-rule",
        choices=["most-fractional", "first-fractional"],
        default="most-fractional",
        help="Branching variable selection rule (default: most-fractional)"
    )
    parser.add_argument(
        "--cut-type",
        choices=["auto", "cg", "gmi"],
        default="auto",
        help="Cut generator type: auto (CG when valid, GMI for MILP/fractional), cg, or gmi"
    )
    parser.add_argument("--max-cuts", type=int, default=15, help="Maximum cuts in cutting-plane mode (default: 15)")
    parser.add_argument("--cuts-per-node", type=int, default=2, help="Max cut rounds per tree node in branch-and-cut mode (default: 2)")
    parser.add_argument("--cut-depth", type=int, default=5, help="Max tree depth to generate cuts in branch-and-cut mode (default: 5)")
    parser.add_argument("--max-nodes", type=int, default=200, help="Maximum nodes to explore in branch-and-bound / branch-and-cut (default: 200)")
    parser.add_argument("--visualize-2d", type=str, help="Output path for standalone 2D interactive HTML visualizer")
    parser.add_argument("--output-md", type=str, help="Output path for Markdown summary report")
    parser.add_argument("--format", choices=["text", "json"], default="text", help="Terminal output format (default: text)")
    parser.add_argument("--verify-scipy", action="store_true", help="Cross-verify with scipy.optimize.milp if installed")

    args = parser.parse_args()

    # Load problem text
    if args.example:
        lp_text = get_builtin_example(args.example)
    else:
        if os.path.exists(args.lp):
            with open(args.lp, "r", encoding="utf-8") as f:
                lp_text = f.read()
        else:
            lp_text = args.lp

    lp = parse_lp_input(lp_text)

    # Normalize method aliases
    method_map = {
        "cutting-plane": "cutting-plane",
        "cg": "cutting-plane",
        "branch-and-bound": "branch-and-bound",
        "bnb": "branch-and-bound",
        "branch-and-cut": "branch-and-cut",
        "bnc": "branch-and-cut"
    }
    method = method_map[args.method]

    # 1. Total Unimodularity Analysis
    matrix_A = [coeffs for coeffs, _, _ in lp.constraints]
    vector_b = [rhs for _, _, rhs in lp.constraints]
    tum_report = analyze_tum(matrix_A, vector_b, lp.var_names)

    # 2. Solve with requested algorithm
    if method == "cutting-plane":
        result = solve_chvatal_gomory(lp, max_cuts=args.max_cuts, cut_type=args.cut_type)
    elif method == "branch-and-bound":
        result = solve_branch_and_bound(
            lp, strategy=args.strategy, branch_rule=args.branch_rule, max_nodes=args.max_nodes
        )
    else:
        result = solve_branch_and_cut(
            lp,
            strategy=args.strategy,
            branch_rule=args.branch_rule,
            cut_type=args.cut_type,
            cuts_per_node=args.cuts_per_node,
            cut_depth=args.cut_depth,
            max_nodes=args.max_nodes
        )

    # 3. Optional SciPy Verification
    scipy_res = None
    if args.verify_scipy:
        scipy_res = verify_with_scipy_milp(lp, result)

    # 4. Terminal / JSON Output
    if args.format == "json":
        json_out: Dict[str, Any] = {
            "method": result["method"],
            "status": result["status"],
            "obj_type": lp.obj_type,
            "is_tum": tum_report["is_tum"],
            "is_ideal_polyhedron": tum_report["is_ideal"],
            "relax_opt_val": str(result["relax_opt_val"]) if result["relax_opt_val"] is not None else None,
            "relax_solution": {k: str(v) for k, v in result["relax_solution"].items()},
            "total_cuts": result["total_cuts"],
            "cuts": [
                {
                    "cut_index": c.cut_index,
                    "kind": c.cut_kind,
                    "node_id": c.node_id,
                    "generating_var": c.generating_var,
                    "tableau_ineq": c.tableau_inequality_str,
                    "orig_ineq": c.orig_vars_inequality_str
                }
                for c in result["cuts_history"]
            ],
            "final_opt_val": str(result["final_opt_val"]) if result["final_opt_val"] is not None else None,
            "final_solution": {k: str(v) for k, v in result["final_solution"].items()},
        }
        if "nodes" in result:
            json_out["strategy"] = result["strategy"]
            json_out["branch_rule"] = result["branch_rule"]
            json_out["total_nodes"] = result["total_nodes"]
            json_out["stats"] = result["stats"]
            json_out["nodes"] = [
                {
                    "node_id": n.node_id,
                    "parent_id": n.parent_id,
                    "depth": n.depth,
                    "branch_desc": n.branch_desc,
                    "lp_status": n.lp_status,
                    "lp_val": str(n.lp_val) if n.lp_val is not None else None,
                    "outcome": n.outcome,
                    "cuts_generated": len(n.cuts_generated)
                }
                for n in result["nodes"]
            ]
        if scipy_res is not None:
            json_out["scipy_verification"] = scipy_res
        print(json.dumps(json_out, indent=2))
    else:
        print_terminal_report(lp, tum_report, result, scipy_res)

    # 5. 2D Visualization
    if args.visualize_2d:
        if lp.is_2d():
            html_path = generate_2d_html_visualizer(lp, result, tum_report, args.visualize_2d)
            if args.format != "json":
                print(f"[Visualizer Exported] Interactive 2D visualizer written to: {html_path}")
        else:
            if args.format != "json":
                print(f"[Visualizer Warning] Problem has {lp.num_vars} variables; 2D visualizer requires exactly 2 variables.")

    # 6. Markdown Report
    if args.output_md:
        md_report = generate_markdown_report(lp, tum_report, result, scipy_res)
        with open(args.output_md, "w", encoding="utf-8") as f:
            f.write(md_report)
        if args.format != "json":
            print(f"[Report Exported] Markdown report written to: {args.output_md}")


if __name__ == "__main__":
    main()
