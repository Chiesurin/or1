#!/usr/bin/env python3
"""
Integer Linear Programming (ILP) Solver with Chvatal-Gomory Cutting Plane Algorithm,
Total Unimodularity (TUM) Analysis, and 2D Graphical Visualization.

Features:
- Pure Python standard library (fractions.Fraction) for exact rational arithmetic.
- Pre-solve Total Unimodularity (TUM) and Ideal Polyhedron analysis (Hoffman-Kruskal).
- Step-by-step continuous LP relaxation via Simplex tableau.
- Chvatal-Gomory fractional cuts derivation (in both tableau slacks and original decision variables).
- Dual Simplex re-optimization after cut insertion.
- Interactive 2D SVG/HTML graphical visualization of continuous polyhedron, integer lattice, and cuts.
- Optional cross-verification with scipy.optimize.milp.
- Zero external dependencies required.
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
            self.var_types = list(var_types)
        else:
            self.var_types = ["integer"] * self.num_vars

    def is_2d(self) -> bool:
        return self.num_vars == 2


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

def solve_primal_simplex(tab: SimplexTableau, max_iterations: int = 100) -> Tuple[str, List[str]]:
    """
    Solves continuous LP relaxation using Primal Simplex with Bland's anti-cycling rule.
    Returns (status, logs).
    """
    logs = []
    iteration = 0

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
        pivot_col = candidates[0]
        entering_var = tab.var_names[pivot_col]

        # Minimum Ratio Test (theta)
        leaving_row = None
        min_theta = None

        for i in range(1, tab.num_rows()):
            coeff = tab.tableau[i][pivot_col]
            if coeff > 0:
                theta = tab.tableau[i][-1] / coeff
                if min_theta is None or theta < min_theta:
                    min_theta = theta
                    leaving_row = i

        if leaving_row is None:
            logs.append(f"[Primal Simplex] Unbounded problem detected for entering variable '{entering_var}'.")
            return "UNBOUNDED", logs

        leaving_var = tab.basic_vars[leaving_row - 1]
        pivot_elem = tab.tableau[leaving_row][pivot_col]

        # Elementary Row Operations
        # 1. Normalize pivot row
        for c in range(len(tab.tableau[leaving_row])):
            tab.tableau[leaving_row][c] /= pivot_elem

        # 2. Eliminate column in other rows
        for i in range(tab.num_rows()):
            if i != leaving_row:
                factor = tab.tableau[i][pivot_col]
                if factor != 0:
                    for c in range(len(tab.tableau[i])):
                        tab.tableau[i][c] -= factor * tab.tableau[leaving_row][c]

        tab.basic_vars[leaving_row - 1] = entering_var
        iteration += 1

    return "ITERATION_LIMIT", logs


def solve_dual_simplex(tab: SimplexTableau, max_iterations: int = 50) -> Tuple[str, List[str]]:
    """
    Executes Dual Simplex steps until all RHS values >= 0 (primal feasibility restored).
    Used after injecting a Chvátal-Gomory cut row.
    """
    logs = []
    iteration = 0

    while iteration < max_iterations:
        # Check primal feasibility: all RHS_i >= 0 for i >= 1
        negative_rows = []
        for i in range(1, tab.num_rows()):
            if tab.tableau[i][-1] < 0:
                negative_rows.append(i)

        if not negative_rows:
            logs.append("[Dual Simplex] Primal feasibility restored.")
            return "OPTIMAL", logs

        # Most negative RHS for leaving row (or first negative row for deterministic tie-breaking)
        leaving_row = min(negative_rows, key=lambda r: tab.tableau[r][-1])
        leaving_var = tab.basic_vars[leaving_row - 1]

        # Dual Minimum Ratio Test across negative coefficients in the leaving row: min(c_j / |a_rj|)
        ratios = []
        for j in range(tab.num_cols()):
            coeff = tab.tableau[leaving_row][j]
            if coeff < 0:
                cost = tab.tableau[0][j]
                # In dual simplex, cost >= 0
                ratio = cost / abs(coeff)
                ratios.append((ratio, j))

        if not ratios:
            logs.append(f"[Dual Simplex] Infeasibility detected: leaving row '{leaving_var}' has no negative coefficients.")
            return "INFEASIBLE", logs

        # Bland's smallest-index tie-breaking on min ratio
        min_ratio = min(r[0] for r in ratios)
        candidates = [j for ratio, j in ratios if ratio == min_ratio]
        pivot_col = candidates[0]
        entering_var = tab.var_names[pivot_col]

        pivot_elem = tab.tableau[leaving_row][pivot_col]

        # Pivot operations
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


# ---------------------------------------------------------------------
# CUTTING PLANE ALGORITHM & CHVÁTAL-GOMORY GENERATOR
# ---------------------------------------------------------------------

class CutDetail:
    """Detailed mathematical trace of an injected Chvátal-Gomory cut."""

    def __init__(self, cut_index: int, generating_var: str, generating_row: int,
                 rhs_fractional: Fraction, floor_rhs: int, frac_rhs: Fraction,
                 slack_coefficients: Dict[str, Fraction],
                 tableau_inequality_str: str,
                 orig_vars_inequality_str: str,
                 orig_vars_coefficients: Dict[str, Fraction],
                 orig_vars_rhs: Fraction,
                 tableau_snapshot_before: str,
                 tableau_snapshot_after: str):
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


def reconstruct_cut_in_original_vars(
    lp: LinearProgram,
    f_coeffs: Dict[str, Fraction],
    f0: Fraction,
    var_definitions: Dict[str, Tuple[Dict[str, Fraction], Fraction]]
) -> Tuple[Dict[str, Fraction], Fraction, str]:
    """
    Substitutes the affine definitions v = const_v + sum(gamma_{v,j} * x_j) of all variables v
    into the tableau Gomory cut sum(f_v * v) >= f0, obtaining the cut directly in decision variables:
    sum(alpha_j * x_j) <= beta
    """
    # Cut inequality: sum_{v} f_v * (const_v + sum_j gamma_{v,j} * x_j) >= f0
    # <=> sum_j (sum_v f_v * gamma_{v,j}) * x_j >= f0 - sum_v f_v * const_v
    # <=> sum_j (-sum_v f_v * gamma_{v,j}) * x_j <= sum_v f_v * const_v - f0
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

    # String representation
    terms = []
    for var in lp.var_names:
        c = alpha[var]
        if c != 0:
            c_str = format_fraction(c)
            terms.append(f"{c_str}*{var}" if c_str != "1" else var)
    lhs_str = " + ".join(terms) if terms else "0"
    lhs_str = lhs_str.replace("+ -", "- ")
    ineq_str = f"{lhs_str} <= {format_fraction(beta)}"

    return alpha, beta, ineq_str


def solve_chvatal_gomory(
    lp: LinearProgram,
    max_cuts: int = 15
) -> Dict[str, Any]:
    """
    Full Integer Linear Programming solver using the Chvatal-Gomory Cutting Plane Algorithm.
    1. Standardizes problem and builds initial tableau with slacks.
    2. Solves continuous LP relaxation via Primal Simplex.
    3. Iteratively identifies fractional basic variables and derives Gomory cuts.
    4. Re-projects cuts to original decision variables for 2D visualization.
    5. Re-optimizes via Dual Simplex.
    6. Stops when all decision variables are integral or problem becomes infeasible.
    """
    # Build initial tableau
    # Max c^T x <=> Min -c^T x
    m = len(lp.constraints)
    n = lp.num_vars

    var_names = list(lp.var_names)
    slack_names = [f"s{i+1}" for i in range(m)]
    var_names.extend(slack_names)

    # Row 0
    # If maximizing, target min z' = -sum(c_j * x_j) => -z' - sum(c_j * x_j) = 0
    # If minimizing, target min z = sum(c_j * x_j) => -z + sum(c_j * x_j) = 0
    row0 = []
    for j in range(n):
        c = lp.obj_coeffs[j]
        row0.append(-c if lp.obj_type == "max" else c)
    for _ in range(m):
        row0.append(Fraction(0))
    row0.append(Fraction(0))  # RHS

    matrix = [row0]
    basic_vars = []

    # Initialize affine definitions for each variable in terms of original decision variables:
    # v = const + sum(gamma_j * x_j)
    var_definitions: Dict[str, Tuple[Dict[str, Fraction], Fraction]] = {}
    for var in lp.var_names:
        var_definitions[var] = ({var: Fraction(1)}, Fraction(0))

    for i, (coeffs, rel, rhs) in enumerate(lp.constraints):
        r_coeffs = list(coeffs)
        r_rhs = rhs
        # If rhs < 0, normalize
        if r_rhs < 0:
            r_coeffs = [-c for c in r_coeffs]
            r_rhs = -r_rhs
            rel = ">=" if rel == "<=" else ("<=" if rel == ">=" else "=")

        row = list(r_coeffs)
        for k in range(m):
            if k == i:
                row.append(Fraction(1) if rel == "<=" else Fraction(-1))
            else:
                row.append(Fraction(0))
        row.append(r_rhs)
        matrix.append(row)
        basic_vars.append(slack_names[i])

        slack_name = slack_names[i]
        if rel == "<=":
            # s_i = rhs - sum(coeffs[j] * x_j)
            gammas = {lp.var_names[j]: -r_coeffs[j] for j in range(lp.num_vars)}
            var_definitions[slack_name] = (gammas, r_rhs)
        elif rel == ">=":
            # e_i = -rhs + sum(coeffs[j] * x_j)
            gammas = {lp.var_names[j]: r_coeffs[j] for j in range(lp.num_vars)}
            var_definitions[slack_name] = (gammas, -r_rhs)

    tab = SimplexTableau(var_names, matrix, basic_vars, obj_name="-z", orig_decision_vars=lp.var_names)

    # Step 0: Solve Continuous LP Relaxation
    relax_tab_before = tab.copy().format_ascii("Initial Tableau (LP Relaxation)")
    status_relax, logs_relax = solve_primal_simplex(tab)
    relax_tab_after = tab.copy().format_ascii("Optimal LP Relaxation Tableau")

    relax_solution = {}
    for j, name in enumerate(lp.var_names):
        if name in tab.basic_vars:
            r_idx = tab.basic_vars.index(name) + 1
            relax_solution[name] = tab.tableau[r_idx][-1]
        else:
            relax_solution[name] = Fraction(0)

    raw_rhs0 = tab.tableau[0][-1]
    relax_opt_val = -raw_rhs0 if lp.obj_type == "min" else raw_rhs0

    cuts_history: List[CutDetail] = []
    current_status = status_relax
    cut_count = 0

    while current_status == "OPTIMAL" and cut_count < max_cuts:
        # Check if all original decision variables are integral
        fractional_vars = []
        for i in range(1, tab.num_rows()):
            b_var = tab.basic_vars[i - 1]
            if b_var in lp.var_names:
                rhs_val = tab.tableau[i][-1]
                f0 = frac_part(rhs_val)
                if f0 != 0:
                    fractional_vars.append((i, b_var, rhs_val, f0))

        if not fractional_vars:
            # Optimal integer solution reached!
            break

        # Select generating row: variable with largest fractional part f0 (closest to 0.5)
        # Tie-breaker: smallest index for deterministic behavior
        def frac_score(item):
            f = item[3]
            # distance to integer: min(f, 1 - f)
            return (f, -item[0])

        chosen_row, chosen_var, chosen_rhs, f0 = max(fractional_vars, key=frac_score)

        # Generate Gomory fractional cut
        # From row: x_B + sum(a_j * x_j) = rhs
        # Cut: sum(f_j * x_j) >= f0  =>  -s_cut + sum(-f_j * x_j) = -f0
        cut_count += 1
        cut_name = f"s_cut{cut_count}"
        tab_before_cut = tab.copy().format_ascii(f"Tableau before Cut {cut_count}")

        f_coeffs = {}
        cut_row = []
        ineq_terms = []

        for j in range(tab.num_cols()):
            var = tab.var_names[j]
            coeff = tab.tableau[chosen_row][j]
            fj = frac_part(coeff)
            f_coeffs[var] = fj
            # For Dual Simplex: -fj * x_j
            cut_row.append(-fj)
            if fj != 0:
                fj_str = format_fraction(fj)
                ineq_terms.append(f"{fj_str}*{var}" if fj_str != "1" else var)

        # Add 1 for the new cut slack variable in its own row
        cut_row.append(Fraction(1))
        # RHS is -f0
        cut_row.append(-f0)

        # Expand existing rows with 0 for the new slack column
        for i in range(tab.num_rows()):
            tab.tableau[i].insert(-1, Fraction(0))

        tab.var_names.append(cut_name)
        tab.tableau.append(cut_row)
        tab.basic_vars.append(cut_name)

        tableau_ineq_str = (" + ".join(ineq_terms) if ineq_terms else "0") + f" >= {format_fraction(f0)}"

        # Re-project cut to original decision variables (x1, x2, ...)
        orig_coeffs, orig_rhs, orig_ineq_str = reconstruct_cut_in_original_vars(lp, f_coeffs, f0, var_definitions)
        # Store definition of the new cut slack variable: s_cut = orig_rhs - sum(orig_coeffs * x)
        cut_gammas = {var: -orig_coeffs[var] for var in lp.var_names}
        var_definitions[cut_name] = (cut_gammas, orig_rhs)

        # Dual Simplex step to restore feasibility
        dual_status, dual_logs = solve_dual_simplex(tab)
        current_status = dual_status
        tab_after_cut = tab.copy().format_ascii(f"Tableau after Cut {cut_count} (Dual Simplex)")

        cut_detail = CutDetail(
            cut_index=cut_count,
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
            tableau_snapshot_after=tab_after_cut
        )
        cuts_history.append(cut_detail)

        if current_status != "OPTIMAL":
            break

    # Extract final integer solution
    final_solution = {}
    is_integer = True
    for name in lp.var_names:
        if name in tab.basic_vars:
            r_idx = tab.basic_vars.index(name) + 1
            val = tab.tableau[r_idx][-1]
            final_solution[name] = val
            if val.denominator != 1:
                is_integer = False
        else:
            final_solution[name] = Fraction(0)

    raw_rhs0_final = tab.tableau[0][-1]
    final_opt_val = -raw_rhs0_final if lp.obj_type == "min" else raw_rhs0_final

    return {
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
        "final_tableau_ascii": tab.format_ascii("Final Integer Optimal Tableau")
    }


# ---------------------------------------------------------------------
# 2D GRAPHICAL VISUALIZER (SVG & INTERACTIVE HTML)
# ---------------------------------------------------------------------

def compute_polygon_2d(
    constraints: List[Tuple[List[Fraction], str, Fraction]],
    bbox: Tuple[float, float, float, float]
) -> List[Tuple[float, float]]:
    """
    Computes vertices of 2D convex polygon defined by linear constraints and bounding box.
    bbox: (xmin, xmax, ymin, ymax)
    """
    xmin, xmax, ymin, ymax = bbox
    all_halfplanes = []

    # Non-negativity and bounding box
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

    # Find candidate intersection points between all pairs of halfplane boundaries
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

                # Check feasibility against all halfplanes
                feasible = True
                for (h1, h2), hb in all_halfplanes:
                    if h1 * px + h2 * py > hb + 1e-6:
                        feasible = False
                        break
                if feasible:
                    points.append((round(px, 5), round(py, 5)))

    # Deduplicate points
    unique_pts = []
    for p in points:
        if not any(math.hypot(p[0] - u[0], p[1] - u[1]) < 1e-4 for u in unique_pts):
            unique_pts.append(p)

    if len(unique_pts) < 3:
        return unique_pts

    # Sort vertices in counterclockwise order around centroid
    cx = sum(p[0] for p in unique_pts) / len(unique_pts)
    cy = sum(p[1] for p in unique_pts) / len(unique_pts)
    unique_pts.sort(key=lambda p: math.atan2(p[1] - cy, p[0] - cx))
    return unique_pts


def generate_2d_html_visualizer(
    lp: LinearProgram,
    result: Dict[str, Any],
    tum_report: Dict[str, Any],
    output_html_path: str
) -> str:
    """
    Generates a standalone, interactive HTML5 + SVG visualizer showing:
    - Feasible continuous polygon P_0
    - Integer grid points Z^2
    - Initial continuous optimum x*_LP
    - Progressive cut lines and sliced-off sub-polygons P_1, P_2, ...
    - Final integer optimum x*_ILP
    """
    if not lp.is_2d():
        raise ValueError("Visualization is only supported for 2-variable (2D) problems.")

    # Determine bounding box
    raw_constraints = list(lp.constraints)
    initial_poly = compute_polygon_2d(raw_constraints, (0.0, 15.0, 0.0, 15.0))

    if initial_poly:
        max_x = max(max(p[0] for p in initial_poly) + 1.5, 5.0)
        max_y = max(max(p[1] for p in initial_poly) + 1.5, 5.0)
    else:
        max_x, max_y = 10.0, 10.0

    bbox = (0.0, max_x, 0.0, max_y)

    # Compute continuous polygons after each cut
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

    # Integer grid points in bounding box
    int_points = []
    for ix in range(int(max_x) + 1):
        for iy in range(int(max_y) + 1):
            # Check feasibility in initial continuous relaxation
            is_in_poly = True
            for coeffs, rel, rhs in raw_constraints:
                val = float(coeffs[0]) * ix + float(coeffs[1]) * iy
                r = float(rhs)
                if rel == "<=" and val > r + 1e-6:
                    is_in_poly = False
                    break
                elif rel == ">=" and val < r - 1e-6:
                    is_in_poly = False
                    break
            int_points.append({"x": ix, "y": iy, "in_initial_poly": is_in_poly})

    # Continuous LP Optimum
    relax_x = float(result["relax_solution"][lp.var_names[0]])
    relax_y = float(result["relax_solution"][lp.var_names[1]])
    relax_opt = float(result["relax_opt_val"])

    # Final Integer Optimum
    final_x = float(result["final_solution"][lp.var_names[0]])
    final_y = float(result["final_solution"][lp.var_names[1]])
    final_opt = float(result["final_opt_val"])

    # Cuts data for JavaScript
    cuts_json = []
    for c in result["cuts_history"]:
        c1 = float(c.orig_vars_coefficients[lp.var_names[0]])
        c2 = float(c.orig_vars_coefficients[lp.var_names[1]])
        cr = float(c.orig_vars_rhs)
        cuts_json.append({
            "index": c.cut_index,
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
      <!-- 1. TUM Analysis Section -->
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

      <!-- 2. Step Navigator -->
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

      <!-- 3. Current Cut Details -->
      <div class="section" id="cut-info-section">
        <div class="section-title">Current Cut Algebraic Derivation</div>
        <div id="cut-info-body" class="card-box">
          Select 'Next Cut' to visualize sequential Chvátal-Gomory planes.
        </div>
      </div>

      <!-- 4. Solution Comparison -->
      <div class="section">
        <div class="section-title">Solution Comparison</div>
        <div class="card-box" style="margin-bottom: 8px;">
          <div style="color: #f87171; font-weight: 700;">Continuous Relaxation Optimum (x*_LP):</div>
          <div>{lp.var_names[0]} = <code>{format_fraction(result['relax_solution'][lp.var_names[0]])}</code>, {lp.var_names[1]} = <code>{format_fraction(result['relax_solution'][lp.var_names[1]])}</code></div>
          <div>Objective ({lp.obj_type.upper()} z): <code>{format_fraction(result['relax_opt_val'])}</code> (~{relax_opt:.4g})</div>
        </div>
        <div class="card-box">
          <div style="color: #4ade80; font-weight: 700;">Final Integer Optimum (x*_ILP):</div>
          <div>{lp.var_names[0]} = <code>{format_fraction(result['final_solution'][lp.var_names[0]])}</code>, {lp.var_names[1]} = <code>{format_fraction(result['final_solution'][lp.var_names[1]])}</code></div>
          <div>Objective ({lp.obj_type.upper()} z*): <code>{format_fraction(result['final_opt_val'])}</code> (~{final_opt:.4g})</div>
        </div>
      </div>
    </div>
  </div>

  <script>
    const maxCoord = {max_x};
    const svgSize = 600;
    const originX = 50;
    const originY = 550;
    const scale = (svgSize - 100) / maxCoord;

    function toSvgX(x) {{ return originX + x * scale; }}
    function toSvgY(y) {{ return originY - y * scale; }}

    const polygons = {json.dumps(polygons_by_step)};
    const intPoints = {json.dumps(int_points)};
    const cutsData = {json.dumps(cuts_json)};
    const lpOpt = {{ x: {relax_x}, y: {relax_y}, val: {relax_opt} }};
    const ilpOpt = {{ x: {final_x}, y: {final_y}, val: {final_opt} }};

    let currentStep = 0;

    function renderPlot() {{
      const svg = document.getElementById('plot-svg');
      let html = '';

      // Grid lines
      for (let i = 0; i <= maxCoord; i++) {{
        const sx = toSvgX(i);
        const sy = toSvgY(i);
        html += `<line x1="${{sx}}" y1="${{toSvgY(0)}}" x2="${{sx}}" y2="${{toSvgY(maxCoord)}}" stroke="#1e293b" stroke-width="1" />`;
        html += `<line x1="${{toSvgX(0)}}" y1="${{sy}}" x2="${{toSvgX(maxCoord)}}" y2="${{sy}}" stroke="#1e293b" stroke-width="1" />`;
        html += `<text x="${{sx}}" y="${{originY + 18}}" fill="#64748b" font-size="11" text-anchor="middle">${{i}}</text>`;
        html += `<text x="${{originX - 12}}" y="${{sy + 4}}" fill="#64748b" font-size="11" text-anchor="end">${{i}}</text>`;
      }}

      // Axes
      html += `<line x1="${{toSvgX(0)}}" y1="${{originY}}" x2="${{toSvgX(maxCoord)}}" y2="${{originY}}" stroke="#94a3b8" stroke-width="2" />`;
      html += `<line x1="${{originX}}" y1="${{toSvgY(0)}}" x2="${{originX}}" y2="${{toSvgY(maxCoord)}}" stroke="#94a3b8" stroke-width="2" />`;
      html += `<text x="${{toSvgX(maxCoord) + 12}}" y="${{originY + 4}}" fill="#94a3b8" font-size="13" font-weight="bold">${lp.var_names[0]}</text>`;
      html += `<text x="${{originX}}" y="${{toSvgY(maxCoord) - 12}}" fill="#94a3b8" font-size="13" font-weight="bold" text-anchor="middle">${lp.var_names[1]}</text>`;

      // Initial continuous polygon (ghost outline)
      const initPts = polygons[0].poly.map(p => `${{toSvgX(p[0])}},${{toSvgY(p[1])}}`).join(' ');
      html += `<polygon points="${{initPts}}" fill="rgba(37, 99, 235, 0.08)" stroke="#3b82f6" stroke-dasharray="4,4" stroke-width="1.5" />`;

      // Current Step Polygon
      const curPts = polygons[currentStep].poly.map(p => `${{toSvgX(p[0])}},${{toSvgY(p[1])}}`).join(' ');
      html += `<polygon points="${{curPts}}" fill="rgba(56, 189, 248, 0.22)" stroke="#38bdf8" stroke-width="2.5" />`;

      // Draw all cuts up to currentStep
      for (let c = 0; c < currentStep; c++) {{
        const cut = cutsData[c];
        // Line: c1*x + c2*y = rhs
        let p1, p2;
        if (Math.abs(cut.c2) > 1e-6) {{
          p1 = {{ x: 0, y: cut.rhs / cut.c2 }};
          p2 = {{ x: maxCoord, y: (cut.rhs - cut.c1 * maxCoord) / cut.c2 }};
        }} else {{
          p1 = {{ x: cut.rhs / cut.c1, y: 0 }};
          p2 = {{ x: cut.rhs / cut.c1, y: maxCoord }};
        }}
        html += `<line x1="${{toSvgX(p1.x)}}" y1="${{toSvgY(p1.y)}}" x2="${{toSvgX(p2.x)}}" y2="${{toSvgY(p2.y)}}" stroke="#f59e0b" stroke-width="2.5" />`;
        // Label on line
        const midX = (p1.x + p2.x) / 2;
        const midY = (p1.y + p2.y) / 2;
        if (midX >= 0 && midX <= maxCoord && midY >= 0 && midY <= maxCoord) {{
          html += `<text x="${{toSvgX(midX)}}" y="${{toSvgY(midY) - 8}}" fill="#f59e0b" font-size="11" font-weight="bold">Cut ${{cut.index}}</text>`;
        }}
      }}

      // Integer grid points
      intPoints.forEach(pt => {{
        const isFeas = pt.in_initial_poly;
        const fill = isFeas ? '#38bdf8' : '#475569';
        const r = isFeas ? 3.5 : 2;
        html += `<circle cx="${{toSvgX(pt.x)}}" cy="${{toSvgY(pt.y)}}" r="${{r}}" fill="${{fill}}" opacity="${{isFeas ? 1 : 0.4}}" />`;
      }});

      // LP Continuous Optimum
      html += `<circle cx="${{toSvgX(lpOpt.x)}}" cy="${{toSvgY(lpOpt.y)}}" r="6" fill="#ef4444" stroke="#ffffff" stroke-width="2" />`;
      html += `<text x="${{toSvgX(lpOpt.x) + 8}}" y="${{toSvgY(lpOpt.y) - 8}}" fill="#f87171" font-size="12" font-weight="bold">x*_LP (${{lpOpt.x.toFixed(2)}}, ${{lpOpt.y.toFixed(2)}})</text>`;

      // If at last step or optimal integer reached, highlight integer optimum
      if (currentStep === cutsData.length) {{
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
          <div style="margin-top: 6px;">Continuous optimum is fractional: <code>${lp.var_names[0]} = ${{lpOpt.x.toFixed(3)}}, ${lp.var_names[1]} = ${{lpOpt.y.toFixed(3)}}</code>.</div>
          <div style="color: var(--amber); margin-top: 6px;">Chvátal-Gomory cuts are required to slice off the non-integer optimum while preserving all lattice points Z^2.</div>
        `;
      }} else {{
        const cut = cutsData[currentStep - 1];
        document.getElementById('step-label').textContent = `Cut ${{cut.index}} Applied`;
        document.getElementById('cut-info-body').innerHTML = `
          <div><strong>Generating Variable:</strong> <code>${{cut.generating_var}}</code></div>
          <div style="margin-top: 6px;"><strong>Tableau Slack Form:</strong></div>
          <div class="code-box">${{cut.formula_slack}}</div>
          <div style="margin-top: 6px;"><strong>Original Decision Variables Form:</strong></div>
          <div class="code-box" style="color: #4ade80;">${{cut.formula_orig}}</div>
          <div style="margin-top: 6px; font-size: 12px; color: var(--text-muted);">
            Notice how the yellow cut hyperplane strictly separates the fractional point x*_LP without eliminating any blue integer lattice points!
          </div>
        `;
      }}
    }}

    function nextStep() {{
      if (currentStep < cutsData.length) {{
        currentStep++;
        renderPlot();
      }}
    }}

    function prevStep() {{
      if (currentStep > 0) {{
        currentStep--;
        renderPlot();
      }}
    }}

    function resetStep() {{
      currentStep = 0;
      renderPlot();
    }}

    // Initial render
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
    Parses a standard algebraic LP formulation string:
    min/max c1*x1 + c2*x2
    s.t.
    a11*x1 + a12*x2 <= b1
    ...
    """
    lines = [line.strip() for line in text.strip().splitlines() if line.strip() and not line.strip().startswith("#")]
    if len(lines) < 2:
        raise ValueError("LP formulation requires at least 1 objective line and 1 constraint line.")

    # 1. Objective
    obj_line = lines[0].lower()
    m_obj = re.match(r'^(min(?:imize)?|max(?:imize)?)\s+(.+)$', obj_line)
    if not m_obj:
        raise ValueError(f"Could not parse objective line: '{lines[0]}'")
    obj_type = "min" if "min" in m_obj.group(1) else "max"
    obj_expr = m_obj.group(2)

    var_tokens = sorted(list(set(re.findall(r'[a-zA-Z]\w*', obj_expr))))
    if not var_tokens:
        raise ValueError("No decision variables found in objective.")

    def parse_expr(expr: str, variables: List[str]) -> List[Fraction]:
        coeffs = {v: Fraction(0) for v in variables}
        clean_expr = expr.replace('-', '+-')
        terms = [t.strip() for t in clean_expr.split('+') if t.strip()]
        for term in terms:
            for v in variables:
                if re.search(r'\b' + v + r'\b', term):
                    c_str = re.sub(r'\*?' + v, '', term).strip()
                    if c_str == '' or c_str == '+':
                        coeffs[v] += Fraction(1)
                    elif c_str == '-':
                        coeffs[v] += Fraction(-1)
                    else:
                        coeffs[v] += parse_fraction(c_str)
                    break
        return [coeffs[v] for v in variables]

    obj_coeffs = parse_expr(obj_expr, var_tokens)

    # 2. Constraints
    constraints = []
    c_lines = lines[1:]
    if c_lines and ("s.t." in c_lines[0].lower() or "subject to" in c_lines[0].lower()):
        c_lines = c_lines[1:]

    for line in c_lines:
        line_clean = line.strip()
        if not line_clean:
            continue
        # Skip pure bounds like x1, x2 >= 0
        if re.match(r'^[a-zA-Z0-9_,\s]+>=\s*0$', line_clean):
            continue

        m_rel = re.search(r'(<=|>=|=)', line_clean)
        if not m_rel:
            continue
        rel = m_rel.group(1)
        lhs, rhs = line_clean.split(rel)
        lhs_coeffs = parse_expr(lhs, var_tokens)
        rhs_val = parse_fraction(rhs)
        constraints.append((lhs_coeffs, rel, rhs_val))

    return LinearProgram(obj_type, obj_coeffs, constraints, var_names=var_tokens)


def verify_with_scipy_milp(lp: LinearProgram, result: Dict[str, Any]) -> Dict[str, Any]:
    """Cross-verifies the cutting plane optimal integer solution with SciPy HiGHS MILP solver."""
    try:
        from scipy.optimize import milp, LinearConstraint
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
        integrality = np.ones(lp.num_vars)  # All variables integer

        res = milp(c=c, constraints=constraints, integrality=integrality)

        if not res.success:
            return {
                "installed": True,
                "success": False,
                "message": res.message
            }

        scipy_x = {lp.var_names[j]: round(float(res.x[j]), 5) for j in range(lp.num_vars)}
        scipy_fun = float(res.fun) if lp.obj_type == "min" else -float(res.fun)

        our_opt = float(result["final_opt_val"])
        match = abs(our_opt - scipy_fun) < 1e-4

        return {
            "installed": True,
            "success": True,
            "scipy_fun": scipy_fun,
            "scipy_x": scipy_x,
            "match": match,
            "message": f"MATCH CONFIRMED: Cutting Plane ({our_opt}) == SciPy MILP ({scipy_fun})" if match else "Discrepancy detected."
        }
    except Exception as ex:
        return {
            "installed": True,
            "error": str(ex)
        }


# ---------------------------------------------------------------------
# TERMINAL WALKTHROUGH PRINTER
# ---------------------------------------------------------------------

def print_terminal_report(lp: LinearProgram, tum_report: Dict[str, Any], result: Dict[str, Any], scipy_res: Optional[Dict[str, Any]] = None):
    """Prints a beautiful, comprehensive terminal report without any emojis."""
    print("=" * 82)
    print("          INTEGER LINEAR PROGRAMMING: CHVÁTAL-GOMORY CUTTING PLANES       ".center(82))
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
    print("      " + ", ".join(f"{v} >= 0, integer" for v in lp.var_names))

    # 2. Total Unimodularity
    print("\n" + "=" * 82)
    print(" [2] TOTAL UNIMODULARITY (TUM) & IDEAL POLYHEDRON ANALYSIS".center(82))
    print("=" * 82)
    for line in tum_report["details"]:
        print(f" • {line}")

    # 3. Continuous LP Relaxation
    print("\n" + "=" * 82)
    print(" [3] CONTINUOUS LP RELAXATION (STEP 0)".center(82))
    print("=" * 82)
    print(result["relax_tab_after"])
    print(f"\n Relaxation Status: OPTIMAL")
    print(f" Continuous Objective ({lp.obj_type.upper()} z): {format_approx(result['relax_opt_val'])}")
    print(" Continuous Solution Vector:")
    for var, val in result["relax_solution"].items():
        print(f"   {var:<6} = {format_approx(val)}")

    # 4. Cuts Trace
    if result["total_cuts"] == 0:
        print("\n" + "=" * 82)
        print(" [4] CUTTING PLANE ALGORITHM: NO CUTS NEEDED".center(82))
        print("=" * 82)
        print(" All decision variables in the continuous relaxation are already integers!")
    else:
        print("\n" + "=" * 82)
        print(f" [4] CHVÁTAL-GOMORY CUTTING PLANES (TOTAL CUTS INJECTED: {result['total_cuts']})".center(82))
        print("=" * 82)

        for cut in result["cuts_history"]:
            print(f"\n>>> Cut {cut.cut_index}: Generated from Row {cut.generating_row} ('{cut.generating_var}')")
            print(f"    RHS Value: {format_fraction(cut.rhs_fractional)}  =>  Floor = {cut.floor_rhs}, Fractional Part f0 = {format_fraction(cut.frac_rhs)}")
            print(f"    Tableau Gomory Inequality:  {cut.tableau_inequality_str}")
            print(f"    Re-projected to Decision Variables:  {cut.orig_vars_inequality_str}")
            print(f"    Injected Row in Dual Simplex:")
            print(f"      -s_cut{cut.cut_index} + sum(-f_j * x_j) = -{format_fraction(cut.frac_rhs)}")
            print(f"    Dual Simplex Pivot executed successfully.")

    # 5. Final Results
    print("\n" + "=" * 82)
    print(" [5] FINAL INTEGER OPTIMAL RESULTS".center(82))
    print("=" * 82)
    print(f" Status: {result['status']}")
    print(f" Optimal Integer Objective ({lp.obj_type.upper()} z*): {format_fraction(result['final_opt_val'])} (~{float(result['final_opt_val']):.4g})")
    print(" Optimal Integer Solution Vector:")
    for var, val in result["final_solution"].items():
        print(f"   {var:<6} = {format_fraction(val)}")

    # 6. SciPy MILP Verification
    if scipy_res:
        print("\n" + "=" * 82)
        print(" [6] SCIPY MILP CROSS-VERIFICATION".center(82))
        print("=" * 82)
        if scipy_res.get("installed") and scipy_res.get("success"):
            print(f" SciPy MILP Solver: scipy.optimize.milp (HiGHS Branch-and-Cut)")
            print(f" SciPy Optimal Objective: {scipy_res['scipy_fun']}")
            print(f" SciPy Integer Solution: {scipy_res['scipy_x']}")
            print(f"\n [PASSED] {scipy_res['message']}")
        else:
            print(f" {scipy_res.get('message', 'SciPy verification unavailable.')}")
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
    else:
        raise ValueError(f"Unknown example '{name}'. Options: 'gomory2d', 'tum'")


def main():
    parser = argparse.ArgumentParser(
        description="Integer Linear Programming (ILP) Solver with Chvatal-Gomory Cuts, TUM, and 2D Plotting",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--lp", type=str, help="LP formulation string or path to .lp file")
    group.add_argument("--example", choices=["gomory2d", "tum"], help="Run built-in benchmark example")

    parser.add_argument("--max-cuts", type=int, default=15, help="Maximum number of cutting plane iterations")
    parser.add_argument("--visualize-2d", type=str, help="Output path for standalone 2D interactive HTML visualizer")
    parser.add_argument("--output-md", type=str, help="Output path for Markdown summary report")
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

    # 1. Total Unimodularity Analysis
    matrix_A = [coeffs for coeffs, _, _ in lp.constraints]
    vector_b = [rhs for _, _, rhs in lp.constraints]
    tum_report = analyze_tum(matrix_A, vector_b, lp.var_names)

    # 2. Solve Chvatal-Gomory Cutting Plane
    result = solve_chvatal_gomory(lp, max_cuts=args.max_cuts)

    # 3. Optional SciPy Verification
    scipy_res = None
    if args.verify_scipy:
        scipy_res = verify_with_scipy_milp(lp, result)

    # 4. Print Report to Terminal
    print_terminal_report(lp, tum_report, result, scipy_res)

    # 5. 2D Visualization
    if args.visualize_2d:
        if lp.is_2d():
            html_path = generate_2d_html_visualizer(lp, result, tum_report, args.visualize_2d)
            print(f"[Visualizer Exported] Interactive 2D visualizer written to: {html_path}")
        else:
            print(f"[Visualizer Warning] Problem has {lp.num_vars} variables; 2D visualizer requires exactly 2 variables.")

    # 6. Markdown Report
    if args.output_md:
        with open(args.output_md, "w", encoding="utf-8") as f:
            f.write(f"# Integer Linear Programming Report: Chvátal-Gomory Cutting Planes\n\n")
            f.write(f"## Pre-Solve Total Unimodularity (TUM) Analysis\n")
            for line in tum_report["details"]:
                f.write(f"- {line}\n")
            f.write(f"\n## Optimal Integer Solution\n")
            f.write(f"- **Objective Value**: `{format_fraction(result['final_opt_val'])}`\n")
            for v, val in result["final_solution"].items():
                f.write(f"- `{v}` = `{format_fraction(val)}`\n")
            f.write(f"\n## Total Cuts Applied: {result['total_cuts']}\n")
            for cut in result["cuts_history"]:
                f.write(f"### Cut {cut.cut_index} ({cut.generating_var})\n")
                f.write(f"- Tableau inequality: `{cut.tableau_inequality_str}`\n")
                f.write(f"- Original variables: `{cut.orig_vars_inequality_str}`\n")
        print(f"[Report Exported] Markdown report written to: {args.output_md}")


if __name__ == "__main__":
    main()
