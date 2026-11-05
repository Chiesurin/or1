#!/usr/bin/env python3
"""
Simplex Tableau Solver with Step-by-Step Mathematical Calculations
Part of the Tableau Solver Skill for LLMs & AI Agents.

Features:
- Pure Python Standard Library (zero external dependencies).
- Exact rational arithmetic (fractions.Fraction) eliminating floating-point errors.
- Universal Minimization (min z) with Row 0 identity: -z = RHS_0 => z* = -RHS_0.
- Anti-cycling: Bland's Smallest-Subscript Rule for both entering and leaving variable ties.
- Full Two-Phase Simplex support (automatic or forced via --two-phase).
- Step-by-step calculation trace showing exact elementary row operations (R_i <- R_i - c * R_pivot).
- Flexible input: CSV file (--csv), raw tableau string (--raw), stdin pipe (--stdin), or LP text (--lp).
- Output modes: Terminal text, full Markdown report (--output-md), or structured JSON (--format json).
- Optional verification with scipy.optimize.linprog if installed (--verify-scipy).
"""

import sys
import os
import re
import csv
import json
import argparse
from fractions import Fraction
from typing import List, Tuple, Dict, Optional, Set, Any


def parse_fraction(val_str: str) -> Fraction:
    """Parses a string into an exact Fraction. Supports integers, fractions (a/b), and decimals (. or ,)."""
    s = val_str.strip()
    if not s:
        raise ValueError("Cannot parse empty value.")
    # Support European decimal comma e.g. '0,5'
    if re.match(r'^-?\d+,\d+$', s):
        s = s.replace(',', '.')
    if '/' in s:
        parts = s.split('/')
        if len(parts) != 2:
            raise ValueError(f"Invalid fraction format: '{s}'")
        return Fraction(int(parts[0].strip()), int(parts[1].strip()))
    return Fraction(s)


def format_fraction(f: Fraction) -> str:
    """Formats a Fraction as an integer if denominator is 1, else as 'num/den'."""
    if f.denominator == 1:
        return str(f.numerator)
    return f"{f.numerator}/{f.denominator}"


def format_approx(f: Fraction) -> str:
    """Formats a Fraction with its decimal approximation if not an integer."""
    s = format_fraction(f)
    if f.denominator != 1:
        s += f" (~{float(f):.4g})"
    return s


def get_display_obj_val(raw_val: Fraction, obj_type: str = "min", obj_name: str = "-z") -> Fraction:
    """Returns the true objective value (for -z = raw_val, returns z = -raw_val when minimizing)."""
    if obj_name.startswith("-") or obj_type == "min":
        return -raw_val
    return raw_val


class SimplexTableau:
    """
    Represents a Simplex Tableau.
    Variables are uniquely indexed with 1-based indexing for Bland's Rule:
    x1, x2, ..., xn, s1, s2, ..., e1, ..., y1, ...
    """
    def __init__(self, var_names: List[str], var_indices: Dict[str, int],
                 tableau: List[List[Fraction]], basic_vars: List[str], obj_name: str = "-z",
                 orig_row0: Optional[List[Fraction]] = None):
        self.var_names = list(var_names)
        self.var_indices = dict(var_indices)
        self.tableau = [list(row) for row in tableau]
        self.basic_vars = list(basic_vars)
        self.obj_name = obj_name
        self.orig_row0 = list(orig_row0) if orig_row0 is not None else None

    def num_rows(self) -> int:
        return len(self.tableau)

    def num_cols(self) -> int:
        return len(self.var_names)

    def to_markdown_table(self) -> str:
        """Renders the current tableau as a standard Markdown table."""
        headers = ["Basis", "RHS"] + self.var_names
        header_line = "| " + " | ".join(headers) + " |"
        sep_line = "| " + " | ".join(["---"] * len(headers)) + " |"
        lines = [header_line, sep_line]

        # Row 0
        raw_obj_rhs = self.tableau[0][-1]
        r0_cells = [self.obj_name, format_fraction(raw_obj_rhs)] + [
            format_fraction(self.tableau[0][j]) for j in range(self.num_cols())
        ]
        lines.append("| " + " | ".join(r0_cells) + " |")

        # Constraint rows
        for i in range(1, self.num_rows()):
            b_var = self.basic_vars[i - 1]
            row_cells = [b_var, format_fraction(self.tableau[i][-1])] + [
                format_fraction(self.tableau[i][j]) for j in range(self.num_cols())
            ]
            lines.append("| " + " | ".join(row_cells) + " |")

        return "\n".join(lines)

    def to_ascii_table(self) -> str:
        """Renders the current tableau as an ASCII box table."""
        headers = ["Basis", "RHS"] + self.var_names
        table_cells: List[List[str]] = [headers]

        raw_obj_rhs = self.tableau[0][-1]
        r0_cells = [self.obj_name, format_fraction(raw_obj_rhs)] + [
            format_fraction(self.tableau[0][j]) for j in range(self.num_cols())
        ]
        table_cells.append(r0_cells)

        for i in range(1, self.num_rows()):
            b_var = self.basic_vars[i - 1]
            row_cells = [b_var, format_fraction(self.tableau[i][-1])] + [
                format_fraction(self.tableau[i][j]) for j in range(self.num_cols())
            ]
            table_cells.append(row_cells)

        num_cols = len(headers)
        col_widths = [max(len(table_cells[r][c]) for r in range(len(table_cells))) for c in range(num_cols)]
        col_widths = [max(w, 5) for w in col_widths]

        sep_line = "+" + "+".join("-" * (w + 2) for w in col_widths) + "+"
        header_str = "| " + " | ".join(table_cells[0][c].center(col_widths[c]) for c in range(num_cols)) + " |"

        lines = [sep_line, header_str, sep_line]
        for r in range(1, len(table_cells)):
            row_str = "| " + " | ".join(
                (table_cells[r][c].center(col_widths[c]) if c == 0 else table_cells[r][c].rjust(col_widths[c]))
                for c in range(num_cols)
            ) + " |"
            lines.append(row_str)
            if r == 1:
                lines.append(sep_line)
        lines.append(sep_line)
        return "\n".join(lines)


# ---------------------------------------------------------------------
# PARSERS: CSV AND STRING
# ---------------------------------------------------------------------

def parse_tableau_text(raw_text: str) -> SimplexTableau:
    """Parses raw text (from CSV, multiline prompt, or stdin) into a SimplexTableau."""
    lines = [line.strip() for line in raw_text.splitlines() if line.strip() and not line.strip().startswith("#")]
    if not lines:
        raise ValueError("Provided tableau text is empty.")

    # Detect delimiter: comma, semicolon, tab, or whitespace
    first_line = lines[0]
    count_semi = first_line.count(';')
    count_comma = first_line.count(',')
    count_tab = first_line.count('\t')

    if count_semi > count_comma and count_semi > count_tab:
        delimiter = ';'
    elif count_tab > count_comma and count_tab > count_semi:
        delimiter = '\t'
    elif count_comma > 0:
        delimiter = ','
    else:
        delimiter = None  # Fallback to whitespace splitting

    if delimiter:
        reader = csv.reader(lines, delimiter=delimiter)
        raw_rows = [[cell.strip() for cell in row] for row in reader]
    else:
        raw_rows = [line.split() for line in lines]

    raw_rows = [row for row in raw_rows if any(cell for cell in row)]
    if len(raw_rows) < 2:
        raise ValueError("Tableau must contain at least an objective row (Row 0) and 1 constraint row.")

    def is_number(val: str) -> bool:
        try:
            parse_fraction(val)
            return True
        except:
            return False

    # Check if first row contains column headers
    has_headers = any(not is_number(cell) for cell in raw_rows[0])
    if has_headers:
        header_row = raw_rows[0]
        data_rows = raw_rows[1:]
    else:
        header_row = None
        data_rows = raw_rows

    # Check if first column is the Basis column
    has_basis_col = False
    if header_row and header_row[0].lower() in ('basis', 'base', 'basic', 'b', 'var', 'variabile', ''):
        has_basis_col = True
    elif not is_number(data_rows[0][0]) and not is_number(data_rows[1][0]):
        has_basis_col = True

    basic_vars: List[str] = []
    obj_name = "-z"
    if has_basis_col:
        row0_lbl = data_rows[0][0].strip().lower()
        if "w" in row0_lbl:
            obj_name = "-w"
        else:
            obj_name = "-z"
        for r in data_rows[1:]:
            b_name = r[0].strip()
            if not b_name:
                raise ValueError(f"Missing basic variable name in row: {r}")
            basic_vars.append(b_name)
        active_headers = header_row[1:] if header_row else None
        active_data = [r[1:] for r in data_rows]
    else:
        active_headers = header_row
        active_data = data_rows

    # Check position of RHS column (left vs right)
    rhs_col_idx = -1
    if active_headers:
        lower_headers = [h.lower() for h in active_headers]
        if lower_headers and lower_headers[0] in ('rhs', 'term', 'bi', 'b', 'noto', 'termini'):
            rhs_col_idx = 0
        elif lower_headers and lower_headers[-1] in ('rhs', 'term', 'bi', 'b', 'noto', 'termini'):
            rhs_col_idx = len(lower_headers) - 1

    # Check if there is an explicit 'z' column (e.g. z = 1, 0, 0...)
    z_col_idx = None
    if active_headers:
        for idx, h in enumerate(active_headers):
            if idx == rhs_col_idx:
                continue
            if h.lower() in ('z', '-z', 'w', '-w'):
                z_col_idx = idx
                break

    # Extract numerical matrix and variable names
    var_cols: List[int] = []
    var_names: List[str] = []
    num_total_cols = len(active_data[0])

    if rhs_col_idx == -1:
        rhs_col_idx = num_total_cols - 1  # Default RHS to last column if unspecified

    for col in range(num_total_cols):
        if col == rhs_col_idx or col == z_col_idx:
            continue
        var_cols.append(col)
        if active_headers and col < len(active_headers):
            var_names.append(active_headers[col])
        else:
            var_names.append(f"x{len(var_names) + 1}")

    # Build tableau rows: [var_col_1, var_col_2, ..., RHS]
    tableau_data: List[List[Fraction]] = []
    for r_idx, row in enumerate(active_data):
        row_fractions: List[Fraction] = []
        for col in var_cols:
            row_fractions.append(parse_fraction(row[col]))
        # Append RHS at the end
        row_fractions.append(parse_fraction(row[rhs_col_idx]))
        tableau_data.append(row_fractions)

    # If basic_vars were not specified in the first column, automatically detect from identity matrix
    m = len(tableau_data) - 1
    if not basic_vars:
        for i in range(1, m + 1):
            matched_var = None
            for c_idx, v_name in enumerate(var_names):
                is_unit = (tableau_data[i][c_idx] == 1 and
                           all(tableau_data[other_r][c_idx] == 0 for other_r in range(len(tableau_data)) if other_r != i))
                if is_unit:
                    matched_var = v_name
                    break
            basic_vars.append(matched_var if matched_var else f"s{i}")

    # Assign 1-based Bland indices
    var_indices = {name: idx + 1 for idx, name in enumerate(var_names)}
    for b_var in basic_vars:
        if b_var not in var_indices:
            var_indices[b_var] = len(var_indices) + 1

    return SimplexTableau(var_names, var_indices, tableau_data, basic_vars, obj_name=obj_name)


def load_tableau_from_file(file_path: str) -> SimplexTableau:
    """Loads and parses a tableau from a CSV or text file, searching common folders."""
    candidates = [
        file_path,
        os.path.join("csv", file_path),
        os.path.join(os.path.dirname(__file__), "..", "..", "..", "csv", file_path),
        os.path.join(os.path.dirname(__file__), "..", "..", "csv", file_path),
        os.path.join(os.path.dirname(__file__), "..", "examples", file_path)
    ]
    resolved = None
    for cand in candidates:
        if os.path.isfile(cand):
            resolved = cand
            break

    if not resolved:
        raise FileNotFoundError(f"Could not locate tableau file '{file_path}'. Checked candidates: {candidates[:3]}")

    with open(resolved, "r", encoding="utf-8-sig") as f:
        content = f.read()

    return parse_tableau_text(content)


# ---------------------------------------------------------------------
# LINEAR PROGRAM FORMULATION PARSER
# ---------------------------------------------------------------------

class LinearProgram:
    def __init__(self, obj_type: str, obj_coeffs: List[Fraction],
                 constraints: List[Tuple[List[Fraction], str, Fraction]],
                 var_names: Optional[List[str]] = None, force_two_phase: bool = False):
        self.obj_type = obj_type.lower()
        self.obj_coeffs = list(obj_coeffs)
        self.constraints = list(constraints)
        self.num_vars = len(obj_coeffs)
        self.force_two_phase = force_two_phase
        self.var_names = list(var_names) if var_names else [f"x{j+1}" for j in range(self.num_vars)]

    def to_tableau(self) -> Tuple[SimplexTableau, bool]:
        """Converts the Linear Program to standard form and creates the initial SimplexTableau."""
        norm_constraints = []
        for coeffs, rel, rhs in self.constraints:
            c_list = list(coeffs)
            r_rel = rel
            r_rhs = rhs
            if r_rhs < 0:
                c_list = [-c for c in c_list]
                r_rhs = -r_rhs
                if r_rel == "<=":
                    r_rel = ">="
                elif r_rel == ">=":
                    r_rel = "<="
            norm_constraints.append((c_list, r_rel, r_rhs))

        needs_phase1 = self.force_two_phase or any(rel in (">=", "=") for _, rel, _ in norm_constraints)

        # Name slacks as x, continuing from highest existing decision variable index
        existing_x_indices = [
            int(v[1:]) for v in self.var_names if v.startswith("x") and v[1:].isdigit()
        ]
        next_x_idx = max(existing_x_indices, default=len(self.var_names)) + 1
        surplus_counter = 1
        artificial_counter = 1

        slacks: List[str] = []
        surpluses: List[str] = []
        artificials: List[str] = []
        constraint_aux: List[Dict[str, Optional[str]]] = []

        for coeffs, rel, rhs in norm_constraints:
            aux: Dict[str, Optional[str]] = {"slack": None, "surplus": None, "artificial": None}
            if rel == "<=":
                s_name = f"x{next_x_idx}"
                next_x_idx += 1
                slacks.append(s_name)
                aux["slack"] = s_name
                if needs_phase1:
                    y_name = f"y{artificial_counter}"
                    artificial_counter += 1
                    artificials.append(y_name)
                    aux["artificial"] = y_name
            elif rel == ">=":
                e_name = f"e{surplus_counter}"
                surplus_counter += 1
                surpluses.append(e_name)
                aux["surplus"] = e_name
                y_name = f"y{artificial_counter}"
                artificial_counter += 1
                artificials.append(y_name)
                aux["artificial"] = y_name
            elif rel == "=":
                y_name = f"y{artificial_counter}"
                artificial_counter += 1
                artificials.append(y_name)
                aux["artificial"] = y_name
            constraint_aux.append(aux)

        full_var_names = self.var_names + slacks + surpluses + artificials
        var_indices = {name: idx + 1 for idx, name in enumerate(full_var_names)}

        n_full = len(full_var_names)
        matrix_A: List[List[Fraction]] = []
        basic_vars: List[str] = []

        for i, (coeffs, rel, rhs) in enumerate(norm_constraints):
            row = [Fraction(0)] * (n_full + 1)
            for j, c in enumerate(coeffs):
                row[j] = c
            aux = constraint_aux[i]
            if aux["slack"]:
                col = full_var_names.index(aux["slack"])
                row[col] = Fraction(1)
                if not needs_phase1:
                    basic_vars.append(aux["slack"])
            if aux["surplus"]:
                col = full_var_names.index(aux["surplus"])
                row[col] = Fraction(-1)
            if aux["artificial"]:
                col = full_var_names.index(aux["artificial"])
                row[col] = Fraction(1)
                basic_vars.append(aux["artificial"])
            row[-1] = rhs
            matrix_A.append(row)

        orig_obj_row = [Fraction(0)] * (n_full + 1)
        for j, name in enumerate(self.var_names):
            c_idx = full_var_names.index(name)
            if self.obj_type == "max":
                orig_obj_row[c_idx] = -self.obj_coeffs[j]
            else:
                orig_obj_row[c_idx] = self.obj_coeffs[j]

        if needs_phase1:
            phase1_row0 = [Fraction(0)] * (n_full + 1)
            for a_name in artificials:
                phase1_row0[full_var_names.index(a_name)] = Fraction(1)
            for i, aux in enumerate(constraint_aux):
                if aux["artificial"]:
                    for c_idx in range(n_full + 1):
                        phase1_row0[c_idx] -= matrix_A[i][c_idx]

            phase1_tableau_data = [phase1_row0] + matrix_A
            return SimplexTableau(full_var_names, var_indices, phase1_tableau_data, basic_vars, obj_name="-w", orig_row0=orig_obj_row), True
        else:
            # Phase 2 directly
            phase2_row0 = list(orig_obj_row)

            # Canonicalize Row 0
            for r_idx, b_var in enumerate(basic_vars):
                c_idx = full_var_names.index(b_var)
                factor = phase2_row0[c_idx]
                if factor != 0:
                    for col in range(len(phase2_row0)):
                        phase2_row0[col] -= factor * matrix_A[r_idx][col]

            phase2_tableau_data = [phase2_row0] + matrix_A
            return SimplexTableau(full_var_names, var_indices, phase2_tableau_data, basic_vars, obj_name="-z", orig_row0=orig_obj_row), False


def parse_lp_text(text: str, force_two_phase: bool = False) -> LinearProgram:
    """Parses standard LP formulation text into a LinearProgram object."""
    lines = [l.strip() for l in text.strip().splitlines() if l.strip() and not l.strip().startswith("#")]
    if len(lines) < 2:
        raise ValueError("LP formulation requires at least an objective function and 1 constraint.")

    # 1. Objective line (e.g. 'min -3*x1 - 5*x2' or 'max 10*x1 + 20*x2')
    obj_line = lines[0].lower()
    obj_type = "min" if "min" in obj_line else "max"
    obj_expr = re.sub(r'^(min|max|minimize|maximize)\s*([a-z]\s*=)?', '', obj_line).strip()

    # Find all variables in obj expression
    var_names_found = []
    for match in re.finditer(r'([+-]?\s*[\d\./,]*)\s*\*?\s*([a-zA-Z][a-zA-Z0-9_]*)', obj_expr):
        v_name = match.group(2)
        if v_name not in var_names_found:
            var_names_found.append(v_name)

    # Sort variables nicely (e.g. x1, x2, ...)
    def var_sort_key(name: str):
        match = re.match(r'([a-zA-Z]+)(\d+)?', name)
        if match:
            prefix, num = match.groups()
            return (prefix, int(num) if num else 0)
        return (name, 0)

    var_names = sorted(var_names_found, key=var_sort_key) if var_names_found else ["x1", "x2"]

    # Parse objective coefficients
    obj_coeffs = [Fraction(0)] * len(var_names)
    for idx, v_name in enumerate(var_names):
        pattern = rf'([+-]?\s*[\d\./,]*)\s*\*?\s*{re.escape(v_name)}(?![a-zA-Z0-9_])'
        match = re.search(pattern, obj_expr)
        if match:
            raw_c = match.group(1).replace(" ", "")
            if not raw_c or raw_c == "+":
                coeff = Fraction(1)
            elif raw_c == "-":
                coeff = Fraction(-1)
            else:
                coeff = parse_fraction(raw_c)
            obj_coeffs[idx] = coeff

    # Parse constraint lines
    constraints = []
    for c_line in lines[1:]:
        c_clean = c_line.strip()
        if re.search(r's\.?t\.?', c_clean, re.I):
            c_clean = re.sub(r's\.?t\.?:?', '', c_clean, flags=re.I).strip()
        if not c_clean:
            continue
        # Skip non-negativity constraint like 'x1, x2 >= 0'
        if all(v in c_clean for v in var_names) and ">=" in c_clean and "0" in c_clean:
            continue

        rel = None
        for cand in ("<=", ">=", "="):
            if cand in c_clean:
                rel = cand
                break
        if not rel:
            continue

        lhs_str, rhs_str = c_clean.split(rel, 1)
        rhs_val = parse_fraction(rhs_str.strip())
        lhs_coeffs = [Fraction(0)] * len(var_names)

        for idx, v_name in enumerate(var_names):
            pattern = rf'([+-]?\s*[\d\./,]*)\s*\*?\s*{re.escape(v_name)}(?![a-zA-Z0-9_])'
            match = re.search(pattern, lhs_str)
            if match:
                raw_c = match.group(1).replace(" ", "")
                if not raw_c or raw_c == "+":
                    coeff = Fraction(1)
                elif raw_c == "-":
                    coeff = Fraction(-1)
                else:
                    coeff = parse_fraction(raw_c)
                lhs_coeffs[idx] = coeff

        constraints.append((lhs_coeffs, rel, rhs_val))

    return LinearProgram(obj_type, obj_coeffs, constraints, var_names, force_two_phase=force_two_phase)


# ---------------------------------------------------------------------
# DETAILED SIMPLEX SOLVER ENGINE
# ---------------------------------------------------------------------

class SimplexStepDetail:
    """Holds detailed calculations and snapshots for one iteration."""
    def __init__(self, iteration: int, phase_name: str, tableau_before: SimplexTableau):
        self.iteration = iteration
        self.phase_name = phase_name
        self.tableau_before = SimplexTableau(
            tableau_before.var_names, tableau_before.var_indices,
            tableau_before.tableau, tableau_before.basic_vars, tableau_before.obj_name
        )
        self.status: str = "CONTINUE"
        self.candidates: List[Dict[str, Any]] = []
        self.entering_var: Optional[str] = None
        self.entering_col: Optional[int] = None
        self.ratios: List[Dict[str, Any]] = []
        self.min_ratio: Optional[Fraction] = None
        self.leaving_var: Optional[str] = None
        self.pivot_row: Optional[int] = None
        self.pivot_element: Optional[Fraction] = None
        self.tie_detected: bool = False
        self.bland_tie_note: Optional[str] = None
        self.pivot_row_normalization_formula: str = ""
        self.elimination_formulas: List[str] = []
        self.tableau_after: Optional[SimplexTableau] = None


def solve_simplex_with_calculations(
    initial_tab: SimplexTableau,
    phase_name: str = "",
    forbidden_vars: Optional[Set[str]] = None,
    obj_type: str = "min",
    max_iterations: int = 100
) -> Tuple[str, List[SimplexStepDetail], SimplexTableau]:
    """
    Solves the Simplex tableau using Bland's rule and records step-by-step mathematical calculations.
    """
    if forbidden_vars is None:
        forbidden_vars = set()

    tab = SimplexTableau(
        initial_tab.var_names, initial_tab.var_indices,
        initial_tab.tableau, initial_tab.basic_vars, initial_tab.obj_name
    )

    history: List[SimplexStepDetail] = []
    iteration = 0

    while iteration < max_iterations:
        step = SimplexStepDetail(iteration, phase_name, tab)

        # -------------------------------------------------------------
        # STEP 1: Entering Variable Selection (Bland's Rule)
        # -------------------------------------------------------------
        candidates = []
        for col_idx, var_name in enumerate(tab.var_names):
            if var_name in forbidden_vars:
                continue
            coeff = tab.tableau[0][col_idx]
            if coeff < 0:
                v_idx = tab.var_indices[var_name]
                candidates.append({
                    "name": var_name,
                    "index": v_idx,
                    "col": col_idx,
                    "coeff": coeff
                })

        step.candidates = candidates

        if not candidates:
            step.status = "OPTIMAL"
            step.tableau_after = SimplexTableau(tab.var_names, tab.var_indices, tab.tableau, tab.basic_vars, tab.obj_name)
            history.append(step)
            return "OPTIMAL", history, tab

        # Bland's Rule: choose smallest variable index
        candidates.sort(key=lambda c: c["index"])
        entering = candidates[0]
        step.entering_var = entering["name"]
        step.entering_col = entering["col"]

        # -------------------------------------------------------------
        # STEP 2: Minimum Ratio Test (Leaving Variable)
        # -------------------------------------------------------------
        ratios = []
        for r in range(1, tab.num_rows()):
            b_var = tab.basic_vars[r - 1]
            b_idx = tab.var_indices[b_var]
            a_ij = tab.tableau[r][step.entering_col]
            rhs = tab.tableau[r][-1]

            if a_ij > 0:
                ratio = rhs / a_ij
                ratios.append({
                    "row": r,
                    "basic": b_var,
                    "index": b_idx,
                    "ratio": ratio,
                    "valid": True,
                    "rhs": rhs,
                    "coeff": a_ij
                })
            else:
                reason = "coeff = 0" if a_ij == 0 else "coeff < 0"
                ratios.append({
                    "row": r,
                    "basic": b_var,
                    "index": b_idx,
                    "ratio": None,
                    "valid": False,
                    "reason": reason,
                    "rhs": rhs,
                    "coeff": a_ij
                })

        step.ratios = ratios
        valid_ratios = [r for r in ratios if r["valid"]]

        if not valid_ratios:
            step.status = "UNBOUNDED"
            step.tableau_after = SimplexTableau(tab.var_names, tab.var_indices, tab.tableau, tab.basic_vars, tab.obj_name)
            history.append(step)
            return "UNBOUNDED", history, tab

        min_ratio_val = min(r["ratio"] for r in valid_ratios)
        tied = [r for r in valid_ratios if r["ratio"] == min_ratio_val]
        step.min_ratio = min_ratio_val

        if len(tied) > 1:
            step.tie_detected = True
            tied.sort(key=lambda r: r["index"])
            leaving_item = tied[0]
            step.bland_tie_note = (
                f"Degenerate tie between rows {[t['row'] for t in tied]} (min θ = {format_fraction(min_ratio_val)}). "
                f"Bland's Rule selects basic variable with smallest index: '{leaving_item['basic']}' (Index {leaving_item['index']})."
            )
        else:
            leaving_item = tied[0]

        step.leaving_var = leaving_item["basic"]
        step.pivot_row = leaving_item["row"]
        pivot_elem = tab.tableau[step.pivot_row][step.entering_col]
        step.pivot_element = pivot_elem

        # -------------------------------------------------------------
        # STEP 3: Elementary Row Operations
        # -------------------------------------------------------------
        # 1. Normalize pivot row
        step.pivot_row_normalization_formula = (
            f"R_{step.pivot_row} \\leftarrow \\frac{{1}}{{{format_fraction(pivot_elem)}}} \\times R_{step.pivot_row}"
            if pivot_elem != 1 else f"R_{step.pivot_row} \\text{{ already normalized (pivot = 1)}}"
        )
        tab.tableau[step.pivot_row] = [x / pivot_elem for x in tab.tableau[step.pivot_row]]

        # 2. Eliminate other rows
        for r in range(tab.num_rows()):
            if r != step.pivot_row:
                factor = tab.tableau[r][step.entering_col]
                if factor != 0:
                    tab.tableau[r] = [
                        tab.tableau[r][c] - factor * tab.tableau[step.pivot_row][c]
                        for c in range(len(tab.tableau[r]))
                    ]
                    row_name = tab.obj_name if r == 0 else f"R_{r}"
                    sign = "-" if factor > 0 else "+"
                    abs_factor_str = format_fraction(abs(factor))
                    if abs_factor_str == "1":
                        formula = f"{row_name} \\leftarrow {row_name} {sign} R_{step.pivot_row}"
                    else:
                        formula = f"{row_name} \\leftarrow {row_name} {sign} ({abs_factor_str}) \\times R_{step.pivot_row}"
                    step.elimination_formulas.append(formula)

        # Update basis
        tab.basic_vars[step.pivot_row - 1] = step.entering_var
        step.tableau_after = SimplexTableau(tab.var_names, tab.var_indices, tab.tableau, tab.basic_vars, tab.obj_name)
        history.append(step)
        iteration += 1

    return "ITERATION_LIMIT", history, tab


def run_full_solver(tab: SimplexTableau, force_two_phase: bool = False, obj_type: str = "min") -> Dict[str, Any]:
    """Executes Simplex or Two-Phase Simplex and returns complete calculations and report data."""
    m = tab.num_rows() - 1
    needs_phase1 = force_two_phase or tab.obj_name == "-w" or any(v.startswith("y") for v in tab.var_names)

    # Check if constraints have negative RHS
    has_negative_rhs = any(tab.tableau[r][-1] < 0 for r in range(1, tab.num_rows()))
    if has_negative_rhs:
        needs_phase1 = True

    if not needs_phase1:
        status, steps, final_tab = solve_simplex_with_calculations(tab, obj_type=obj_type)
        return build_solution_dict(status, [("Simplex", steps)], final_tab, tab.var_names, obj_type)

    # -----------------------------------------------------------------
    # TWO-PHASE SIMPLEX
    # -----------------------------------------------------------------
    original_var_names = [v for v in tab.var_names if not v.startswith("y")]
    n_orig = len(original_var_names)

    # 1. Normalize constraint rows so all RHS >= 0
    norm_matrix_A: List[List[Fraction]] = []
    for r in range(1, tab.num_rows()):
        row = list(tab.tableau[r][:n_orig]) + [tab.tableau[r][-1]]
        if row[-1] < 0:
            row = [-val for val in row]
        norm_matrix_A.append(row)

    if tab.orig_row0 is not None:
        orig_row0 = [Fraction(0)] * (n_orig + 1)
        for j, name in enumerate(original_var_names):
            if name in tab.var_names:
                c_idx = tab.var_names.index(name)
                orig_row0[j] = tab.orig_row0[c_idx]
        orig_row0[-1] = tab.orig_row0[-1]
    else:
        orig_row0 = list(tab.tableau[0][:n_orig]) + [tab.tableau[0][-1]]

    # 2. Setup Phase 1 with artificial variables y1, ..., ym
    artificials = [f"y{i+1}" for i in range(m)]
    phase1_var_names = original_var_names + artificials
    phase1_var_indices = {name: idx + 1 for idx, name in enumerate(phase1_var_names)}
    n_full = len(phase1_var_names)

    phase1_rows: List[List[Fraction]] = []
    for i in range(m):
        row = [Fraction(0)] * (n_full + 1)
        for j in range(n_orig):
            row[j] = norm_matrix_A[i][j]
        row[n_orig + i] = Fraction(1)  # Identity for artificial variable
        row[-1] = norm_matrix_A[i][-1] # RHS
        phase1_rows.append(row)

    # Row 0: -w = -sum(y_i)
    phase1_row0 = [Fraction(0)] * (n_full + 1)
    for col in range(n_full + 1):
        for i in range(m):
            phase1_row0[col] -= phase1_rows[i][col]
    for i in range(m):
        phase1_row0[n_orig + i] = Fraction(0)

    phase1_tableau_data = [phase1_row0] + phase1_rows
    phase1_basic_vars = list(artificials)
    phase1_tab = SimplexTableau(phase1_var_names, phase1_var_indices, phase1_tableau_data,
                                phase1_basic_vars, obj_name="-w")

    p1_status, p1_steps, solved_p1_tab = solve_simplex_with_calculations(
        phase1_tab, phase_name="Phase 1", forbidden_vars=set(artificials), obj_type="min"
    )

    if p1_status != "OPTIMAL":
        return build_solution_dict(p1_status, [("Phase 1", p1_steps)], solved_p1_tab, original_var_names, obj_type)

    w_opt = solved_p1_tab.tableau[0][-1]
    if w_opt != 0:
        res = build_solution_dict("INFEASIBLE", [("Phase 1", p1_steps)], solved_p1_tab, original_var_names, obj_type)
        res["infeasible_reason"] = f"Minimum sum of artificial variables W = {format_fraction(-w_opt)} > 0"
        return res

    # Drive degenerate artificial variables out of basis
    redundant_rows = []
    for r_idx in range(1, solved_p1_tab.num_rows()):
        b_var = solved_p1_tab.basic_vars[r_idx - 1]
        if b_var in artificials:
            driven_out = False
            for c_idx, v_name in enumerate(solved_p1_tab.var_names):
                if v_name not in artificials and solved_p1_tab.tableau[r_idx][c_idx] != 0:
                    p_val = solved_p1_tab.tableau[r_idx][c_idx]
                    solved_p1_tab.tableau[r_idx] = [x / p_val for x in solved_p1_tab.tableau[r_idx]]
                    for r_other in range(solved_p1_tab.num_rows()):
                        if r_other != r_idx:
                            factor = solved_p1_tab.tableau[r_other][c_idx]
                            solved_p1_tab.tableau[r_other] = [
                                solved_p1_tab.tableau[r_other][c] - factor * solved_p1_tab.tableau[r_idx][c]
                                for c in range(len(solved_p1_tab.tableau[r_other]))
                            ]
                    solved_p1_tab.basic_vars[r_idx - 1] = v_name
                    driven_out = True
                    break
            if not driven_out:
                redundant_rows.append(r_idx)

    # 3. Setup Phase 2
    non_art_cols = [c for c, name in enumerate(solved_p1_tab.var_names) if name not in artificials]
    phase2_var_names = [solved_p1_tab.var_names[c] for c in non_art_cols]
    phase2_var_indices = {name: idx + 1 for idx, name in enumerate(phase2_var_names)}

    phase2_rows = []
    phase2_basic_vars = []
    for r in range(1, solved_p1_tab.num_rows()):
        if r in redundant_rows:
            continue
        row = [solved_p1_tab.tableau[r][c] for c in non_art_cols] + [solved_p1_tab.tableau[r][-1]]
        phase2_rows.append(row)
        phase2_basic_vars.append(solved_p1_tab.basic_vars[r - 1])

    # Reconstitute original Row 0 and canonicalize
    phase2_row0 = [Fraction(0)] * (len(phase2_var_names) + 1)
    for j, name in enumerate(phase2_var_names):
        if name in original_var_names:
            orig_col = original_var_names.index(name)
            phase2_row0[j] = orig_row0[orig_col]
    phase2_row0[-1] = orig_row0[-1]

    for r_idx, b_var in enumerate(phase2_basic_vars):
        c_idx = phase2_var_names.index(b_var)
        factor = phase2_row0[c_idx]
        if factor != 0:
            for col in range(len(phase2_row0)):
                phase2_row0[col] -= factor * phase2_rows[r_idx][col]

    phase2_tableau_data = [phase2_row0] + phase2_rows
    phase2_tab = SimplexTableau(phase2_var_names, phase2_var_indices, phase2_tableau_data,
                                phase2_basic_vars, obj_name="-z")

    p2_status, p2_steps, final_tab = solve_simplex_with_calculations(
        phase2_tab, phase_name="Phase 2", obj_type=obj_type
    )

    all_phases = [("Phase 1", p1_steps), ("Phase 2", p2_steps)]
    return build_solution_dict(p2_status, all_phases, final_tab, original_var_names, obj_type)


def build_solution_dict(status: str, phase_steps: List[Tuple[str, List[SimplexStepDetail]]],
                        final_tab: SimplexTableau, original_var_names: List[str],
                        obj_type: str) -> Dict[str, Any]:
    """Assembles the final dictionary with solution vector, objective value, and full steps."""
    result: Dict[str, Any] = {
        "status": status,
        "obj_type": obj_type,
        "phases": phase_steps,
        "final_tableau": final_tab,
        "solution": None
    }

    if status == "OPTIMAL":
        raw_rhs = final_tab.tableau[0][-1]
        opt_z = get_display_obj_val(raw_rhs, obj_type, final_tab.obj_name)
        sol_vector = {name: Fraction(0) for name in original_var_names}
        for r_idx in range(1, final_tab.num_rows()):
            b_var = final_tab.basic_vars[r_idx - 1]
            if b_var in sol_vector:
                sol_vector[b_var] = final_tab.tableau[r_idx][-1]

        result["solution"] = {
            "optimal_z": opt_z,
            "optimal_z_approx": format_approx(opt_z),
            "raw_rhs": raw_rhs,
            "solution_vector": sol_vector,
            "basic_vars": list(final_tab.basic_vars)
        }
    return result


# ---------------------------------------------------------------------
# MARKDOWN REPORT GENERATOR
# ---------------------------------------------------------------------

def generate_markdown_report(result: Dict[str, Any], title: str = "Simplex Tableau Solution Report") -> str:
    """Generates a detailed, publication-ready Markdown file with all iterations and calculations."""
    lines: List[str] = [
        f"# {title}",
        "",
        "> [!NOTE]",
        "> Calculations performed using exact rational arithmetic (`fractions.Fraction`) with Bland's Smallest-Subscript Anti-Cycling Rule.",
        f"> **Objective Formulation**: Universal Minimization ($\\min z$). Row 0 Identity: $-z = \\text{{RHS}}_0 \\implies z^* = -\\text{{RHS}}_0$.",
        ""
    ]

    # Iteration Walkthrough
    for phase_name, steps in result["phases"]:
        lines.append(f"## {phase_name} Walkthrough")
        lines.append("")

        for step in steps:
            lines.append(f"### Iteration {step.iteration} ({phase_name})")
            lines.append("")
            lines.append("#### Tableau State:")
            lines.append(step.tableau_before.to_markdown_table())
            lines.append("")

            # Basic Solution
            basic_parts = [
                f"{step.tableau_before.basic_vars[i-1]} = {format_fraction(step.tableau_before.tableau[i][-1])}"
                for i in range(1, step.tableau_before.num_rows())
            ]
            raw_rhs = step.tableau_before.tableau[0][-1]
            z_curr = -raw_rhs
            lines.append(f"**Current Basic Solution**: {', '.join(basic_parts)}")
            lines.append(f"**Objective Equation**: ${step.tableau_before.obj_name} = {format_fraction(raw_rhs)} \\implies z = {format_fraction(z_curr)}$")
            lines.append("")

            if step.status == "OPTIMAL":
                lines.append("> [!TIP]")
                lines.append(f"> **Optimality Reached**: All reduced costs in Row 0 are $\\ge 0$. The current tableau is optimal for {phase_name}.")
                lines.append("")
                continue

            if step.status == "UNBOUNDED":
                lines.append("> [!WARNING]")
                lines.append(f"> **Unbounded Problem**: All constraint coefficients in entering column `{step.entering_var}` are $\\le 0$. The objective can be improved indefinitely.")
                lines.append("")
                continue

            # Entering Variable
            lines.append("#### 1. Entering Variable Selection (Bland's Rule)")
            cand_str = ", ".join(f"`{c['name']}` (cost: ${format_fraction(c['coeff'])}$, Bland index: {c['index']})" for c in step.candidates)
            lines.append(f"- **Candidates with negative reduced costs ($c_j < 0$)**: {cand_str}")
            lines.append(f"- **Selected Entering Variable**: `{step.entering_var}` (smallest Bland index).")
            lines.append("")

            # Ratio Test
            lines.append(f"#### 2. Minimum Ratio Test (Leaving Variable for Column `{step.entering_var}`)")
            lines.append("| Row | Basis | RHS | Pivot Col Entry | Ratio $\\theta = \\text{RHS} / a_{ij}$ | Valid? |")
            lines.append("| --- | --- | --- | --- | --- | --- |")
            for r in step.ratios:
                if r["valid"]:
                    ratio_str = f"$\\frac{{{format_fraction(r['rhs'])}}}{{{format_fraction(r['coeff'])}}} = {format_fraction(r['ratio'])}$"
                    valid_str = "Yes"
                else:
                    ratio_str = f"Skipped ({r['reason']})"
                    valid_str = "No"
                lines.append(f"| $R_{r['row']}$ | `{r['basic']}` | ${format_fraction(r['rhs'])}$ | ${format_fraction(r['coeff'])}$ | {ratio_str} | {valid_str} |")
            lines.append("")
            lines.append(f"- **Minimum Ratio**: $\\theta^* = {format_fraction(step.min_ratio)}$")
            lines.append(f"- **Selected Leaving Variable**: `{step.leaving_var}` ($R_{step.pivot_row}$)")
            if step.tie_detected:
                lines.append(f"- > [!IMPORTANT]\n  > {step.bland_tie_note}")
            lines.append(f"- **Pivot Element**: $(R_{step.pivot_row}, \\text{{{step.entering_var}}}) = {format_fraction(step.pivot_element)}$")
            lines.append("")

            # Row Operations
            lines.append("#### 3. Elementary Row Operations (Pivot)")
            lines.append(f"- **Normalize Pivot Row**: ${step.pivot_row_normalization_formula}$")
            if step.elimination_formulas:
                lines.append("- **Row Eliminations**:")
                for elim in step.elimination_formulas:
                    lines.append(f"  - ${elim}$")
            lines.append(f"- **Basis Updated**: Row $R_{step.pivot_row}$ basis changes from `{step.leaving_var}` $\\rightarrow$ `{step.entering_var}`.")
            lines.append("")
            lines.append("---")
            lines.append("")

    # Final Summary
    lines.append("## Final Solution Summary")
    lines.append("")
    lines.append(f"- **Termination Status**: `{result['status']}`")
    if result["status"] == "OPTIMAL":
        sol = result["solution"]
        lines.append(f"- **Optimal Value**: $\\min z^* = {sol['optimal_z_approx']}$")
        lines.append(f"- **Row 0 RHS**: $-z = {format_fraction(sol['raw_rhs'])}$")
        lines.append("")
        lines.append("### Optimal Variable Assignment:")
        lines.append("| Variable | Type | Value |")
        lines.append("| --- | --- | --- |")
        for v_name, val in sol["solution_vector"].items():
            v_type = "Basic" if v_name in sol["basic_vars"] else "Non-Basic"
            lines.append(f"| `{v_name}` | {v_type} | ${format_approx(val)}$ |")
        lines.append("")
    elif result["status"] == "INFEASIBLE":
        lines.append(f"- **Reason**: {result.get('infeasible_reason', 'Artificial variables could not be driven to zero.')}")
    elif result["status"] == "UNBOUNDED":
        lines.append("- **Reason**: An entering column has no positive constraint entries.")
    lines.append("")

    # SciPy Integration Section in Markdown
    scipy_info = result.get("scipy")
    if scipy_info:
        lines.append("## Python Library Integration (SciPy `scipy.optimize.linprog`)")
        lines.append("")
        if scipy_info.get("installed"):
            lines.append(f"- **Library Version**: SciPy v{scipy_info.get('version')}")
            lines.append(f"- **Solver Method**: `{scipy_info.get('method')}`")
            lines.append(f"- **Termination Status**: `{scipy_info.get('message')}` (Code {scipy_info.get('status_code')})")
            lines.append(f"- **SciPy Iterations**: `{scipy_info.get('nit')}`")
            if scipy_info.get("fun") is not None:
                lines.append(f"- **SciPy Optimal Value (Min $z^*$)**: `{scipy_info.get('fun')}`")
            if scipy_info.get("x"):
                x_str = ", ".join(f"`{k}` = {v}" for k, v in scipy_info["x"].items())
                lines.append(f"- **SciPy Optimal Vector**: {x_str}")
        else:
            lines.append("> [!NOTE]")
            lines.append(f"> **SciPy Library Status**: {scipy_info.get('message')}")
            lines.append(f"> **How to Install**:")
            lines.append(f"> ```bash\n> {scipy_info.get('install_guide')}\n> ```")
        lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------
# TERMINAL TEXT FORMATTER
# ---------------------------------------------------------------------

def print_terminal_walkthrough(result: Dict[str, Any]):
    """Prints a clear, structured terminal walkthrough of all iterations."""
    print("=" * 82)
    print(" SIMPLEX TABLEAU STEP-BY-STEP CALCULATION ENGINE".center(82))
    print("=" * 82)

    for phase_name, steps in result["phases"]:
        for step in steps:
            print(f"\n>>> [{phase_name}] Iteration {step.iteration}")
            print(step.tableau_before.to_ascii_table())

            basic_parts = [
                f"{step.tableau_before.basic_vars[i-1]} = {format_fraction(step.tableau_before.tableau[i][-1])}"
                for i in range(1, step.tableau_before.num_rows())
            ]
            raw_rhs = step.tableau_before.tableau[0][-1]
            print(f"    Basis: {', '.join(basic_parts)}")
            print(f"    Objective: {step.tableau_before.obj_name} = {format_fraction(raw_rhs)}  =>  z = {format_fraction(-raw_rhs)}")

            if step.status == "OPTIMAL":
                print(f"\n[OPTIMAL] All Row 0 reduced costs are >= 0. {phase_name} reached optimality!")
                continue
            if step.status == "UNBOUNDED":
                print(f"\n[UNBOUNDED] Column '{step.entering_var}' has no positive constraint entries.")
                continue

            # Entering
            print(f"\n[1] Bland Entering Variable: '{step.entering_var}' (Index {step.tableau_before.var_indices[step.entering_var]})")
            cands_desc = [f"{c['name']} (cost: {format_fraction(c['coeff'])})" for c in step.candidates]
            print(f"    Candidates: {', '.join(cands_desc)}")

            # Leaving
            print(f"[2] Minimum Ratio Test for '{step.entering_var}':")
            for r in step.ratios:
                if r["valid"]:
                    print(f"    Row {r['row']} ({r['basic']}): {format_fraction(r['rhs'])} / {format_fraction(r['coeff'])} = {format_fraction(r['ratio'])}")
                else:
                    print(f"    Row {r['row']} ({r['basic']}): Skipped ({r['reason']})")
            print(f"    Selected Leaving Variable: '{step.leaving_var}' (Row {step.pivot_row}, min θ = {format_fraction(step.min_ratio)})")
            if step.tie_detected:
                print(f"    [Bland Tie Break] {step.bland_tie_note}")
            print(f"    Pivot Element: (R{step.pivot_row}, {step.entering_var}) = {format_fraction(step.pivot_element)}")

            # Row ops
            print("[3] Elementary Row Operations:")
            print(f"    Normalize: {step.pivot_row_normalization_formula}")
            for elim in step.elimination_formulas:
                print(f"    Eliminate: {elim}")
            print(f"    Updated Basis: Row R{step.pivot_row} is now '{step.entering_var}'")

    print("\n" + "=" * 82)
    print(" FINAL RESULTS".center(82))
    print("=" * 82)
    print(f" Status: {result['status']}")
    if result["status"] == "OPTIMAL":
        sol = result["solution"]
        print(f" Optimal Objective (Min z*): {sol['optimal_z_approx']}")
        print(f" Row 0 RHS (-z): {format_fraction(sol['raw_rhs'])}")
        print(" Optimal Solution Vector:")
        for v_name, val in sol["solution_vector"].items():
            b_tag = "BASIC" if v_name in sol["basic_vars"] else "NON-BASIC"
            print(f"   {v_name:<6} = {format_approx(val):<12} ({b_tag})")
    elif result["status"] == "INFEASIBLE":
        print(f" Infeasible: {result.get('infeasible_reason', '')}")

    scipy_info = result.get("scipy")
    if scipy_info:
        print("\n" + "=" * 82)
        print(" PYTHON LIBRARY (SciPy linprog) STATUS".center(82))
        print("=" * 82)
        if scipy_info.get("installed"):
            print(f" SciPy Version: v{scipy_info.get('version')}")
            print(f" Solver Method: {scipy_info.get('method')}")
            print(f" Status Message: {scipy_info.get('message')} (Code {scipy_info.get('status_code')})")
            print(f" SciPy Iterations: {scipy_info.get('nit')}")
            if scipy_info.get("fun") is not None:
                print(f" SciPy Optimal Value (Min z*): {scipy_info.get('fun')}")
            if scipy_info.get("x"):
                x_parts = [f"{k} = {v}" for k, v in scipy_info["x"].items()]
                print(f" SciPy Solution Vector: {', '.join(x_parts)}")
        else:
            print(f" [!] {scipy_info.get('message')}")
            print(f" To install SciPy in your Python environment, run:")
            print(f"   {scipy_info.get('install_guide')}")


# ---------------------------------------------------------------------
# SCIPY INTEGRATION & VERIFICATION
# ---------------------------------------------------------------------

def solve_with_scipy(tab: SimplexTableau, lp: Optional[LinearProgram] = None, obj_type: str = "min") -> Dict[str, Any]:
    """
    Executes linear programming via SciPy's official library (scipy.optimize.linprog)
    and returns a structured comparison dictionary.
    """
    try:
        import scipy
        from scipy.optimize import linprog
        import numpy as np
    except ImportError:
        return {
            "installed": False,
            "message": "SciPy is not installed in the active Python environment.",
            "install_guide": "pip3 install --user scipy numpy"
        }

    try:
        if lp is not None:
            c = [float(val) if lp.obj_type == "min" else -float(val) for val in lp.obj_coeffs]
            A_ub, b_ub, A_eq, b_eq = [], [], [], []
            for coeffs, rel, rhs in lp.constraints:
                c_float = [float(v) for v in coeffs]
                r_float = float(rhs)
                if rel == "<=":
                    A_ub.append(c_float)
                    b_ub.append(r_float)
                elif rel == ">=":
                    A_ub.append([-v for v in c_float])
                    b_ub.append(-r_float)
                elif rel == "=":
                    A_eq.append(c_float)
                    b_eq.append(r_float)

            res = linprog(
                c=c,
                A_ub=A_ub if A_ub else None,
                b_ub=b_ub if b_ub else None,
                A_eq=A_eq if A_eq else None,
                b_eq=b_eq if b_eq else None,
                bounds=(0, None),
                method="highs"
            )
            var_names = lp.var_names
            raw_rhs_0 = 0.0
        else:
            # Formulate from Tableau constraint system Ax = b
            n_vars = tab.num_cols()
            if tab.orig_row0 is not None:
                c = [float(tab.orig_row0[j]) for j in range(n_vars)]
            else:
                c = [float(tab.tableau[0][j]) for j in range(n_vars)]

            A_eq = [[float(tab.tableau[i][j]) for j in range(n_vars)] for i in range(1, tab.num_rows())]
            b_eq = [float(tab.tableau[i][-1]) for i in range(1, tab.num_rows())]
            raw_rhs_0 = float(tab.tableau[0][-1])

            res = linprog(
                c=c,
                A_eq=A_eq,
                b_eq=b_eq,
                bounds=(0, None),
                method="highs"
            )
            var_names = tab.var_names

        sol_dict = {}
        if res.x is not None:
            for idx, name in enumerate(var_names):
                if idx < len(res.x):
                    sol_dict[name] = round(float(res.x[idx]), 6)

        calc_fun = float(res.fun) if res.fun is not None else None
        # In a canonical tableau, -z = RHS_0 => true z = fun - RHS_0
        if lp is None and calc_fun is not None and raw_rhs_0 != 0:
            calc_fun = calc_fun - raw_rhs_0

        return {
            "installed": True,
            "version": scipy.__version__,
            "success": res.success,
            "status_code": res.status,
            "message": res.message,
            "fun": calc_fun,
            "nit": int(res.nit),
            "x": sol_dict,
            "method": "scipy.optimize.linprog (method='highs')"
        }
    except Exception as ex:
        return {
            "installed": True,
            "error": str(ex)
        }


def verify_with_scipy(tab: SimplexTableau, result: Dict[str, Any]):
    """Compares the custom Simplex solution against SciPy linprog output."""
    scipy_info = result.get("scipy")
    if not scipy_info or not scipy_info.get("installed"):
        print("\n[Cross-Verification] SciPy is not installed; verification skipped.")
        return
    if not scipy_info.get("success"):
        print(f"\n[Cross-Verification] SciPy could not find an optimal solution: {scipy_info.get('message')}")
        return
    if result.get("status") == "OPTIMAL" and result.get("solution"):
        our_z = float(result["solution"]["optimal_z"])
        scipy_z = scipy_info.get("fun")
        if scipy_z is not None:
            if abs(our_z - scipy_z) < 1e-5:
                print(f"\n[Cross-Verification] MATCH CONFIRMED: Simplex Engine ({our_z:.6g}) == SciPy ({scipy_z:.6g})")
            else:
                print(f"\n[Cross-Verification] DISCREPANCY: Simplex Engine ({our_z}) vs SciPy ({scipy_z})")


# ---------------------------------------------------------------------
# CLI ENTRY POINT
# ---------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Simplex Tableau Solver with Step-by-Step Calculations & Markdown Export",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Solve from CSV file
  python3 solve_tableau.py --csv tableau.csv

  # Solve from raw text passed inline
  python3 solve_tableau.py --raw "Basis,RHS,x1,x2\\n-z,0,-3,-5\\nx3,4,1,0\\nx4,12,0,2\\nx5,18,3,2"

  # Force Two-Phase Simplex and output a Markdown report file
  python3 solve_tableau.py --csv tableau.csv --two-phase --output-md solution.md

  # Solve LP formulation from prompt text
  python3 solve_tableau.py --lp "min -3*x1 - 5*x2 s.t. x1 <= 4; 2*x2 <= 12; 3*x1 + 2*x2 <= 18; x1, x2 >= 0"
        """
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--csv", type=str, help="Path to input CSV file")
    group.add_argument("--raw", type=str, help="Raw multiline tableau CSV/text string")
    group.add_argument("--stdin", action="store_true", help="Read raw tableau from standard input")
    group.add_argument("--lp", type=str, help="Linear Program formulation text")

    parser.add_argument("--two-phase", action="store_true", help="Force Two-Phase Simplex even if an initial basis exists")
    parser.add_argument("--output-md", type=str, help="Path to write the complete Markdown report (.md)")
    parser.add_argument("--format", choices=["text", "markdown", "json"], default="text", help="Terminal output format")
    parser.add_argument("--verify-scipy", action="store_true", help="Cross-verify against scipy.optimize.linprog if installed")

    args = parser.parse_args()

    # Load / Parse Tableau
    try:
        if args.csv:
            tableau = load_tableau_from_file(args.csv)
            obj_type = "min"
        elif args.raw:
            tableau = parse_tableau_text(args.raw)
            obj_type = "min"
        elif args.stdin:
            tableau = parse_tableau_text(sys.stdin.read())
            obj_type = "min"
        elif args.lp:
            lp = parse_lp_text(args.lp, force_two_phase=args.two_phase)
            tableau, forced_p1 = lp.to_tableau()
            obj_type = lp.obj_type
            if forced_p1:
                args.two_phase = True
    except Exception as e:
        print(f"Error loading input tableau: {e}", file=sys.stderr)
        sys.exit(1)

    # Solve
    result = run_full_solver(tableau, force_two_phase=args.two_phase, obj_type=obj_type)

    # Execute SciPy linear programming solver / status check
    lp_obj = lp if args.lp else None
    result["scipy"] = solve_with_scipy(tableau, lp=lp_obj, obj_type=obj_type)

    # Output according to format
    if args.format == "text":
        print_terminal_walkthrough(result)
    elif args.format == "markdown":
        md_text = generate_markdown_report(result)
        print(md_text)
    elif args.format == "json":
        # Serialize JSON summary
        json_obj = {
            "status": result["status"],
            "obj_type": result["obj_type"],
            "solution": {
                "optimal_z": str(result["solution"]["optimal_z"]),
                "solution_vector": {k: str(v) for k, v in result["solution"]["solution_vector"].items()},
                "basic_vars": result["solution"]["basic_vars"]
            } if result["solution"] else None,
            "num_phases": len(result["phases"]),
            "iterations_per_phase": {p[0]: len(p[1]) for p in result["phases"]}
        }
        print(json.dumps(json_obj, indent=2))

    # Optional Markdown file export
    if args.output_md:
        md_content = generate_markdown_report(result, title=f"Simplex Tableau Solution: {os.path.basename(args.csv or 'Custom Tableau')}")
        with open(args.output_md, "w", encoding="utf-8") as f:
            f.write(md_content)
        print(f"\n[Report Exported] Solution markdown report successfully written to: {args.output_md}")

    # Optional Scipy Verification
    if args.verify_scipy:
        verify_with_scipy(tableau, result)


if __name__ == "__main__":
    main()
