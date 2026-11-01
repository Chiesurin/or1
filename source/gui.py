#!/usr/bin/env python3
"""
Simplex Method with Bland's Rule - Standalone Graphical Web Application
========================================================================
Interactive graphical solver for Linear Programming with:
- Exact rational fraction arithmetic (fractions.Fraction)
- Bland's smallest-subscript anti-cycling rule
- Two-Phase Simplex with artificial variables y1, y2, ... in yellow
- Full step-by-step tableau inspection with custom color coding:
    * Basic columns: Red
    * Non-basic columns: Blue
    * Artificial y variables: Yellow
    * Reduced costs: Green (Row 0)
    * RHS: Violet / Purple
    * Pivot element: Gold / Yellow
- -z = RHS convention for minimization
- Zero external dependencies (Python standard library only).

Usage:
  python3 source/gui.py
  python3 source/gui.py csv/tableau.csv --cli
"""

import os
import sys
import json
import glob
import socket
import urllib.parse
import webbrowser
from http.server import HTTPServer, BaseHTTPRequestHandler
from fractions import Fraction
from typing import List, Tuple, Dict, Optional, Set
import re
import csv

# Directory and Path Configuration for restructured folders (source/ and csv/)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR) if os.path.basename(SCRIPT_DIR) == "source" else SCRIPT_DIR

if os.path.isdir(os.path.join(PROJECT_ROOT, "csv")):
    CSV_DIR = os.path.join(PROJECT_ROOT, "csv")
elif os.path.isdir(os.path.join(SCRIPT_DIR, "csv")):
    CSV_DIR = os.path.join(SCRIPT_DIR, "csv")
else:
    CSV_DIR = PROJECT_ROOT


def resolve_csv_path(filename: str, for_writing: bool = False) -> str:
    """Resolves a CSV filename or path to an absolute path across csv/ and project root."""
    if not filename:
        filename = "tableau.csv"

    if os.path.isabs(filename):
        return filename

    # Check directly from current working directory
    if os.path.exists(filename):
        return os.path.abspath(filename)

    # Check in CSV_DIR
    in_csv_dir = os.path.join(CSV_DIR, filename)
    if os.path.exists(in_csv_dir):
        return in_csv_dir

    base_in_csv = os.path.join(CSV_DIR, os.path.basename(filename))
    if os.path.exists(base_in_csv):
        return base_in_csv

    # Check in PROJECT_ROOT
    in_root = os.path.join(PROJECT_ROOT, filename)
    if os.path.exists(in_root):
        return in_root

    base_in_root = os.path.join(PROJECT_ROOT, os.path.basename(filename))
    if os.path.exists(base_in_root):
        return base_in_root

    # For saving new files, prefer CSV_DIR
    if for_writing:
        if os.path.isdir(CSV_DIR):
            return os.path.join(CSV_DIR, os.path.basename(filename))
        return os.path.join(PROJECT_ROOT, os.path.basename(filename))

    return in_csv_dir if os.path.isdir(CSV_DIR) else filename


def get_available_csv_files() -> List[str]:
    """Finds all available CSV files in csv/ directory, project root, and cwd."""
    found: Dict[str, str] = {}

    search_dirs = [CSV_DIR, PROJECT_ROOT, os.getcwd()]
    for d in search_dirs:
        if os.path.isdir(d):
            for path in sorted(glob.glob(os.path.join(d, "*.csv"))):
                base = os.path.basename(path)
                if base not in found:
                    found[base] = path

    file_list = sorted(list(found.keys()))
    if "tableau.csv" in file_list:
        file_list.remove("tableau.csv")
        file_list.insert(0, "tableau.csv")
    return file_list


# =====================================================================
# SIMPLEX ALGORITHM & LINEAR PROGRAM SOLVER ENGINE
# =====================================================================


def parse_fraction(val_str: str) -> Fraction:
    """Parses a string into an exact Fraction. Supports integers, fractions (a/b), and decimals (with . or ,)."""
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
                 tableau: List[List[Fraction]], basic_vars: List[str], obj_name: str = "-z"):
        self.var_names = list(var_names)
        self.var_indices = dict(var_indices)
        self.tableau = [list(row) for row in tableau]
        self.basic_vars = list(basic_vars)
        self.obj_name = obj_name

    def num_rows(self) -> int:
        return len(self.tableau)

    def num_cols(self) -> int:
        return len(self.var_names)

    def print_tableau(self, title: str = "", note: str = "", obj_type: str = "min"):
        """Prints the tableau in a beautifully formatted ASCII table with RHS on the left and no z column."""
        print("\n" + "=" * 82)
        if title:
            print(f" {title.center(80)}")
            print("=" * 82)
        if note:
            print(f" Note: {note}")
            print("-" * 82)

        # Place RHS on the left: Basis | RHS | x1 | x2 | ... (No z column)
        col_headers = ["Basis", "RHS"] + self.var_names
        
        # Build cell contents as strings
        table_cells: List[List[str]] = []
        table_cells.append(col_headers)

        # Row 0 (Cost / Objective row)
        # Keeps exact raw RHS: -z = raw_rhs (or -w = raw_rhs)
        raw_obj_rhs = self.tableau[0][-1]
        row0_cells = [self.obj_name, format_fraction(raw_obj_rhs)]
        for j in range(self.num_cols()):
            row0_cells.append(format_fraction(self.tableau[0][j]))
        table_cells.append(row0_cells)

        # Rows 1 to m (Constraint rows)
        for i in range(1, self.num_rows()):
            b_var = self.basic_vars[i - 1]
            row_cells = [b_var, format_fraction(self.tableau[i][-1])]
            for j in range(self.num_cols()):
                row_cells.append(format_fraction(self.tableau[i][j]))
            table_cells.append(row_cells)

        # Calculate optimal column widths (based on unformatted text length)
        num_columns = len(col_headers)
        col_widths = []
        for c in range(num_columns):
            max_len = max(len(table_cells[r][c]) for r in range(len(table_cells)))
            col_widths.append(max(max_len, 5))

        sep_line = "+" + "+".join("-" * (w + 2) for w in col_widths) + "+"

        # Detect terminal color support
        use_color = (
            os.environ.get("NO_COLOR") is None
            and hasattr(sys.stdout, "isatty")
            and sys.stdout.isatty()
        )

        def format_cell_colored(raw_text: str, r: int, c: int, width: int) -> str:
            # First, pad the text with exact plain width
            if c == 0:
                padded = raw_text.center(width)
            else:
                padded = raw_text.rjust(width) if r > 0 else raw_text.center(width)

            if not use_color:
                return padded

            # Apply ANSI colors based on rules:
            # c == 0: Basis
            # c == 1: RHS (Violet/Magenta: \033[1;35m)
            # c >= 2: variable columns
            if c == 0:
                if r == 0:
                    return f"\033[1m{padded}\033[0m"
                elif r == 1:
                    return f"\033[1;33m{padded}\033[0m"  # Objective basis (-z, -w)
                else:
                    b_var = table_cells[r][0].strip()
                    if b_var.lower().startswith("y"):
                        return f"\033[1;33m{padded}\033[0m"  # Yellow for y artificial basis
                    else:
                        return f"\033[1;31m{padded}\033[0m"  # Basic variable name in red
            elif c == 1:
                # RHS: PURPLE / MAGENTA
                return f"\033[1;35m{padded}\033[0m"
            else:
                var_idx = c - 2
                var_name = self.var_names[var_idx]
                is_y = var_name.lower().startswith("y")
                is_basic = var_name in self.basic_vars
                if is_y:
                    if r == 0:  # Header
                        return f"\033[1;33m{padded}\033[0m"  # Yellow bold header
                    elif r == 1:  # Row 0 (Reduced Costs)
                        return f"\033[1;33m{padded}\033[0m"  # Yellow reduced cost
                    else:  # Constraint rows
                        return f"\033[33m{padded}\033[0m"    # Yellow cell for y
                else:
                    if r == 0:  # Header
                        if is_basic:
                            return f"\033[1;31m{padded}\033[0m"  # RED for basic column
                        else:
                            return f"\033[1;34m{padded}\033[0m"  # BLUE for non-basic column
                    elif r == 1:  # Row 0 (Reduced Costs)
                        return f"\033[1;32m{padded}\033[0m"      # GREEN for reduced costs!
                    else:  # Constraint rows
                        if is_basic:
                            return f"\033[31m{padded}\033[0m"    # RED for basic column
                        else:
                            return f"\033[34m{padded}\033[0m"    # BLUE for non-basic column

        print(sep_line)
        header_str = "| " + " | ".join(format_cell_colored(table_cells[0][c], 0, c, col_widths[c]) for c in range(num_columns)) + " |"
        print(header_str)
        print(sep_line)

        for r in range(1, len(table_cells)):
            row_str = "| " + " | ".join(
                format_cell_colored(table_cells[r][c], r, c, col_widths[c])
                for c in range(num_columns)
            ) + " |"
            print(row_str)
            if r == 1:
                # Separator after objective row
                print(sep_line)
        print(sep_line)
        if use_color:
            print("  Color Scheme: \033[31m[Red: Basic Cols]\033[0m \033[34m[Blue: Non-Basic Cols]\033[0m \033[33m[Yellow: y Vars]\033[0m \033[32m[Green: Reduced Costs (Row 0)]\033[0m \033[35m[Magenta: RHS]\033[0m")

        # Print current basic solution
        print("\n  Current Basic Solution:")
        if self.obj_name.startswith("-") or obj_type == "min":
            clean_name = self.obj_name.lstrip("-")
            z_val = -raw_obj_rhs
            print(f"    {self.obj_name} = {format_approx(raw_obj_rhs)}   ==>   {clean_name} = {format_approx(z_val)}")
        else:
            print(f"    {self.obj_name} = {format_approx(raw_obj_rhs)}")
        basic_sol_parts = []
        for i in range(1, self.num_rows()):
            b_var = self.basic_vars[i - 1]
            val = self.tableau[i][-1]
            basic_sol_parts.append(f"{b_var} = {format_fraction(val)}")
        print("    " + ", ".join(basic_sol_parts))
        non_basics = [v for v in self.var_names if v not in self.basic_vars]
        if non_basics:
            print("    Non-basic (0): " + ", ".join(non_basics))
        print()


def solve_simplex(tableau: SimplexTableau, phase_name: str = "",
                  forbidden_vars: Optional[set] = None,
                  obj_type: str = "min") -> Tuple[str, SimplexTableau]:
    """
    Executes the Simplex Method with Bland's Rule on the given tableau.
    
    Returns:
        (status, final_tableau)
        where status is 'OPTIMAL', 'UNBOUNDED', or 'ITERATION_LIMIT'.
    """
    if forbidden_vars is None:
        forbidden_vars = set()
    iteration = 0
    max_iterations = 200

    while iteration < max_iterations:
        title = f"{phase_name} - Tableau {iteration}" if phase_name else f"Tableau {iteration}"
        if iteration == 0:
            title += " (Initial Tableau)"
        
        tableau.print_tableau(title=title, obj_type=obj_type)

        # -------------------------------------------------------------
        # STEP 1: BLAND'S RULE - Entering Variable Selection
        # -------------------------------------------------------------
        # In Row 0 (objective row: z - sum c_j x_j = RHS),
        # an entering variable is a non-basic candidate with coefficient < 0.
        candidates: List[Tuple[int, str, int, Fraction]] = []
        # (var_index, var_name, col_index, coeff)
        for col_idx, var_name in enumerate(tableau.var_names):
            if var_name in forbidden_vars:
                continue
            coeff = tableau.tableau[0][col_idx]
            if coeff < 0:
                v_idx = tableau.var_indices[var_name]
                candidates.append((v_idx, var_name, col_idx, coeff))

        print("-" * 82)
        print(" [STEP 1] Optimality Check & Entering Variable Selection (Bland's Rule)")
        print("-" * 82)

        if not candidates:
            print(" • All coefficients in Row 0 are >= 0.")
            print(" • OPTIMALITY CRITERION SATISFIED! Current tableau is optimal.")
            return "OPTIMAL", tableau

        print(f" • Found {len(candidates)} candidate(s) with strictly negative Row 0 coefficient (improving):")
        for v_idx, v_name, col_idx, coeff in candidates:
            print(f"    - Variable {v_name:<5} | Bland Index: {v_idx:<3} | Row 0 Coeff: {format_fraction(coeff)}")

        # Bland's Rule: Choose the variable with the SMALLEST index
        candidates.sort(key=lambda item: item[0])
        entering_idx, entering_name, entering_col, entering_coeff = candidates[0]

        print("\n • BLAND'S RULE SELECTION:")
        print(f"   Entering variable chosen: '{entering_name}' (Bland Index {entering_idx})")
        if len(candidates) > 1:
            print(f"   (Selected '{entering_name}' because its index {entering_idx} is the smallest among all candidates)")

        # -------------------------------------------------------------
        # STEP 2: BLAND'S RULE - Leaving Variable Selection (Ratio Test)
        # -------------------------------------------------------------
        print("\n" + "-" * 82)
        print(f" [STEP 2] Minimum Ratio Test for Column '{entering_name}' (Leaving Variable)")
        print("-" * 82)

        ratios: List[Tuple[Fraction, int, str, int]] = []
        # (ratio, basic_var_index, basic_var_name, row_idx)

        for row_idx in range(1, tableau.num_rows()):
            b_var = tableau.basic_vars[row_idx - 1]
            b_idx = tableau.var_indices[b_var]
            a_ij = tableau.tableau[row_idx][entering_col]
            rhs = tableau.tableau[row_idx][-1]

            if a_ij > 0:
                ratio = rhs / a_ij
                ratios.append((ratio, b_idx, b_var, row_idx))
                print(f"   Row {row_idx} (Basis: {b_var:<5}, Index: {b_idx:<3}): "
                      f"RHS / Coeff = ({format_fraction(rhs)}) / ({format_fraction(a_ij)}) = {format_fraction(ratio)}")
            else:
                reason = "denominator = 0" if a_ij == 0 else f"denominator < 0 ({format_fraction(a_ij)})"
                print(f"   Row {row_idx} (Basis: {b_var:<5}, Index: {b_idx:<3}): "
                      f"No valid ratio ({reason}) -> Skipped")

        if not ratios:
            print("\n [!] UNBOUNDED PROBLEM DETECTED:")
            print(f"     All constraint coefficients in column '{entering_name}' are <= 0.")
            print(f"     Variable '{entering_name}' can be increased indefinitely, driving {tableau.obj_name} to infinity.")
            return "UNBOUNDED", tableau

        min_ratio = min(item[0] for item in ratios)
        tied_candidates = [item for item in ratios if item[0] == min_ratio]

        print(f"\n • Minimum ratio: {format_fraction(min_ratio)}")

        if len(tied_candidates) > 1:
            print(f" • DEGENERACY / TIE DETECTED between {len(tied_candidates)} rows having ratio {format_fraction(min_ratio)}:")
            for r_val, b_idx, b_var, r_idx in tied_candidates:
                print(f"    - Row {r_idx}: Basic variable {b_var} (Bland Index {b_idx})")
            
            # Bland's Rule: Choose the tied basic variable with the SMALLEST index
            tied_candidates.sort(key=lambda item: item[1])
            leaving_ratio, leaving_b_idx, leaving_name, pivot_row = tied_candidates[0]
            print("\n • BLAND'S TIE-BREAKING RULE APPLIED:")
            print(f"   Tied basic variables: {[item[2] for item in tied_candidates]}")
            print(f"   Smallest basic variable index is {leaving_b_idx} for '{leaving_name}'.")
            print(f"   Leaving variable chosen: '{leaving_name}' (Row {pivot_row})")
        else:
            leaving_ratio, leaving_b_idx, leaving_name, pivot_row = tied_candidates[0]
            print(f" • Unique minimum ratio achieved by Row {pivot_row} (Basis: '{leaving_name}', Index: {leaving_b_idx}).")
            print(f"   Leaving variable chosen: '{leaving_name}' (Row {pivot_row})")

        pivot_element = tableau.tableau[pivot_row][entering_col]
        print(f"\n • PIVOT DETAILS:")
        print(f"   Pivot Row:    Row {pivot_row} (Leaving: '{leaving_name}')")
        print(f"   Pivot Column: Column '{entering_name}' (Index: {entering_col + 1})")
        print(f"   Pivot Value:  {format_fraction(pivot_element)}")

        # -------------------------------------------------------------
        # STEP 3: ELEMENTARY ROW OPERATIONS (Gaussian Pivot)
        # -------------------------------------------------------------
        print("\n" + "-" * 82)
        print(" [STEP 3] Performing Elementary Row Operations (Pivot)")
        print("-" * 82)

        # 1. Normalize pivot row
        tableau.tableau[pivot_row] = [x / pivot_element for x in tableau.tableau[pivot_row]]
        print(f" • R{pivot_row} -> R{pivot_row} / ({format_fraction(pivot_element)})")

        # 2. Eliminate other rows
        for r in range(tableau.num_rows()):
            if r != pivot_row:
                factor = tableau.tableau[r][entering_col]
                if factor != 0:
                    tableau.tableau[r] = [
                        tableau.tableau[r][c] - factor * tableau.tableau[pivot_row][c]
                        for c in range(len(tableau.tableau[r]))
                    ]
                    row_name = tableau.obj_name if r == 0 else f"R{r}"
                    sign = "-" if factor > 0 else "+"
                    abs_factor = abs(factor)
                    print(f" • {row_name:<2} -> {row_name:<2} {sign} ({format_fraction(abs_factor)}) * R{pivot_row}")

        # Update basis
        tableau.basic_vars[pivot_row - 1] = entering_name
        print(f" • Updated Basis for Row {pivot_row}: '{leaving_name}' -> '{entering_name}'")

        iteration += 1

    print("\n [!] Reached maximum iteration limit (200). Stopping.")
    return "ITERATION_LIMIT", tableau

def make_step_snapshot(
    tab: SimplexTableau,
    iteration: int,
    substep_idx: int,
    substep_total: int,
    substep_type: str,
    substep_title: str,
    substep_desc: str = "",
    formula: str = "",
    target_row: Optional[int] = None,
    pivot_row: Optional[int] = None,
    pivot_col: Optional[int] = None,
    eliminated_cell: Optional[Dict] = None,
    active_rows: Optional[List[int]] = None,
    active_cols: Optional[List[int]] = None,
    entering_var: Optional[Dict] = None,
    leaving_var: Optional[Dict] = None,
    candidates: Optional[List[Dict]] = None,
    ratios: Optional[List[Dict]] = None,
    min_ratio: Optional[str] = None,
    tied: bool = False,
    bland_tie_note: Optional[str] = None,
    status: str = "CONTINUE",
    message: str = "",
    phase_name: str = "",
    obj_type: str = "min"
) -> Dict:
    raw_obj_rhs = tab.tableau[0][-1]
    display_cost = get_display_obj_val(raw_obj_rhs, obj_type, tab.obj_name)

    rows = []
    basic_sol = {}
    for r in range(tab.num_rows()):
        b_var = tab.obj_name if r == 0 else tab.basic_vars[r - 1]
        cells = [format_fraction(val) for val in tab.tableau[r][:-1]]
        rhs_val = tab.tableau[r][-1]
        rhs_str = format_fraction(rhs_val)
        rows.append({
            "basic": b_var,
            "rhs": rhs_str,
            "cells": cells
        })
        if r > 0:
            basic_sol[b_var] = rhs_str
        else:
            basic_sol[tab.obj_name] = rhs_str

    pivot_dict = None
    if (pivot_row is not None and pivot_col is not None
            and 0 <= pivot_row < tab.num_rows()
            and 0 <= pivot_col < len(tab.var_names)):
        pivot_dict = {
            "row": pivot_row,
            "col": pivot_col,
            "val": format_fraction(tab.tableau[pivot_row][pivot_col])
        }

    return {
        "iteration": iteration,
        "substep_idx": substep_idx,
        "substep_total": substep_total,
        "substep_type": substep_type,
        "title": substep_title,
        "substep_desc": substep_desc,
        "formula": formula,
        "target_row": target_row,
        "pivot_row": pivot_row,
        "pivot_col": pivot_col,
        "eliminated_cell": eliminated_cell,
        "active_rows": active_rows if active_rows is not None else [],
        "active_cols": active_cols if active_cols is not None else [],
        "phase_name": phase_name,
        "is_phase1": (phase_name == "Phase 1" or "w" in tab.obj_name),
        "is_phase2": (phase_name == "Phase 2"),
        "var_names": list(tab.var_names),
        "var_indices": dict(tab.var_indices),
        "basic_vars": list(tab.basic_vars),
        "obj_name": tab.obj_name,
        "obj_type": obj_type,
        "rows": rows,
        "basic_solution": basic_sol,
        "raw_rhs": format_fraction(raw_obj_rhs),
        "obj_val": format_approx(display_cost),
        "raw_obj_val": format_fraction(display_cost),
        "status": status,
        "message": message,
        "pivot": pivot_dict,
        "entering_var": entering_var,
        "leaving_var": leaving_var,
        "candidates": candidates or [],
        "ratios": ratios or [],
        "min_ratio": min_ratio,
        "tied": tied,
        "bland_tie_note": bland_tie_note
    }


def record_simplex_steps(
    initial_tableau: SimplexTableau,
    phase_name: str = "",
    forbidden_vars: Optional[Set[str]] = None,
    obj_type: str = "min"
) -> Tuple[str, List[Dict], SimplexTableau]:
    """
    Solves a tableau using Bland's rule and records step-by-step sub-operations:
    1. Pivot Selection (Bland candidate selection + Minimum Ratio theta test)
    2. Normalize Pivot Row (Make Pivot Element = 1)
    3. Row-by-Row Gauss-Jordan Elimination (Drive column entries to 0)
    4. Iteration Complete (Basis updated)
    """
    if forbidden_vars is None:
        forbidden_vars = set()

    tab = SimplexTableau(
        initial_tableau.var_names,
        initial_tableau.var_indices,
        initial_tableau.tableau,
        initial_tableau.basic_vars,
        initial_tableau.obj_name
    )

    steps: List[Dict] = []
    iteration = 0
    max_iterations = 200

    while iteration < max_iterations:
        # 1. Candidate entering variables
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
                    "coeff": format_fraction(coeff)
                })

        # Optimality criterion check
        if not candidates:
            opt_title = f"{phase_name} - Optimal Solution Reached" if phase_name else "Optimal Solution Reached"
            opt_step = make_step_snapshot(
                tab,
                iteration=iteration,
                substep_idx=1,
                substep_total=1,
                substep_type="optimal",
                substep_title=opt_title,
                substep_desc="All relative costs in Row 0 are >= 0. Optimality criterion satisfied!",
                formula="Optimal: All c̄_j >= 0",
                status="OPTIMAL",
                message="All coefficients in Row 0 are >= 0. Optimality criterion satisfied!",
                phase_name=phase_name,
                obj_type=obj_type
            )
            steps.append(opt_step)
            return "OPTIMAL", steps, tab

        # Bland's Rule: Choose candidate with smallest index
        candidates.sort(key=lambda c: c["index"])
        entering = candidates[0]
        entering_col = entering["col"]

        # Minimum Ratio test theta (with bracketed fractions)
        ratios = []
        for row_idx in range(1, tab.num_rows()):
            b_var = tab.basic_vars[row_idx - 1]
            b_idx = tab.var_indices[b_var]
            a_ij = tab.tableau[row_idx][entering_col]
            rhs = tab.tableau[row_idx][-1]
            calc_str = f"[{format_fraction(rhs)}] / [{format_fraction(a_ij)}]"

            if a_ij > 0:
                ratio = rhs / a_ij
                ratios.append({
                    "row": row_idx,
                    "basic": b_var,
                    "index": b_idx,
                    "ratio": format_fraction(ratio),
                    "theta": format_fraction(ratio),
                    "ratio_val": ratio,
                    "valid": True,
                    "rhs": format_fraction(rhs),
                    "coeff": format_fraction(a_ij),
                    "calc_str": calc_str
                })
            else:
                reason = f"coeff <= 0 ([{format_fraction(a_ij)}])"
                ratios.append({
                    "row": row_idx,
                    "basic": b_var,
                    "index": b_idx,
                    "ratio": None,
                    "theta": None,
                    "valid": False,
                    "reason": reason,
                    "rhs": format_fraction(rhs),
                    "coeff": format_fraction(a_ij),
                    "calc_str": calc_str
                })

        valid_ratios = [r for r in ratios if r["valid"]]
        if not valid_ratios:
            unb_title = f"{phase_name} - Problem Unbounded" if phase_name else "Problem Unbounded"
            unb_step = make_step_snapshot(
                tab,
                iteration=iteration,
                substep_idx=1,
                substep_total=1,
                substep_type="unbounded",
                substep_title=unb_title,
                substep_desc=f"All coefficients in column '{entering['name']}' are <= 0. The problem is UNBOUNDED!",
                formula=f"Unbounded: all a_i,{entering['name']} <= 0",
                status="UNBOUNDED",
                message=f"All coefficients in column '{entering['name']}' are <= 0. The problem is UNBOUNDED!",
                entering_var=entering,
                candidates=candidates,
                ratios=ratios,
                phase_name=phase_name,
                obj_type=obj_type
            )
            steps.append(unb_step)
            return "UNBOUNDED", steps, tab

        min_ratio_val = min(r["ratio_val"] for r in valid_ratios)
        tied_candidates = [r for r in valid_ratios if r["ratio_val"] == min_ratio_val]
        tied = len(tied_candidates) > 1

        if tied:
            tied_candidates.sort(key=lambda r: r["index"])
            leaving_item = tied_candidates[0]
            bland_tie_note = (
                f"Tie detected between {len(tied_candidates)} rows with minimum θ = {leaving_item['theta']}. "
                f"Bland's Rule selects basic variable with smallest index: '{leaving_item['basic']}' (Index {leaving_item['index']})."
            )
        else:
            leaving_item = tied_candidates[0]
            bland_tie_note = None

        pivot_row = leaving_item["row"]
        pivot_element = tab.tableau[pivot_row][entering_col]

        for r in ratios:
            r.pop("ratio_val", None)
        leaving = dict(leaving_item)
        leaving.pop("ratio_val", None)

        rows_to_eliminate = [r for r in range(tab.num_rows()) if r != pivot_row and tab.tableau[r][entering_col] != 0]
        substep_total = 1 + 1 + len(rows_to_eliminate) + 1

        # SUB-STEP 1: Pivot Selection (Bland's Rule & Minimum theta Test)
        substep_1 = make_step_snapshot(
            tab,
            iteration=iteration,
            substep_idx=1,
            substep_total=substep_total,
            substep_type="pivot_selection",
            substep_title=f"{phase_name} Iteration {iteration} • Step 1: Select Pivot" if phase_name else f"Iteration {iteration} • Step 1: Select Pivot",
            substep_desc=f"Entering: '{entering['name']}' (smallest index with negative reduced cost). Leaving: '{leaving['basic']}' (minimum θ = {format_fraction(min_ratio_val)}). Pivot is at (Row R{pivot_row}, {entering['name']}) = {format_fraction(pivot_element)}.",
            formula=f"Pivot: a(R{pivot_row}, {entering['name']}) = {format_fraction(pivot_element)}, min θ = {format_fraction(min_ratio_val)}",
            pivot_row=pivot_row,
            pivot_col=entering_col,
            active_rows=[pivot_row],
            active_cols=[entering_col],
            entering_var=entering,
            leaving_var=leaving,
            candidates=candidates,
            ratios=ratios,
            min_ratio=format_fraction(min_ratio_val),
            tied=tied,
            bland_tie_note=bland_tie_note,
            phase_name=phase_name,
            obj_type=obj_type
        )
        steps.append(substep_1)

        # SUB-STEP 2: Normalize Pivot Row
        tab.tableau[pivot_row] = [x / pivot_element for x in tab.tableau[pivot_row]]
        substep_2 = make_step_snapshot(
            tab,
            iteration=iteration,
            substep_idx=2,
            substep_total=substep_total,
            substep_type="normalize_pivot",
            substep_title=f"{phase_name} Iteration {iteration} • Step 2: Normalize Pivot Row R{pivot_row}" if phase_name else f"Iteration {iteration} • Step 2: Normalize Pivot Row R{pivot_row}",
            substep_desc=f"Dividing pivot row R{pivot_row} by pivot element ({format_fraction(pivot_element)}) so that pivot entry becomes 1.",
            formula=f"R{pivot_row} ← R{pivot_row} / ({format_fraction(pivot_element)})",
            target_row=pivot_row,
            pivot_row=pivot_row,
            pivot_col=entering_col,
            active_rows=[pivot_row],
            active_cols=[entering_col],
            entering_var=entering,
            leaving_var=leaving,
            candidates=candidates,
            ratios=ratios,
            min_ratio=format_fraction(min_ratio_val),
            phase_name=phase_name,
            obj_type=obj_type
        )
        steps.append(substep_2)

        # SUB-STEPS 3+: Gauss-Jordan Elimination Row by Row
        curr_substep = 3
        for r in rows_to_eliminate:
            factor = tab.tableau[r][entering_col]
            tab.tableau[r] = [
                tab.tableau[r][c] - factor * tab.tableau[pivot_row][c]
                for c in range(len(tab.tableau[r]))
            ]

            r_name = tab.obj_name if r == 0 else f"R{r}"
            sign = "-" if factor > 0 else "+"
            abs_factor_str = format_fraction(abs(factor))
            op_formula = f"{r_name} ← {r_name} {sign} ({abs_factor_str}) * R{pivot_row}" if abs_factor_str != "1" else f"{r_name} ← {r_name} {sign} R{pivot_row}"

            substep_elim = make_step_snapshot(
                tab,
                iteration=iteration,
                substep_idx=curr_substep,
                substep_total=substep_total,
                substep_type="eliminate_row",
                substep_title=f"{phase_name} Iteration {iteration} • Step {curr_substep}: Eliminate in {r_name}" if phase_name else f"Iteration {iteration} • Step {curr_substep}: Eliminate in {r_name}",
                substep_desc=f"Subtracting ({format_fraction(factor)}) * R{pivot_row} from {r_name} to drive the coefficient of '{entering['name']}' in {r_name} to 0.",
                formula=op_formula,
                target_row=r,
                pivot_row=pivot_row,
                pivot_col=entering_col,
                eliminated_cell={"row": r, "col": entering_col},
                active_rows=[r, pivot_row],
                active_cols=[entering_col],
                entering_var=entering,
                leaving_var=leaving,
                candidates=candidates,
                ratios=ratios,
                min_ratio=format_fraction(min_ratio_val),
                phase_name=phase_name,
                obj_type=obj_type
            )
            steps.append(substep_elim)
            curr_substep += 1

        # SUB-STEP FINAL: Iteration Complete
        tab.basic_vars[pivot_row - 1] = entering["name"]
        substep_complete = make_step_snapshot(
            tab,
            iteration=iteration,
            substep_idx=curr_substep,
            substep_total=substep_total,
            substep_type="iteration_complete",
            substep_title=f"{phase_name} Iteration {iteration} Complete: Tableau Updated" if phase_name else f"Iteration {iteration} Complete: Tableau Updated",
            substep_desc=f"Variable '{entering['name']}' has entered the basis, replacing '{leaving['basic']}'. Column '{entering['name']}' is now a unit vector.",
            formula=f"Basis Updated: R{pivot_row} is now '{entering['name']}' (replacing '{leaving['basic']}')",
            target_row=pivot_row,
            pivot_row=pivot_row,
            pivot_col=entering_col,
            active_rows=[],
            active_cols=[],
            entering_var=entering,
            leaving_var=leaving,
            candidates=candidates,
            ratios=ratios,
            min_ratio=format_fraction(min_ratio_val),
            phase_name=phase_name,
            obj_type=obj_type
        )
        steps.append(substep_complete)

        iteration += 1

    return "ITERATION_LIMIT", steps, tab


def record_tableau_steps(tab: SimplexTableau, obj_type: str = "min") -> Dict:
    """Solves a tableau and returns complete trajectory and final solution data for GUI."""
    status, steps, final_tab = record_simplex_steps(tab, obj_type=obj_type)
    result = {
        "status": status,
        "obj_type": obj_type,
        "steps": steps,
        "final_solution": None
    }
    if status == "OPTIMAL":
        raw_z = final_tab.tableau[0][-1]
        opt_z = get_display_obj_val(raw_z, obj_type, final_tab.obj_name)
        sol_dict = {name: format_fraction(Fraction(0)) for name in final_tab.var_names}
        for r_idx in range(1, final_tab.num_rows()):
            b_var = final_tab.basic_vars[r_idx - 1]
            sol_dict[b_var] = format_fraction(final_tab.tableau[r_idx][-1])

        alt_optima = [
            name for c_idx, name in enumerate(final_tab.var_names)
            if name not in final_tab.basic_vars and final_tab.tableau[0][c_idx] == 0
        ]

        result["final_solution"] = {
            "obj_val": format_approx(opt_z),
            "raw_obj_val": format_fraction(opt_z),
            "raw_rhs": format_fraction(raw_z),
            "obj_name": final_tab.obj_name,
            "solution_vector": sol_dict,
            "basic_vars": list(final_tab.basic_vars),
            "alternative_optima": alt_optima
        }
    return result


def record_tableau_steps_twophase(tab: SimplexTableau, obj_type: str = "min") -> Dict:
    """
    Forces Two-Phase Simplex on an imported Simplex Tableau:
    1. Ensures all constraint rows have RHS >= 0 (multiplying negative RHS rows by -1).
    2. Phase 1: Introduces artificial variables y1, ..., ym as the starting basis.
       Constructs Phase 1 objective Row 0 (-w = -sum(y_i)) and canonicalizes.
       Solves Phase 1 using Bland's rule until w = 0 (feasible) or min w > 0 (infeasible).
    3. Phase 2: Drops artificial columns, restores original objective Row 0, canonicalizes
       with respect to the Phase 1 final basis, and solves to optimality.
    """
    m = tab.num_rows() - 1
    if m <= 0:
        return record_tableau_steps(tab, obj_type=obj_type)

    # If the tableau is already an active Phase 1 tableau (-w or contains artificial variables)
    if tab.obj_name == "-w" or any(v.startswith("y") for v in tab.var_names):
        return record_tableau_steps(tab, obj_type=obj_type)

    original_var_names = list(tab.var_names)
    n_orig = len(original_var_names)

    # 1. Normalize constraint rows so all RHS >= 0
    norm_matrix_A: List[List[Fraction]] = []
    for r in range(1, tab.num_rows()):
        row = list(tab.tableau[r])
        if row[-1] < 0:
            row = [-val for val in row]
        norm_matrix_A.append(row)

    # Save original objective row
    orig_row0 = list(tab.tableau[0])
    orig_obj_name = tab.obj_name
    orig_obj_type = obj_type

    # 2. Phase 1 Setup
    artificials = [f"y{i+1}" for i in range(m)]
    phase1_var_names = original_var_names + artificials
    phase1_var_indices = {name: idx + 1 for idx, name in enumerate(phase1_var_names)}
    n_full = len(phase1_var_names)

    # Build Phase 1 constraint rows (original coeffs + identity block for y + RHS)
    phase1_rows: List[List[Fraction]] = []
    for i in range(m):
        row = [Fraction(0)] * (n_full + 1)
        for j in range(n_orig):
            row[j] = norm_matrix_A[i][j]
        row[n_orig + i] = Fraction(1)
        row[-1] = norm_matrix_A[i][-1]
        phase1_rows.append(row)

    # Build Phase 1 Row 0: -w = -sum(y_i)
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

    combined_steps: List[Dict] = []
    p1_status, p1_steps, solved_phase1_tab = record_simplex_steps(
        phase1_tab,
        phase_name="Phase 1",
        forbidden_vars=set(artificials),
        obj_type="min"
    )
    combined_steps.extend(p1_steps)

    if p1_status != "OPTIMAL":
        return {
            "status": p1_status,
            "obj_type": orig_obj_type,
            "steps": combined_steps,
            "final_solution": None
        }

    w_opt = solved_phase1_tab.tableau[0][-1]
    if w_opt != 0:
        return {
            "status": "INFEASIBLE",
            "obj_type": orig_obj_type,
            "steps": combined_steps,
            "message": f"The linear program is INFEASIBLE (Minimum sum of artificial variables = {format_fraction(-w_opt)} > 0)",
            "final_solution": None
        }

    # Degeneracy check: drive any remaining artificial variables out of the basis
    redundant_rows = []
    for r_idx in range(1, solved_phase1_tab.num_rows()):
        b_var = solved_phase1_tab.basic_vars[r_idx - 1]
        if b_var in artificials:
            driven_out = False
            for c_idx, v_name in enumerate(solved_phase1_tab.var_names):
                if v_name not in artificials and solved_phase1_tab.tableau[r_idx][c_idx] != 0:
                    p_val = solved_phase1_tab.tableau[r_idx][c_idx]
                    solved_phase1_tab.tableau[r_idx] = [x / p_val for x in solved_phase1_tab.tableau[r_idx]]
                    for r_other in range(solved_phase1_tab.num_rows()):
                        if r_other != r_idx:
                            factor = solved_phase1_tab.tableau[r_other][c_idx]
                            solved_phase1_tab.tableau[r_other] = [
                                solved_phase1_tab.tableau[r_other][c] - factor * solved_phase1_tab.tableau[r_idx][c]
                                for c in range(len(solved_phase1_tab.tableau[r_other]))
                            ]
                    solved_phase1_tab.basic_vars[r_idx - 1] = v_name
                    driven_out = True
                    break
            if not driven_out:
                redundant_rows.append(r_idx)

    # Strip artificial columns for Phase 2
    non_art_cols = [c for c, name in enumerate(solved_phase1_tab.var_names) if name not in artificials]
    phase2_var_names = [solved_phase1_tab.var_names[c] for c in non_art_cols]
    phase2_var_indices = {name: idx + 1 for idx, name in enumerate(phase2_var_names)}

    phase2_rows = []
    phase2_basic_vars = []
    for r in range(1, solved_phase1_tab.num_rows()):
        if r in redundant_rows:
            continue
        row = [solved_phase1_tab.tableau[r][c] for c in non_art_cols] + [solved_phase1_tab.tableau[r][-1]]
        phase2_rows.append(row)
        phase2_basic_vars.append(solved_phase1_tab.basic_vars[r - 1])

    # 3. Phase 2 Objective Setup
    phase2_row0 = [Fraction(0)] * (len(phase2_var_names) + 1)
    for j, name in enumerate(phase2_var_names):
        if name in original_var_names:
            orig_col = original_var_names.index(name)
            phase2_row0[j] = orig_row0[orig_col]
    phase2_row0[-1] = orig_row0[-1]

    # Canonicalize Phase 2 Row 0 with respect to the new basis
    for r_idx, b_var in enumerate(phase2_basic_vars):
        c_idx = phase2_var_names.index(b_var)
        factor = phase2_row0[c_idx]
        if factor != 0:
            for col in range(len(phase2_row0)):
                phase2_row0[col] -= factor * phase2_rows[r_idx][col]

    phase2_tableau_data = [phase2_row0] + phase2_rows
    phase2_tab = SimplexTableau(phase2_var_names, phase2_var_indices, phase2_tableau_data,
                                phase2_basic_vars, obj_name=orig_obj_name)

    p2_status, p2_steps, final_tab = record_simplex_steps(
        phase2_tab,
        phase_name="Phase 2",
        obj_type=orig_obj_type
    )
    combined_steps.extend(p2_steps)

    result = {
        "status": p2_status,
        "obj_type": orig_obj_type,
        "needs_phase1": True,
        "steps": combined_steps,
        "final_solution": None
    }

    if p2_status == "OPTIMAL":
        raw_z = final_tab.tableau[0][-1]
        opt_z = get_display_obj_val(raw_z, orig_obj_type, final_tab.obj_name)
        sol_dict = {name: format_fraction(Fraction(0)) for name in final_tab.var_names}
        for r_idx in range(1, final_tab.num_rows()):
            b_var = final_tab.basic_vars[r_idx - 1]
            sol_dict[b_var] = format_fraction(final_tab.tableau[r_idx][-1])

        alt_optima = [
            name for c_idx, name in enumerate(final_tab.var_names)
            if name not in final_tab.basic_vars and final_tab.tableau[0][c_idx] == 0
        ]

        result["final_solution"] = {
            "obj_val": format_approx(opt_z),
            "raw_obj_val": format_fraction(opt_z),
            "raw_rhs": format_fraction(raw_z),
            "obj_name": final_tab.obj_name,
            "solution_vector": sol_dict,
            "basic_vars": list(final_tab.basic_vars),
            "alternative_optima": alt_optima
        }

    return result


def record_lp_steps(lp: 'LinearProgram') -> Dict:
    """Solves a Linear Program (with Two-Phase if needed) and returns complete trajectory for GUI."""
    norm_constraints = []
    for coeffs, rel, rhs in lp.constraints:
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

    needs_phase1 = getattr(lp, "force_two_phase", False) or any(rel in (">=", "=") for _, rel, _ in norm_constraints)

    all_var_names: List[str] = list(lp.var_names)
    slacks: List[str] = []
    surpluses: List[str] = []
    artificials: List[str] = []

    # Slack variables named as x, continuing from highest existing decision variable index
    existing_x_indices = [
        int(v[1:]) for v in lp.var_names if v.startswith("x") and v[1:].isdigit()
    ]
    next_x_idx = max(existing_x_indices, default=len(lp.var_names)) + 1
    surplus_counter = 1
    artificial_counter = 1

    constraint_aux: List[Dict[str, Optional[str]]] = []

    for coeffs, rel, rhs in norm_constraints:
        aux: Dict[str, Optional[str]] = {"slack": None, "surplus": None, "artificial": None}
        if rel == "<=":
            s_name = f"x{next_x_idx}"
            next_x_idx += 1
            slacks.append(s_name)
            aux["slack"] = s_name
            # In Two-Phase, force the use of y as the basis for ALL constraints,
            # even if we already have a basis with slack variables!
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

    full_var_names = lp.var_names + slacks + surpluses + artificials
    var_indices = {name: idx + 1 for idx, name in enumerate(full_var_names)}

    m = len(norm_constraints)
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

    combined_steps: List[Dict] = []

    if needs_phase1:
        phase1_row0 = [Fraction(0)] * (n_full + 1)
        for a_name in artificials:
            phase1_row0[full_var_names.index(a_name)] = Fraction(1)
        for i, aux in enumerate(constraint_aux):
            if aux["artificial"]:
                for c_idx in range(n_full + 1):
                    phase1_row0[c_idx] -= matrix_A[i][c_idx]

        phase1_tableau_data = [phase1_row0] + [list(r) for r in matrix_A]
        phase1_tab = SimplexTableau(full_var_names, var_indices, phase1_tableau_data, basic_vars, obj_name="-w")

        p1_status, p1_steps, solved_phase1_tab = record_simplex_steps(phase1_tab, phase_name="Phase 1",
                                                                       forbidden_vars=set(artificials),
                                                                       obj_type="min")
        combined_steps.extend(p1_steps)

        if p1_status != "OPTIMAL":
            return {"status": p1_status, "obj_type": lp.obj_type, "steps": combined_steps, "final_solution": None}

        w_opt = solved_phase1_tab.tableau[0][-1]
        if w_opt != 0:
            return {
                "status": "INFEASIBLE",
                "obj_type": lp.obj_type,
                "steps": combined_steps,
                "message": f"The linear program is INFEASIBLE (Sum of artificial variables = {format_fraction(-w_opt)} > 0)",
                "final_solution": None
            }

        # Remove artificials from basis if degenerate
        for r_idx in range(1, solved_phase1_tab.num_rows()):
            b_var = solved_phase1_tab.basic_vars[r_idx - 1]
            if b_var in artificials:
                for c_idx, v_name in enumerate(solved_phase1_tab.var_names):
                    if v_name not in artificials and solved_phase1_tab.tableau[r_idx][c_idx] != 0:
                        p_val = solved_phase1_tab.tableau[r_idx][c_idx]
                        solved_phase1_tab.tableau[r_idx] = [x / p_val for x in solved_phase1_tab.tableau[r_idx]]
                        for r_other in range(solved_phase1_tab.num_rows()):
                            if r_other != r_idx:
                                factor = solved_phase1_tab.tableau[r_other][c_idx]
                                solved_phase1_tab.tableau[r_other] = [
                                    solved_phase1_tab.tableau[r_other][c] - factor * solved_phase1_tab.tableau[r_idx][c]
                                    for c in range(len(solved_phase1_tab.tableau[r_other]))
                                ]
                        solved_phase1_tab.basic_vars[r_idx - 1] = v_name
                        break

        # Remove artificial columns for Phase 2
        non_art_cols = [c for c, name in enumerate(solved_phase1_tab.var_names) if name not in artificials]
        phase2_var_names = [solved_phase1_tab.var_names[c] for c in non_art_cols]
        phase2_var_indices = {name: var_indices[name] for name in phase2_var_names}

        phase2_rows = []
        for r in range(1, solved_phase1_tab.num_rows()):
            row = [solved_phase1_tab.tableau[r][c] for c in non_art_cols] + [solved_phase1_tab.tableau[r][-1]]
            phase2_rows.append(row)

        phase2_basic_vars = list(solved_phase1_tab.basic_vars)

    else:
        phase2_var_names = full_var_names
        phase2_var_indices = var_indices
        phase2_rows = [list(r) for r in matrix_A]
        phase2_basic_vars = list(basic_vars)

    phase2_row0 = [Fraction(0)] * (len(phase2_var_names) + 1)
    for j, name in enumerate(lp.var_names):
        c_idx = phase2_var_names.index(name)
        if lp.obj_type == "max":
            phase2_row0[c_idx] = -lp.obj_coeffs[j]
        else:
            phase2_row0[c_idx] = lp.obj_coeffs[j]

    # Canonicalize Row 0
    for r_idx, b_var in enumerate(phase2_basic_vars):
        c_idx = phase2_var_names.index(b_var)
        factor = phase2_row0[c_idx]
        if factor != 0:
            for col in range(len(phase2_row0)):
                phase2_row0[col] -= factor * phase2_rows[r_idx][col]

    phase2_tableau_data = [phase2_row0] + phase2_rows
    p2_obj_name = "-z" if lp.obj_type == "min" else "z"
    phase2_tab = SimplexTableau(phase2_var_names, phase2_var_indices, phase2_tableau_data,
                                phase2_basic_vars, obj_name=p2_obj_name)

    p2_title = "Phase 2" if needs_phase1 else ""
    p2_status, p2_steps, final_tab = record_simplex_steps(phase2_tab, phase_name=p2_title, obj_type=lp.obj_type)
    combined_steps.extend(p2_steps)

    result = {
        "status": p2_status,
        "obj_type": lp.obj_type,
        "needs_phase1": needs_phase1,
        "steps": combined_steps,
        "final_solution": None
    }

    if p2_status == "OPTIMAL":
        raw_z = final_tab.tableau[0][-1]
        opt_z = get_display_obj_val(raw_z, lp.obj_type, final_tab.obj_name)
        sol_dict = {name: format_fraction(Fraction(0)) for name in lp.var_names}
        aux_dict = {name: format_fraction(Fraction(0)) for name in final_tab.var_names if name not in lp.var_names}

        for r_idx in range(1, final_tab.num_rows()):
            b_var = final_tab.basic_vars[r_idx - 1]
            val = final_tab.tableau[r_idx][-1]
            if b_var in sol_dict:
                sol_dict[b_var] = format_fraction(val)
            elif b_var in aux_dict:
                aux_dict[b_var] = format_fraction(val)

        alt_optima = [
            name for c_idx, name in enumerate(final_tab.var_names)
            if name not in final_tab.basic_vars and final_tab.tableau[0][c_idx] == 0
        ]

        result["final_solution"] = {
            "obj_val": format_approx(opt_z),
            "raw_obj_val": format_fraction(opt_z),
            "raw_rhs": format_fraction(raw_z),
            "obj_name": final_tab.obj_name,
            "solution_vector": sol_dict,
            "aux_vector": aux_dict,
            "basic_vars": list(final_tab.basic_vars),
            "alternative_optima": alt_optima
        }

    return result


class LinearProgram:
    """Represents a Linear Programming problem (always minimization of z)."""
    def __init__(self, obj_type: str, obj_coeffs: List[Fraction],
                 constraints: List[Tuple[List[Fraction], str, Fraction]],
                 var_names: Optional[List[str]] = None,
                 force_two_phase: bool = False):
        self.obj_type = "min"  # Always minimization of z
        self.obj_coeffs = [Fraction(c) for c in obj_coeffs]
        self.num_vars = len(self.obj_coeffs)
        self.constraints = constraints
        self.num_constraints = len(constraints)
        self.force_two_phase = force_two_phase
        
        if var_names:
            self.var_names = list(var_names)
        else:
            self.var_names = [f"x{j+1}" for j in range(self.num_vars)]

    def solve(self):
        """Sets up and solves the LP using Bland's Simplex / Two-Phase Simplex."""
        print("\n" + "#" * 82)
        print(" PROBLEM FORMULATION")
        print("#" * 82)
        
        obj_str = " + ".join(f"({format_fraction(c)})*{name}" for c, name in zip(self.obj_coeffs, self.var_names))
        print(f"  MIN {obj_str}")
        print("  Subject to:")
        for idx, (coeffs, rel, rhs) in enumerate(self.constraints):
            lhs_str = " + ".join(f"({format_fraction(c)})*{name}" for c, name in zip(coeffs, self.var_names))
            print(f"    ({idx+1}) {lhs_str} {rel} {format_fraction(rhs)}")
        print("    " + ", ".join(f"{name} >= 0" for name in self.var_names))
        print("#" * 82)

        # Step A: Standardize constraints so RHS >= 0
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

        # Check if artificial variables are needed
        needs_phase1 = self.force_two_phase or any(rel in (">=", "=") for _, rel, _ in norm_constraints)

        # Prepare variable catalog and ordering for Bland's Rule
        # Order: Decision variables (x1..xn), Slacks (s1..), Surpluses (e1..), Artificials (y1..)
        all_var_names: List[str] = list(self.var_names)
        slacks: List[str] = []
        surpluses: List[str] = []
        artificials: List[str] = []

        # Slack variables named as x, continuing from highest existing decision variable index
        existing_x_indices = [
            int(v[1:]) for v in self.var_names if v.startswith("x") and v[1:].isdigit()
        ]
        next_x_idx = max(existing_x_indices, default=len(self.var_names)) + 1
        surplus_counter = 1
        artificial_counter = 1

        constraint_aux: List[Dict[str, Optional[str]]] = []

        for coeffs, rel, rhs in norm_constraints:
            aux: Dict[str, Optional[str]] = {"slack": None, "surplus": None, "artificial": None}
            if rel == "<=":
                s_name = f"x{next_x_idx}"
                next_x_idx += 1
                slacks.append(s_name)
                aux["slack"] = s_name
                # In Two-Phase, force the use of y as the basis for ALL constraints,
                # even if we already have a basis with slack variables!
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

        print("\n Variable Ordering & Bland's Rule Indexing:")
        for name in full_var_names:
            v_type = "Decision" if name in self.var_names else ("Slack" if name in slacks else ("Surplus" if name in surpluses else "Artificial"))
            print(f"   Index {var_indices[name]:<2}: {name:<5} ({v_type})")
        print()

        # Build constraint rows
        m = len(norm_constraints)
        n_full = len(full_var_names)
        matrix_A: List[List[Fraction]] = []
        basic_vars: List[str] = []

        for i, (coeffs, rel, rhs) in enumerate(norm_constraints):
            row = [Fraction(0)] * (n_full + 1)
            # Fill decision variables
            for j, c in enumerate(coeffs):
                row[j] = c
            # Fill auxiliary variables
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

        # -------------------------------------------------------------
        # PHASE 1 (If Artificial Variables are present)
        # -------------------------------------------------------------
        if needs_phase1:
            print("\n" + "=" * 82)
            print(" [PHASE 1] TWO-PHASE SIMPLEX: MINIMIZING ARTIFICIAL VARIABLES (y)")
            print(" Goal: Drive all artificial variables (y1, y2, ...) to 0 to find a feasible basis")
            print(" Objective: Min W = sum y_i (Max W' = - sum y_i)")
            print("=" * 82)

            # Row 0 for Phase 1: W' + sum(y_i) = 0
            # Since y_i are basic, substitute y_i = RHS_i - sum(A_ij * x_j):
            # W' - sum_j (sum_{i in Art} A_ij) * x_j = - sum_{i in Art} RHS_i
            phase1_row0 = [Fraction(0)] * (n_full + 1)
            for a_name in artificials:
                phase1_row0[full_var_names.index(a_name)] = Fraction(1)
            for i, aux in enumerate(constraint_aux):
                if aux["artificial"]:
                    for c_idx in range(n_full + 1):
                        phase1_row0[c_idx] -= matrix_A[i][c_idx]

            phase1_tableau_data = [phase1_row0] + [list(r) for r in matrix_A]
            phase1_tab = SimplexTableau(full_var_names, var_indices, phase1_tableau_data, basic_vars, obj_name="-w")

            status, solved_phase1_tab = solve_simplex(phase1_tab, phase_name="Phase 1",
                                                      forbidden_vars=set(artificials),
                                                      obj_type="min")

            if status != "OPTIMAL":
                print(f"\n [!] Phase 1 terminated with status: {status}")
                return

            w_opt = solved_phase1_tab.tableau[0][-1]
            print("\n" + "=" * 82)
            print(f" PHASE 1 COMPLETED: Optimal W' = {format_fraction(w_opt)} (Sum of artificials = {format_fraction(-w_opt)})")
            print("=" * 82)

            if w_opt != 0:
                print("\n" + "!" * 82)
                print(" [RESULT] THE LINEAR PROGRAM IS INFEASIBLE!")
                print(f" Reason: Artificial variables could not be driven to zero (Optimal sum = {format_fraction(-w_opt)} > 0).")
                print("!" * 82)
                return

            print(" [Phase 1 Success] A feasible basic solution to the original problem has been found!")

            # Clean up: If any artificial variable is still basic at zero, pivot it out
            for r_idx in range(1, solved_phase1_tab.num_rows()):
                b_var = solved_phase1_tab.basic_vars[r_idx - 1]
                if b_var in artificials:
                    print(f" • Artificial variable '{b_var}' remains basic at value 0.")
                    # Find a non-artificial candidate to pivot
                    pivoted = False
                    for c_idx, v_name in enumerate(solved_phase1_tab.var_names):
                        if v_name not in artificials and solved_phase1_tab.tableau[r_idx][c_idx] != 0:
                            print(f"   Pivoting on non-artificial variable '{v_name}' to remove '{b_var}' from basis.")
                            # Pivot
                            p_val = solved_phase1_tab.tableau[r_idx][c_idx]
                            solved_phase1_tab.tableau[r_idx] = [x / p_val for x in solved_phase1_tab.tableau[r_idx]]
                            for r_other in range(solved_phase1_tab.num_rows()):
                                if r_other != r_idx:
                                    factor = solved_phase1_tab.tableau[r_other][c_idx]
                                    solved_phase1_tab.tableau[r_other] = [
                                        solved_phase1_tab.tableau[r_other][c] - factor * solved_phase1_tab.tableau[r_idx][c]
                                        for c in range(len(solved_phase1_tab.tableau[r_other]))
                                    ]
                            solved_phase1_tab.basic_vars[r_idx - 1] = v_name
                            pivoted = True
                            break
                    if not pivoted:
                        print(f"   Row {r_idx} is redundant (all non-artificial coefficients are 0).")

            # Remove artificial columns for Phase 2
            non_art_cols = [c for c, name in enumerate(solved_phase1_tab.var_names) if name not in artificials]
            phase2_var_names = [solved_phase1_tab.var_names[c] for c in non_art_cols]
            phase2_var_indices = {name: var_indices[name] for name in phase2_var_names}

            phase2_rows = []
            for r in range(1, solved_phase1_tab.num_rows()):
                row = [solved_phase1_tab.tableau[r][c] for c in non_art_cols] + [solved_phase1_tab.tableau[r][-1]]
                phase2_rows.append(row)

            phase2_basic_vars = list(solved_phase1_tab.basic_vars)

        else:
            # No Phase 1 required! All constraints are <= with RHS >= 0
            phase2_var_names = full_var_names
            phase2_var_indices = var_indices
            phase2_rows = [list(r) for r in matrix_A]
            phase2_basic_vars = list(basic_vars)

        # -------------------------------------------------------------
        # PHASE 2 (Solve Original Problem)
        # -------------------------------------------------------------
        print("\n" + "=" * 82)
        if needs_phase1:
            print(" [PHASE 2] TWO-PHASE SIMPLEX: OPTIMIZING ORIGINAL OBJECTIVE (z)")
            print(f" Goal: {self.obj_type.upper()} original objective starting from feasible basic solution")
        else:
            print(f" SIMPLEX METHOD ({self.obj_type.upper()} Z)")
        print("=" * 82)

        # Original objective:
        # If MINIMIZE: Min Z = sum c_j x_j  <=>  Max Z' = - sum c_j x_j
        # Row 0 equation: Z' - sum (c_j') x_j = 0
        # If MAXIMIZE: Row 0 equation: Z - sum (c_j) x_j = 0
        phase2_row0 = [Fraction(0)] * (len(phase2_var_names) + 1)
        
        # Set coefficients for decision variables in Row 0
        for j, name in enumerate(self.var_names):
            c_idx = phase2_var_names.index(name)
            if self.obj_type == "max":
                # Z - sum(c_j * x_j) = 0 => coeff of x_j in Row 0 is -c_j
                phase2_row0[c_idx] = -self.obj_coeffs[j]
            else:
                # Max Z' = -Z = sum(-c_j * x_j) => Z' - sum(-c_j * x_j) = 0 => coeff of x_j is +c_j
                phase2_row0[c_idx] = self.obj_coeffs[j]

        # Canonicalize Row 0: Make coefficients of current basic variables equal to 0
        for r_idx, b_var in enumerate(phase2_basic_vars):
            c_idx = phase2_var_names.index(b_var)
            factor = phase2_row0[c_idx]
            if factor != 0:
                for col in range(len(phase2_row0)):
                    phase2_row0[col] -= factor * phase2_rows[r_idx][col]

        phase2_tableau_data = [phase2_row0] + phase2_rows
        p2_obj_name = "-z" if self.obj_type == "min" else "z"
        phase2_tab = SimplexTableau(phase2_var_names, phase2_var_indices, phase2_tableau_data,
                                    phase2_basic_vars, obj_name=p2_obj_name)

        phase_title = "Phase 2" if needs_phase1 else ""
        final_status, final_tab = solve_simplex(phase2_tab, phase_name=phase_title, obj_type=self.obj_type)

        # -------------------------------------------------------------
        # FINAL REPORT
        # -------------------------------------------------------------
        print("\n" + "=" * 82)
        print(" FINAL RESULTS SUMMARY")
        print("=" * 82)

        if final_status == "OPTIMAL":
            final_tab.print_tableau(title="Final Optimal Tableau", obj_type=self.obj_type)
            
            # Extract objective value
            raw_z = final_tab.tableau[0][-1]
            opt_z = get_display_obj_val(raw_z, self.obj_type, final_tab.obj_name)

            print(" [STATUS] OPTIMAL SOLUTION FOUND!")
            if self.obj_type == "min":
                print(f" • Tableau Row 0 Equation:  {final_tab.obj_name}* = {format_fraction(raw_z)}")
                print(f" • Optimal Minimum Value:   Min z* = -({final_tab.obj_name}*) = {format_approx(opt_z)}")
            else:
                print(f" • Optimal Objective Value ({self.obj_type.upper()} Z*): {format_approx(opt_z)}")

            # Extract variable values
            solution_dict: Dict[str, Fraction] = {name: Fraction(0) for name in self.var_names}
            aux_dict: Dict[str, Fraction] = {name: Fraction(0) for name in final_tab.var_names if name not in self.var_names}

            for r_idx in range(1, final_tab.num_rows()):
                b_var = final_tab.basic_vars[r_idx - 1]
                val = final_tab.tableau[r_idx][-1]
                if b_var in solution_dict:
                    solution_dict[b_var] = val
                elif b_var in aux_dict:
                    aux_dict[b_var] = val

            print("\n • Optimal Decision Variables:")
            for name in self.var_names:
                print(f"    {name} = {format_approx(solution_dict[name])}")

            if aux_dict:
                print("\n • Slack / Surplus Variables:")
                for name, val in aux_dict.items():
                    print(f"    {name} = {format_approx(val)}")

            # Alternative optimal solutions check
            alt_candidates = []
            for col_idx, name in enumerate(final_tab.var_names):
                if name not in final_tab.basic_vars and final_tab.tableau[0][col_idx] == 0:
                    alt_candidates.append(name)
            if alt_candidates:
                print(f"\n • Note: Non-basic variable(s) {alt_candidates} have reduced cost = 0.")
                print("   Alternative optimal solutions exist!")

        elif final_status == "UNBOUNDED":
            print(" [STATUS] THE LINEAR PROGRAM IS UNBOUNDED!")
            print(" The objective function can be increased/decreased indefinitely without violating constraints.")
        else:
            print(f" [STATUS] Terminated with status: {final_status}")
        print("=" * 82 + "\n")


# ---------------------------------------------------------------------
# MANUAL INPUT INTERFACE
# ---------------------------------------------------------------------

def prompt_guided_input() -> LinearProgram:
    """Guides the user through entering an LP problem manually."""
    print("\n" + "=" * 60)
    print(" MANUAL INPUT: STEP-BY-STEP GUIDED WIZARD")
    print("=" * 60)

    # 1. Objective Type (Always minimization of z)
    obj_type = "min"
    print("\nOptimization Goal: MINIMIZATION (Min z)")

    # 2. Number of Decision Variables
    while True:
        try:
            n_str = input("\nEnter number of decision variables (n >= 1): ").strip()
            num_vars = int(n_str)
            if num_vars < 1:
                print("Number of variables must be at least 1.")
                continue
            break
        except ValueError:
            print("Please enter a valid integer.")

    var_names = [f"x{j+1}" for j in range(num_vars)]

    # 3. Number of Constraints
    while True:
        try:
            m_str = input("Enter number of constraints (m >= 1): ").strip()
            num_constraints = int(m_str)
            if num_constraints < 1:
                print("Number of constraints must be at least 1.")
                continue
            break
        except ValueError:
            print("Please enter a valid integer.")

    # 4. Objective Coefficients
    print(f"\nObjective Function Coefficients for {var_names}:")
    print(" (Enter numbers, decimals like 0.5, or fractions like 1/2 or -5/3)")
    while True:
        coeffs_str = input(f"Enter {num_vars} coefficient(s) separated by spaces: ").strip()
        parts = coeffs_str.split()
        if len(parts) != num_vars:
            print(f"Expected {num_vars} values, got {len(parts)}. Please try again.")
            continue
        try:
            obj_coeffs = [parse_fraction(p) for p in parts]
            break
        except Exception as e:
            print(f"Error parsing coefficients: {e}. Please try again.")

    # 5. Constraints
    constraints = []
    print("\nEnter each constraint:")
    print(" You can enter in one line, e.g.: '1 2 <= 4' or '1/2 -3/2 >= 0'")
    for i in range(num_constraints):
        while True:
            line = input(f"\nConstraint #{i+1}: ").strip()
            if not line:
                print("Cannot be empty. Please enter constraint.")
                continue
            
            # Detect relation <=, >=, or =
            rel = None
            for candidate in ("<=", ">=", "="):
                if candidate in line:
                    rel = candidate
                    break
            
            if rel:
                lhs_str, rhs_str = line.split(rel, 1)
                lhs_parts = lhs_str.strip().split()
                rhs_val_str = rhs_str.strip()
                if len(lhs_parts) != num_vars:
                    print(f"Expected {num_vars} LHS coefficients, got {len(lhs_parts)}. Try again.")
                    continue
                try:
                    coeffs = [parse_fraction(p) for p in lhs_parts]
                    rhs = parse_fraction(rhs_val_str)
                    constraints.append((coeffs, rel, rhs))
                    break
                except Exception as e:
                    print(f"Error parsing constraint: {e}. Try again.")
            else:
                # Fallback: prompt piece by piece
                print(f"Could not find '<=', '>=', or '=' in '{line}'. Let's enter piece by piece:")
                while True:
                    c_str = input(f"  Enter {num_vars} LHS coefficients for Constraint #{i+1}: ").strip()
                    parts = c_str.split()
                    if len(parts) == num_vars:
                        try:
                            coeffs = [parse_fraction(p) for p in parts]
                            break
                        except Exception as e:
                            print(f"  Error: {e}")
                    else:
                        print(f"  Expected {num_vars} values.")
                while True:
                    r_str = input("  Enter relation (<=, >=, =): ").strip()
                    if r_str in ("<=", ">=", "="):
                        rel = r_str
                        break
                    print("  Please enter <=, >=, or =.")
                while True:
                    rhs_s = input("  Enter RHS value (b): ").strip()
                    try:
                        rhs = parse_fraction(rhs_s)
                        break
                    except Exception as e:
                        print(f"  Error: {e}")
                constraints.append((coeffs, rel, rhs))
                break

    return LinearProgram(obj_type, obj_coeffs, constraints, var_names)


def prompt_multiline_input() -> LinearProgram:
    """Allows pasting a full problem block."""
    print("\n" + "=" * 60)
    print(" PASTE / MULTI-LINE INPUT")
    print("=" * 60)
    print("Format:")
    print("  Line 1: max or min")
    print("  Line 2: Objective coefficients (e.g., '10 -57 -9 -24')")
    print("  Subsequent lines: one constraint per line (e.g., '0.5 -5.5 -2.5 9 <= 0')")
    print("  End input with an empty line or EOF (Ctrl+D / Ctrl+Z).")
    print("-" * 60)

    lines = []
    while True:
        try:
            line = input()
            if not line.strip() and len(lines) >= 3:
                break
            if line.strip():
                lines.append(line.strip())
        except EOFError:
            break

    if len(lines) < 3:
        raise ValueError("Insufficient lines entered. Need objective type, coefficients, and at least 1 constraint.")

    obj_type = lines[0].lower()
    if "min" in obj_type:
        obj_type = "min"
    else:
        obj_type = "max"

    obj_coeffs = [parse_fraction(p) for p in lines[1].split()]
    num_vars = len(obj_coeffs)
    var_names = [f"x{j+1}" for j in range(num_vars)]

    constraints = []
    for c_line in lines[2:]:
        rel = None
        for candidate in ("<=", ">=", "="):
            if candidate in c_line:
                rel = candidate
                break
        if not rel:
            raise ValueError(f"Constraint line missing relation (<=, >=, =): '{c_line}'")
        lhs_str, rhs_str = c_line.split(rel, 1)
        coeffs = [parse_fraction(p) for p in lhs_str.strip().split()]
        if len(coeffs) != num_vars:
            raise ValueError(f"Expected {num_vars} LHS coefficients in '{c_line}', got {len(coeffs)}")
        rhs = parse_fraction(rhs_str.strip())
        constraints.append((coeffs, rel, rhs))

    return LinearProgram(obj_type, obj_coeffs, constraints, var_names)


# ---------------------------------------------------------------------
# CLASSIC BENCHMARK EXAMPLES
# ---------------------------------------------------------------------

def get_beale_cycling_example() -> LinearProgram:
    """
    Beale's Classic Cycling Example (1955):
    Under Dantzig's largest-coefficient rule, this problem CYCLES forever (period 6).
    Under Bland's Smallest-Subscript Rule, cycling is prevented, and the optimal solution is reached!
    
    Minimize z = -10*x1 + 57*x2 + 9*x3 + 24*x4
    s.t.
      0.5*x1 - 5.5*x2 - 2.5*x3 + 9*x4 <= 0
      0.5*x1 - 1.5*x2 - 0.5*x3 + 1*x4 <= 0
      1*x1                            <= 1
      x1, x2, x3, x4 >= 0
    Optimal solution: x1 = 1, x2 = 0, x3 = 1, x4 = 0, z* = -1.
    """
    obj_type = "min"
    obj_coeffs = [Fraction(-10), Fraction(57), Fraction(9), Fraction(24)]
    constraints = [
        ([Fraction(1, 2), Fraction(-11, 2), Fraction(-5, 2), Fraction(9)], "<=", Fraction(0)),
        ([Fraction(1, 2), Fraction(-3, 2), Fraction(-1, 2), Fraction(1)], "<=", Fraction(0)),
        ([Fraction(1), Fraction(0), Fraction(0), Fraction(0)], "<=", Fraction(1)),
    ]
    var_names = ["x1", "x2", "x3", "x4"]
    return LinearProgram(obj_type, obj_coeffs, constraints, var_names)


def get_twophase_example() -> LinearProgram:
    """
    Example requiring Two-Phase Simplex:
    Minimize z = 4*x1 + x2
    s.t.
      3*x1 + x2  = 3
      4*x1 + 3*x2 >= 6
      x1 + 2*x2 <= 4
      x1, x2 >= 0
    Optimal solution: x1 = 2/5, x2 = 9/5, z* = 17/5 = 3.4.
    """
    obj_type = "min"
    obj_coeffs = [Fraction(4), Fraction(1)]
    constraints = [
        ([Fraction(3), Fraction(1)], "=", Fraction(3)),
        ([Fraction(4), Fraction(3)], ">=", Fraction(6)),
        ([Fraction(1), Fraction(2)], "<=", Fraction(4)),
    ]
    var_names = ["x1", "x2"]
    return LinearProgram(obj_type, obj_coeffs, constraints, var_names, force_two_phase=True)


def get_standard_2d_example() -> LinearProgram:
    """
    Standard 2D textbook minimization problem:
    Minimize z = -3*x1 - 5*x2
    s.t.
      x1        <= 4
      2*x2      <= 12
      3*x1 + 2*x2 <= 18
      x1, x2 >= 0
    Optimal: x1 = 2, x2 = 6, z* = -36.
    """
    obj_type = "min"
    obj_coeffs = [Fraction(-3), Fraction(-5)]
    constraints = [
        ([Fraction(1), Fraction(0)], "<=", Fraction(4)),
        ([Fraction(0), Fraction(2)], "<=", Fraction(12)),
        ([Fraction(3), Fraction(2)], "<=", Fraction(18)),
    ]
    var_names = ["x1", "x2"]
    return LinearProgram(obj_type, obj_coeffs, constraints, var_names)


# ---------------------------------------------------------------------
# CSV TABLEAU IMPORT
# ---------------------------------------------------------------------

def load_tableau_from_csv(file_path: str) -> Tuple[SimplexTableau, str]:
    """
    Imports an initial Simplex Tableau from a CSV file.
    
    Supported formats:
    - Comma, semicolon, or tab delimiters.
    - Standard fractions (1/2, -11/2), decimals (0.5 or European 0,5), and integers.
    - RHS column on the left (Basis, RHS, x1, x2, ...) or on the right.
    - With or without an explicit 'z' column (z column is stripped automatically).
    - Automatic basic variable identification via unit columns if not specified in column 0.
    """
    resolved_path = resolve_csv_path(file_path)
    if os.path.exists(resolved_path):
        file_path = resolved_path
    elif not os.path.exists(file_path):
        raise FileNotFoundError(f"File '{file_path}' not found (checked '{file_path}' and '{resolved_path}').")

    with open(file_path, "r", encoding="utf-8-sig") as f:
        lines = [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]

    if not lines:
        raise ValueError(f"CSV file '{file_path}' is empty or contains only comments.")

    # Detect delimiter
    first_line = lines[0]
    count_semi = first_line.count(';')
    count_comma = first_line.count(',')
    count_tab = first_line.count('\t')
    delimiter = ','
    if count_semi > count_comma and count_semi > count_tab:
        delimiter = ';'
    elif count_tab > count_comma and count_tab > count_semi:
        delimiter = '\t'

    reader = csv.reader(lines, delimiter=delimiter)
    raw_rows = [[cell.strip() for cell in row] for row in reader]
    raw_rows = [row for row in raw_rows if any(cell for cell in row)]

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

    if len(data_rows) < 2:
        raise ValueError("Tableau must contain at least an objective row (Row 0) and at least 1 constraint row.")

    # Check if first column is the Basis column
    has_basis_col = False
    if header_row and header_row[0].lower() in ('basis', 'base', 'basic', 'b', 'var', 'variabile', ''):
        has_basis_col = True
    elif not is_number(data_rows[0][0]) and not is_number(data_rows[1][0]):
        has_basis_col = True

    basic_vars: List[str] = []
    obj_name = "-z"
    detected_obj_type = "min"

    if has_basis_col:
        row0_raw = data_rows[0][0].strip()
        row0_lower = row0_raw.lower()
        if row0_lower in ("w", "-w"):
            obj_name = "-w"
            detected_obj_type = "min"
        else:
            obj_name = "-z"
            detected_obj_type = "min"
        for r in data_rows[1:]:
            b_name = r[0].strip()
            if not b_name:
                raise ValueError(f"Missing basic variable name in row: {r}")
            basic_vars.append(b_name)
        col_headers = header_row[1:] if header_row else None
        cleaned_data = [r[1:] for r in data_rows]
    else:
        col_headers = header_row
        cleaned_data = data_rows

    # Check for an explicit z/w column (e.g. [1, 0, 0, ...])
    z_col_idx = None
    if col_headers:
        for idx, h in enumerate(col_headers):
            if h.lower() in ('z', 'w', '-z', 'obj'):
                z_col_idx = idx
                break
    else:
        if all(parse_fraction(cleaned_data[r][0]) == (Fraction(1) if r == 0 else Fraction(0)) for r in range(len(cleaned_data))):
            z_col_idx = 0

    if z_col_idx is not None:
        if col_headers:
            col_headers.pop(z_col_idx)
        cleaned_data = [r[:z_col_idx] + r[z_col_idx+1:] for r in cleaned_data]

    # Check if RHS is at index 0 (on the left: Basis, RHS, x1, x2, ...)
    rhs_on_left = False
    if col_headers and col_headers[0].lower() in ('rhs', 'b', 'val', 'costo', 'term'):
        rhs_on_left = True

    if rhs_on_left:
        var_names = col_headers[1:]
        # Rearrange matrix rows so RHS is stored at index -1 internally
        matrix: List[List[Fraction]] = []
        for row_idx, r in enumerate(cleaned_data):
            rhs_val = parse_fraction(r[0])
            var_vals = [parse_fraction(cell) for cell in r[1:]]
            matrix.append(var_vals + [rhs_val])
    else:
        num_cols = len(cleaned_data[0])
        if col_headers:
            var_names = col_headers[:-1]
        else:
            var_names = [f"x{j+1}" for j in range(num_cols - 1)]

        # Parse numeric matrix directly (RHS is last column)
        matrix = []
        for row_idx, r in enumerate(cleaned_data):
            if len(r) != num_cols:
                raise ValueError(f"Inconsistent column count at row {row_idx}: expected {num_cols}, got {len(r)}")
            matrix.append([parse_fraction(cell) for cell in r])

    # Convert any slack variables starting with 's' to 'x' variables
    s_vars = [name for name in var_names if name.lower().startswith('s') and name[1:].isdigit()]
    if s_vars:
        existing_x = [int(v[1:]) for v in var_names if v.startswith('x') and v[1:].isdigit()]
        curr_x = max(existing_x, default=0) + 1
        rename_map = {}
        for s in s_vars:
            rename_map[s] = f"x{curr_x}"
            curr_x += 1
        var_names = [rename_map.get(name, name) for name in var_names]
        basic_vars = [rename_map.get(name, name) for name in basic_vars]

    # If no basic variables were provided, infer them from unit columns
    if not basic_vars:
        for r in range(1, len(matrix)):
            found = False
            for c in range(len(var_names)):
                # Unit column check
                if matrix[r][c] == 1 and all(matrix[other][c] == 0 for other in range(len(matrix)) if other != r):
                    basic_vars.append(var_names[c])
                    found = True
                    break
            if not found:
                basic_vars.append(f"x{len(var_names) + r}")

    var_indices = {name: idx + 1 for idx, name in enumerate(var_names)}
    tab = SimplexTableau(var_names, var_indices, matrix, basic_vars, obj_name=obj_name)
    return tab, detected_obj_type


def solve_imported_tableau(tab: SimplexTableau, obj_type: str = "min", force_two_phase: bool = True):
    """Solves an imported tableau using Bland's rule and prints the final report."""
    print("\n" + "=" * 82)
    mode_str = "FORCED TWO-PHASE" if force_two_phase else "STANDARD SIMPLEX"
    print(f" SIMPLEX METHOD ON IMPORTED CSV TABLEAU ({mode_str})")
    print("=" * 82)
    is_min = tab.obj_name.startswith("-") or obj_type == "min"
    direction_str = "MINIMIZATION" if is_min else "MAXIMIZATION"
    print(f" Target Direction: {direction_str}")
    print(" Variables and Bland's Rule Index Ordering:")
    for v_name in tab.var_names:
        print(f"   Index {tab.var_indices[v_name]:<2}: {v_name}")
    print(f" Initial Basis: {tab.basic_vars}")
    print("=" * 82)

    if force_two_phase:
        result = record_tableau_steps_twophase(tab, obj_type=obj_type)
        print(f" Status: {result['status']}")
        if result.get('final_solution'):
            fs = result['final_solution']
            print(f" • Optimal Objective Value: {fs['obj_val']}")
            print(f" • Final Basic Variables: {fs['basic_vars']}")
            print(f" • Solution Vector:")
            for k, v in fs['solution_vector'].items():
                print(f"    {k:<5} = {v}")
        return

    status, final_tab = solve_simplex(tab, obj_type=obj_type)

    print("\n" + "=" * 82)
    print(" FINAL RESULTS SUMMARY")
    print("=" * 82)

    if status == "OPTIMAL":
        final_tab.print_tableau(title="Final Optimal Tableau", obj_type=obj_type)
        raw_z = final_tab.tableau[0][-1]
        opt_z = get_display_obj_val(raw_z, obj_type, final_tab.obj_name)

        print(" [STATUS] OPTIMAL SOLUTION FOUND!")
        if is_min:
            clean_name = final_tab.obj_name.lstrip("-")
            print(f" • Tableau Row 0 Equation:  {final_tab.obj_name}* = {format_fraction(raw_z)}")
            print(f" • Optimal Minimum Value:   Min {clean_name}* = -({final_tab.obj_name}*) = {format_approx(opt_z)}")
        else:
            print(f" • Optimal Objective Value ({obj_type.upper()} {final_tab.obj_name}*): {format_approx(opt_z)}")

        sol_dict = {name: Fraction(0) for name in final_tab.var_names}
        for r_idx in range(1, final_tab.num_rows()):
            b_var = final_tab.basic_vars[r_idx - 1]
            sol_dict[b_var] = final_tab.tableau[r_idx][-1]

        print("\n • Solution Vector:")
        for name in final_tab.var_names:
            is_basic = name in final_tab.basic_vars
            b_tag = "(Basic)" if is_basic else "(Non-Basic = 0)"
            print(f"    {name:<5} = {format_approx(sol_dict[name]):<15} {b_tag}")

        alt_candidates = [
            name for c_idx, name in enumerate(final_tab.var_names)
            if name not in final_tab.basic_vars and final_tab.tableau[0][c_idx] == 0
        ]
        if alt_candidates:
            print(f"\n • Note: Non-basic variable(s) {alt_candidates} have reduced cost = 0.")
            print("   Alternative optimal solutions exist!")
    elif status == "UNBOUNDED":
        print(" [STATUS] THE LINEAR PROGRAM IS UNBOUNDED!")
        print(" Objective function can improve indefinitely without violating constraints.")
    else:
        print(f" [STATUS] Terminated with status: {status}")
    print("=" * 82 + "\n")


# Alias for compatibility
get_beale_problem = get_beale_cycling_example

# =====================================================================
# WEB APPLICATION INTERFACE (HTML / CSS / JAVASCRIPT)
# =====================================================================

HTML_PAGE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Simplex Method with Bland's Rule</title>
  <style>
    :root {
      --bg: #f8fafc;
      --card-bg: #ffffff;
      --border: #e2e8f0;
      --text: #0f172a;
      --text-muted: #64748b;
      --primary: #2563eb;
      --primary-hover: #1d4ed8;
      --primary-light: #eff6ff;
      --success: #16a34a;
      --success-light: #dcfce7;
      --warning: #d97706;
      --warning-light: #fef3c7;
      --purple: #7c3aed;
      --purple-light: #f5f3ff;
      --pivot-cell: #fde047;
      --pivot-col: #dbeafe;
      --pivot-row: #f3e8ff;
      --danger: #dc2626;
      --danger-light: #fee2e2;
    }
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      background: var(--bg);
      color: var(--text);
      line-height: 1.5;
      padding-bottom: 60px;
    }
    .header {
      background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%);
      color: white;
      padding: 0 24px;
      min-height: 52px;
      display: flex;
      align-items: center;
      box-shadow: 0 2px 4px rgba(0,0,0,0.15);
      position: sticky;
      top: 0;
      z-index: 100;
    }
    .header-inner {
      width: 100%;
      max-width: 1360px;
      margin: 0 auto;
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 16px;
      flex-wrap: wrap;
      padding: 6px 0;
    }
    .header-brand {
      display: flex;
      align-items: center;
      gap: 10px;
      white-space: nowrap;
    }
    .header-brand h1 {
      font-size: 16px;
      font-weight: 700;
      color: #ffffff;
      margin: 0;
      letter-spacing: -0.2px;
    }
    .header-subtitle {
      font-weight: 500;
      font-size: 12px;
      color: #94a3b8;
      margin-left: 4px;
    }
    .header-nav {
      display: flex;
      align-items: center;
      gap: 6px;
      overflow-x: auto;
    }
    .tab-btn {
      padding: 6px 12px;
      font-size: 12.5px;
      font-weight: 600;
      color: #cbd5e1;
      background: rgba(255, 255, 255, 0.06);
      border: 1px solid rgba(255, 255, 255, 0.12);
      border-radius: 6px;
      cursor: pointer;
      display: flex;
      align-items: center;
      gap: 6px;
      white-space: nowrap;
      transition: all 0.15s ease;
    }
    .tab-btn:hover {
      color: #ffffff;
      background: rgba(255, 255, 255, 0.14);
      border-color: rgba(255, 255, 255, 0.25);
    }
    .tab-btn.active {
      color: #ffffff;
      background: #2563eb;
      border-color: #3b82f6;
      box-shadow: 0 1px 3px rgba(0, 0, 0, 0.25);
    }
    .tab-content { display: none; }
    .tab-content.active { display: block; }
    .container { max-width: 1360px; margin: 0 auto; padding: 20px 24px; }

    /* Welcome / Landing Choice Grid */
    .welcome-hero {
      background: linear-gradient(135deg, #f8fafc 0%, #e2e8f0 100%);
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 24px 28px;
      margin-bottom: 24px;
      text-align: center;
    }
    .welcome-hero h2 {
      font-size: 22px;
      font-weight: 800;
      color: #0f172a;
      margin-bottom: 6px;
    }
    .welcome-hero p {
      font-size: 14px;
      color: #64748b;
    }
    .choice-cards-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
      gap: 20px;
      margin-bottom: 24px;
    }
    .choice-card {
      background: var(--card-bg);
      border: 1.5px solid var(--border);
      border-radius: 14px;
      padding: 24px;
      display: flex;
      flex-direction: column;
      box-shadow: 0 2px 6px rgba(0,0,0,0.04);
      transition: transform 0.15s, box-shadow 0.15s, border-color 0.15s;
    }
    .choice-card:hover {
      transform: translateY(-2px);
      box-shadow: 0 8px 16px rgba(0,0,0,0.06);
      border-color: #93c5fd;
    }
    .choice-header {
      display: flex;
      align-items: center;
      gap: 14px;
      margin-bottom: 12px;
    }
    .choice-icon {
      width: 38px;
      height: 38px;
      border-radius: 8px;
      background: #eff6ff;
      border: 1.5px solid #bfdbfe;
      display: flex;
      align-items: center;
      justify-content: center;
      font-size: 16px;
      font-weight: 800;
      color: #2563eb;
      flex-shrink: 0;
    }
    .choice-card h3 {
      font-size: 16px;
      font-weight: 700;
      color: #0f172a;
    }
    .choice-subtitle {
      font-size: 12px;
      color: #64748b;
    }
    .choice-desc {
      font-size: 13px;
      color: #475569;
      line-height: 1.5;
      flex-grow: 1;
    }
    .benchmark-btn {
      background: #f8fafc;
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 10px 12px;
      text-align: left;
      cursor: pointer;
      transition: all 0.15s;
      width: 100%;
    }
    .benchmark-btn:hover {
      background: #eff6ff;
      border-color: #3b82f6;
    }
    .bench-title {
      font-size: 13px;
      font-weight: 700;
      color: #1e293b;
    }
    .bench-desc {
      font-size: 11px;
      color: #64748b;
      margin-top: 2px;
    }

    /* Cards */
    .card {
      background: var(--card-bg);
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 20px;
      margin-bottom: 20px;
      box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    }
    .card-title {
      font-size: 16px;
      font-weight: 700;
      color: var(--text);
      margin-bottom: 14px;
      display: flex;
      align-items: center;
      justify-content: space-between;
    }

    /* Consolidated Navigation Box (Iteration Row + Sub-step Row) */
    .stepper-box {
      background: #f1f5f9;
      border-radius: 10px;
      padding: 12px 16px;
      margin-bottom: 14px;
      display: flex;
      flex-direction: column;
      gap: 10px;
      border: 1px solid var(--border);
    }
    .stepper-row {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      width: 100%;
    }
    .iteration-side-left {
      display: flex;
      align-items: center;
      gap: 8px;
      min-width: 130px;
    }
    .iteration-side-right {
      display: flex;
      align-items: center;
      gap: 8px;
      min-width: 130px;
      justify-content: flex-end;
    }
    .iteration-controls {
      display: flex;
      align-items: center;
      justify-content: center;
      gap: 8px;
      flex: 1;
    }
    .iter-counter {
      font-weight: 800;
      font-size: 14px;
      color: #0f172a;
      padding: 6px 16px;
      background: #ffffff;
      border: 1px solid var(--border);
      border-radius: 6px;
      min-width: 130px;
      text-align: center;
      display: inline-block;
    }

    /* Sub-step Row (Directly below iteration row in the same box) */
    .stepper-row-substep {
      display: flex;
      align-items: center;
      justify-content: center;
      padding-top: 8px;
      border-top: 1px solid #e2e8f0;
    }
    .substep-controls {
      display: flex;
      align-items: center;
      justify-content: center;
      gap: 10px;
    }
    .substep-counter-badge {
      font-weight: 700;
      font-size: 13px;
      color: #334155;
      padding: 5px 14px;
      background: #ffffff;
      border: 1px solid #cbd5e1;
      border-radius: 6px;
      min-width: 180px;
      text-align: center;
      display: inline-block;
    }
    .btn-group { display: flex; align-items: center; gap: 6px; }
    .btn {
      background: var(--primary);
      color: white;
      border: none;
      padding: 8px 14px;
      border-radius: 6px;
      font-size: 13px;
      font-weight: 600;
      cursor: pointer;
      display: inline-flex;
      align-items: center;
      justify-content: center;
      gap: 6px;
      transition: background 0.15s;
    }
    .btn:hover { background: var(--primary-hover); }
    .btn-secondary {
      background: white;
      color: var(--text);
      border: 1px solid var(--border);
    }
    .btn-secondary:hover { background: #f8fafc; }
    .btn:disabled { opacity: 0.5; cursor: not-allowed; }
    .badge {
      padding: 4px 10px;
      border-radius: 20px;
      font-size: 12px;
      font-weight: 700;
    }
    .badge-blue { background: var(--primary-light); color: var(--primary); }
    .badge-green { background: var(--success-light); color: var(--success); }
    .badge-amber { background: var(--warning-light); color: #92400e; border: 1px solid #f59e0b; }
    .badge-purple { background: var(--purple-light); color: #5b21b6; border: 1px solid #8b5cf6; }
    .badge-red { background: var(--danger-light); color: var(--danger); border: 1px solid #f87171; }

    /* Phase Banners */
    .phase-banner {
      border-radius: 8px;
      padding: 12px 18px;
      margin-bottom: 14px;
      display: flex;
      align-items: center;
      gap: 14px;
      font-size: 13px;
    }
    .phase-banner-p1 {
      background: #fffbeb;
      border: 1.5px solid #f59e0b;
      color: #92400e;
    }
    .phase-banner-p2 {
      background: #faf5ff;
      border: 1.5px solid #8b5cf6;
      color: #5b21b6;
    }
    .phase-banner-std {
      background: #eff6ff;
      border: 1.5px solid #3b82f6;
      color: #1e40af;
    }
    .phase-badge-p1 {
      background: #d97706;
      color: white;
      font-weight: 800;
      font-size: 11px;
      letter-spacing: 0.5px;
      padding: 4px 10px;
      border-radius: 4px;
      text-transform: uppercase;
      white-space: nowrap;
    }
    .phase-badge-p2 {
      background: #7c3aed;
      color: white;
      font-weight: 800;
      font-size: 11px;
      letter-spacing: 0.5px;
      padding: 4px 10px;
      border-radius: 4px;
      text-transform: uppercase;
      white-space: nowrap;
    }
    .phase-badge-std {
      background: #2563eb;
      color: white;
      font-weight: 800;
      font-size: 11px;
      letter-spacing: 0.5px;
      padding: 4px 10px;
      border-radius: 4px;
      text-transform: uppercase;
      white-space: nowrap;
    }

    /* Tableau Table */
    .table-container {
      overflow-x: auto;
      margin-bottom: 16px;
      border: 2px solid #64748b;
      border-radius: 8px;
    }
    table.tableau {
      width: 100%;
      border-collapse: collapse;
      font-size: 13px;
      text-align: right;
    }
    table.tableau th, table.tableau td {
      padding: 10px 14px;
      border: 2px solid #94a3b8;
      white-space: nowrap;
    }
    table.tableau th {
      background: #f8fafc;
      font-weight: 700;
      color: var(--text);
      text-align: center;
      border-bottom: 2.5px solid #64748b;
    }
    table.tableau td.basis-cell {
      text-align: center;
      font-weight: 700;
      background: #f8fafc;
    }

    /* Tableau Color Legend */
    .tableau-legend {
      display: flex;
      flex-wrap: wrap;
      align-items: center;
      gap: 14px;
      margin-bottom: 14px;
      padding: 10px 14px;
      background: #f8fafc;
      border: 1px solid var(--border);
      border-radius: 8px;
      font-size: 12px;
    }
    .legend-item {
      display: inline-flex;
      align-items: center;
      gap: 6px;
      color: #334155;
    }
    .legend-swatch {
      width: 14px;
      height: 14px;
      border-radius: 3px;
      display: inline-block;
      flex-shrink: 0;
    }
    .swatch-basic { background: #fee2e2; border: 1.5px solid #ef4444; }
    .swatch-nonbasic { background: #eff6ff; border: 1.5px solid #3b82f6; }
    .swatch-pivot { background: #fef08a; border: 1.5px solid #ca8a04; }
    .swatch-active-row { background: #eff6ff; border: 1.5px solid #3b82f6; }
    .swatch-elim { background: #dcfce7; border: 1.5px solid #16a34a; }

    /* Header Column Boxes & Tags */
    .col-header-box {
      display: flex;
      flex-direction: column;
      align-items: center;
      gap: 3px;
    }
    .col-name-basic {
      color: #dc2626;
      font-weight: 800;
      font-size: 13px;
    }
    .col-name-nonbasic {
      color: #2563eb;
      font-weight: 800;
      font-size: 13px;
    }
    .col-tag {
      font-size: 9px;
      font-weight: 700;
      padding: 1px 6px;
      border-radius: 4px;
      text-transform: uppercase;
      letter-spacing: 0.3px;
      line-height: 1.2;
    }
    .tag-basic {
      background: #fee2e2;
      color: #991b1b;
      border: 1px solid #f87171;
    }
    .tag-nonbasic {
      background: #eff6ff;
      color: #1e40af;
      border: 1px solid #bfdbfe;
    }
    .tag-rhs {
      background: #ddd6fe;
      color: #5b21b6;
      border: 1px solid #a78bfa;
    }
    .tag-y {
      background: #fef08a;
      color: #854d0e;
      border: 1px solid #eab308;
    }
    .tag-pivot {
      background: #fef3c7;
      color: #92400e;
      border: 1px solid #f59e0b;
    }

    /* Row Index Column */
    table.tableau th.col-row-idx,
    table.tableau td.cell-row-idx {
      width: 65px;
      min-width: 65px;
      text-align: center;
      background: #f8fafc;
      border-right: 2px solid #94a3b8;
    }
    .row-badge {
      display: inline-block;
      font-size: 11px;
      font-weight: 700;
      color: #64748b;
      background: #f1f5f9;
      border: 1px solid #cbd5e1;
      border-radius: 4px;
      padding: 2px 7px;
      letter-spacing: 0.2px;
    }
    .row-badge-pivot {
      color: #92400e;
      background: #fef3c7;
      border: 1.5px solid #f59e0b;
      font-weight: 800;
    }

    /* Basis Column */
    table.tableau th.col-basis,
    table.tableau td.cell-basis {
      width: 75px;
      min-width: 75px;
      text-align: center;
      background: #f8fafc;
      border-right: 2px solid #94a3b8;
    }

    /* RHS Column (Clean border divider, neutral cell) */
    table.tableau th.col-rhs,
    table.tableau td.cell-rhs {
      background: #ffffff;
      color: var(--text);
      font-weight: 700;
      border-left: 2px solid #94a3b8;
      text-align: right;
    }
    table.tableau th.col-rhs {
      background: #f8fafc;
    }

    /* Active Pivot Column Header Highlight */
    table.tableau th.col-pivot-active {
      background: #fffbeb;
      box-shadow: inset 0 0 0 2px #ca8a04;
    }

    /* Prominent dividers: RHS Column right border & -z row bottom border */
    table.tableau th.col-rhs,
    table.tableau td.cell-rhs,
    table.tableau th:nth-child(3),
    table.tableau td:nth-child(3),
    .editor-table th.col-rhs,
    .editor-table td.cell-rhs,
    .editor-table th:nth-child(3),
    .editor-table td:nth-child(3) {
      border-right: 3.5px solid #0f172a !important;
    }

    table.tableau tbody tr.cost-row td,
    table.tableau tbody tr.cost-row th,
    table.tableau tbody tr:first-child td,
    table.tableau tbody tr:first-child th,
    .editor-table tbody tr.cost-row td,
    .editor-table tbody tr.cost-row th,
    .editor-table tbody tr:first-child td,
    .editor-table tbody tr:first-child th {
      border-bottom: 3.5px solid #0f172a !important;
    }

    /* Sub-step Elementary Row Operation Focus Styles */
    table.tableau td.cell-pivot {
      background: #fef08a !important;
      color: #78350f !important;
      font-weight: 900 !important;
      box-shadow: inset 0 0 0 2.5px #ca8a04 !important;
    }
    table.tableau td.cell-dimmed {
      background: #fafafa !important;
      color: #94a3b8 !important;
      opacity: 0.65;
    }
    table.tableau th.col-dimmed {
      opacity: 0.55;
    }
    table.tableau td.cell-active-norm {
      background: #fef9c3 !important;
      color: #854d0e !important;
      font-weight: 600;
    }
    table.tableau td.cell-target-row {
      background: #eff6ff !important;
      color: #1e40af !important;
      font-weight: 600;
    }
    table.tableau td.cell-pivot-row {
      background: #fffbeb !important;
      color: #92400e !important;
      font-weight: 600;
    }
    table.tableau td.cell-eliminated {
      background: #dcfce7 !important;
      color: #15803d !important;
      font-weight: 900 !important;
      box-shadow: inset 0 0 0 2.5px #22c55e !important;
    }

    /* Sub-step Operation Banner */
    .substep-banner {
      background: #ffffff;
      border: 1.5px solid var(--border);
      border-radius: 10px;
      padding: 14px 18px;
      margin-bottom: 16px;
      display: flex;
      flex-direction: column;
      gap: 8px;
      box-shadow: 0 1px 4px rgba(0,0,0,0.04);
    }
    .substep-banner-header {
      display: flex;
      align-items: center;
      justify-content: space-between;
      flex-wrap: wrap;
      gap: 10px;
    }
    .substep-pill {
      font-size: 11px;
      font-weight: 800;
      letter-spacing: 0.5px;
      padding: 4px 12px;
      border-radius: 14px;
      text-transform: uppercase;
    }
    .substep-pill-pivot { background: #fef3c7; color: #92400e; border: 1px solid #f59e0b; }
    .substep-pill-norm { background: #fef9c3; color: #854d0e; border: 1px solid #eab308; }
    .substep-pill-elim { background: #eff6ff; color: #1e40af; border: 1px solid #3b82f6; }
    .substep-pill-done { background: #dcfce7; color: #166534; border: 1px solid #22c55e; }
    .substep-pill-opt { background: #dcfce7; color: #14532d; border: 1.5px solid #16a34a; }
    .substep-pill-unb { background: #fee2e2; color: #991b1b; border: 1px solid #ef4444; }

    .substep-formula {
      font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
      font-size: 14px;
      font-weight: 700;
      color: #0f172a;
      background: #f8fafc;
      padding: 7px 14px;
      border-radius: 6px;
      border-left: 4px solid var(--primary);
      display: inline-block;
    }
    .substep-desc {
      font-size: 13px;
      color: #475569;
      line-height: 1.45;
    }

    /* Analysis Box */
    .analysis-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
      gap: 16px;
      margin-top: 16px;
    }
    .info-card {
      background: #f8fafc;
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 14px;
    }
    .info-card h4 {
      font-size: 13px;
      font-weight: 700;
      margin-bottom: 8px;
      color: var(--text);
      display: flex;
      align-items: center;
      gap: 6px;
    }
    .info-card ul {
      list-style-type: none;
      font-size: 13px;
    }
    .info-card li { margin-bottom: 4px; }

    /* Final Result Banner */
    .result-banner {
      background: var(--success-light);
      border: 1px solid #86efac;
      border-radius: 10px;
      padding: 16px 20px;
      margin-top: 16px;
      display: flex;
      flex-direction: column;
      gap: 8px;
    }
    .result-banner h3 {
      color: #14532d;
      font-size: 16px;
      display: flex;
      align-items: center;
      gap: 8px;
    }

    /* CSV Editor Table */
    .editor-table {
      width: 100%;
      border-collapse: collapse;
      margin-top: 12px;
    }
    .editor-table th, .editor-table td {
      border: 2px solid #94a3b8;
      padding: 6px 8px;
    }
    .editor-table th {
      background: #f8fafc;
      font-weight: 700;
      border-bottom: 2.5px solid #64748b;
    }
    .editor-table input {
      width: 100%;
      border: none;
      padding: 6px;
      font-size: 13px;
      font-family: inherit;
      background: transparent;
      outline: none;
    }
    .editor-table input:focus {
      background: #eff6ff;
      border-radius: 4px;
    }
    .select-input {
      padding: 8px 12px;
      border-radius: 6px;
      border: 1px solid var(--border);
      background: white;
      font-size: 14px;
      font-family: inherit;
      width: 100%;
    }
    .form-group {
      display: flex;
      flex-direction: column;
      gap: 6px;
      margin-bottom: 12px;
    }
    .form-row {
      display: flex;
      gap: 12px;
      align-items: center;
      flex-wrap: wrap;
      margin-bottom: 12px;
    }

    /* Modern Toggle Switch & Card */
    .toggle-card {
      display: inline-flex;
      align-items: center;
      gap: 10px;
      padding: 6px 12px;
      background: #ffffff;
      border: 1.5px solid #cbd5e1;
      border-radius: 8px;
      cursor: pointer;
      user-select: none;
      transition: all 0.15s ease;
      min-height: 38px;
      box-sizing: border-box;
    }
    .toggle-card:hover {
      background: #f8fafc;
      border-color: #94a3b8;
    }
    .toggle-card.active,
    .toggle-card:has(input:checked) {
      background: #eff6ff;
      border-color: #3b82f6;
    }
    .toggle-switch {
      position: relative;
      display: inline-block;
      width: 36px;
      height: 20px;
      flex-shrink: 0;
    }
    .toggle-switch input {
      opacity: 0;
      width: 0;
      height: 0;
      position: absolute;
    }
    .toggle-slider {
      position: absolute;
      cursor: pointer;
      top: 0;
      left: 0;
      right: 0;
      bottom: 0;
      background-color: #cbd5e1;
      transition: transform 0.2s ease, background-color 0.2s ease;
      border-radius: 20px;
    }
    .toggle-slider:before {
      position: absolute;
      content: "";
      height: 14px;
      width: 14px;
      left: 3px;
      bottom: 3px;
      background-color: white;
      transition: transform 0.2s ease;
      border-radius: 50%;
      box-shadow: 0 1px 3px rgba(0, 0, 0, 0.25);
    }
    .toggle-switch input:checked + .toggle-slider {
      background-color: #2563eb;
    }
    .toggle-switch input:checked + .toggle-slider:before {
      transform: translateX(16px);
    }
    .toggle-content {
      display: flex;
      flex-direction: column;
      justify-content: center;
      line-height: 1.25;
    }
    .toggle-title {
      font-size: 13px;
      font-weight: 600;
      color: #1e293b;
    }
    .toggle-card.active .toggle-title,
    .toggle-card:has(input:checked) .toggle-title {
      color: #1d4ed8;
    }
    .toggle-desc {
      font-size: 11px;
      color: #64748b;
      font-weight: 400;
      margin-top: 1px;
    }
  </style>
</head>
<body>

  <header class="header">
    <div class="header-inner">
      <div class="header-brand">
        <h1>Simplex Method <span class="header-subtitle">&bull; Bland's Rule</span></h1>
      </div>
      <div class="header-nav">
        <button id="tab-btn-choose" class="tab-btn active" onclick="switchTab('choose')">Choose Problem</button>
        <button id="tab-btn-visualizer" class="tab-btn" onclick="switchTab('visualizer')">Step-by-Step Tableau</button>
        <button id="tab-btn-csv" class="tab-btn" onclick="switchTab('csv')">CSV Editor</button>
        <button id="tab-btn-manual" class="tab-btn" onclick="switchTab('manual')">Manual LP Input</button>
        <button id="tab-btn-demos" class="tab-btn" onclick="switchTab('demos')">Benchmarks</button>
      </div>
    </div>
  </header>

  <div class="container">

    <!-- TAB 1: CHOOSE PROBLEM (HOME LANDING SCREEN) -->
    <div id="tab-choose" class="tab-content active">
      <div class="welcome-hero">
        <h2>Simplex Method Solver (Bland's Rule)</h2>
        <p>Choose an option below to load a problem and visualize every simplex iteration step-by-step:</p>
      </div>

      <div class="choice-cards-grid">
        
        <!-- Option 1: CSV File Input -->
        <div class="choice-card">
          <div class="choice-header">
            <div class="choice-icon">1</div>
            <div>
              <h3>1. Choose CSV & Run Simplex</h3>
              <p class="choice-subtitle">Select or upload a starting tableau</p>
            </div>
          </div>
          <p class="choice-desc">
            Load an initial tableau formatted as <code>Basis, RHS, x1, x2, ...</code>.
            Select a CSV file in your workspace or upload a file. You can solve directly or inspect and edit the table.
          </p>

          <div class="form-group" style="margin-top: 14px;">
            <label style="font-weight: 700; font-size: 13px;">Select CSV file in workspace:</label>
            <select id="landing-csv-select" class="select-input" onchange="syncLandingCsv(this.value)">
              <option value="tableau.csv" selected>tableau.csv</option>
              <option value="example_standard.csv">example_standard.csv</option>
            </select>
          </div>

          <div style="margin-top: 10px;">
            <label style="font-weight: 600; font-size: 12px; color: var(--text-muted); cursor: pointer;">
              Or upload custom CSV:
              <input type="file" id="landing-file-upload" accept=".csv" class="select-input" style="margin-top: 4px;" onchange="handleLandingFileUpload(event)">
            </label>
          </div>

          <div style="display: flex; gap: 8px; margin-top: 18px;">
            <button type="button" class="btn" style="flex: 1.3;" onclick="solveLandingCsv()">
              Run Simplex on CSV
            </button>
            <button type="button" class="btn btn-secondary" style="flex: 1;" onclick="previewLandingCsv()">
              View / Edit
            </button>
          </div>
        </div>

        <!-- Option 2: Predefined Benchmarks -->
        <div class="choice-card">
          <div class="choice-header">
            <div class="choice-icon">2</div>
            <div>
              <h3>2. Predefined Benchmarks</h3>
              <p class="choice-subtitle">1-Click classic operations research problems</p>
            </div>
          </div>
          <p class="choice-desc">
            Test Bland's smallest-subscript rule on classic textbook problems and observe anti-cycling:
          </p>

          <div style="margin-top: 14px; display: flex; flex-direction: column; gap: 10px;">
            <button type="button" class="benchmark-btn" onclick="loadDemo('beale')">
              <div class="bench-title">Beale's Cycling Problem (1955)</div>
              <div class="bench-desc">Dantzig rule cycles indefinitely with period 6; Bland's rule terminates in 7 steps!</div>
            </button>

            <button type="button" class="benchmark-btn" onclick="loadDemo('twophase')">
              <div class="bench-title">Two-Phase Problem (y Artificials)</div>
              <div class="bench-desc">Mixed &ge;, =, &le; constraints with artificial variables y in yellow. Phase 1 &rarr; Phase 2.</div>
            </button>

            <button type="button" class="benchmark-btn" onclick="loadDemo('standard2d')">
              <div class="bench-title">Standard 2D Textbook LP</div>
              <div class="bench-desc">Min -3x1 - 5x2 with standard inequalities &le;. Pedagogical example.</div>
            </button>
          </div>
        </div>

        <!-- Option 3: Manual LP Input -->
        <div class="choice-card">
          <div class="choice-header">
            <div class="choice-icon">3</div>
            <div>
              <h3>3. Input LP Manually</h3>
              <p class="choice-subtitle">Formulate custom equations from scratch</p>
            </div>
          </div>
          <p class="choice-desc">
            Specify the objective function coefficients (minimizing z), number of decision variables, and custom constraints with &le;, &ge;, or = relations.
          </p>

          <div style="margin-top: 14px; background: #f8fafc; border: 1px solid var(--border); border-radius: 8px; padding: 12px; font-size: 13px; color: var(--text-muted);">
            Automatically assigns artificial variables <code>y1, y2, ...</code> with Two-Phase Simplex if &ge; or = constraints exist.
          </div>

          <div style="margin-top: 24px;">
            <button type="button" class="btn" style="width: 100%;" onclick="openManualInputTab()">
              Open Manual LP Form
            </button>
          </div>
        </div>

      </div>
    </div>

    <!-- TAB 2: STEP-BY-STEP TABLEAU VISUALIZER -->
    <div id="tab-visualizer" class="tab-content">
      <div class="card">
        
        <!-- CONSOLIDATED NAVIGATION BOX: Iteration Row + Sub-step Row -->
        <div class="stepper-box">
          <!-- Row 1: Iteration Controls (Centered) -->
          <div class="stepper-row">
            <div class="iteration-side-left">
              <button id="play-btn" class="btn btn-secondary" onclick="togglePlay()">Auto-play</button>
            </div>

            <div class="iteration-controls">
              <button id="btn-first-step" class="btn btn-secondary" onclick="firstStep()" title="First step (Initial tableau)">First</button>
              <button id="btn-prev-iter" class="btn btn-secondary" onclick="prevIteration()" title="Previous iteration">&larr; Prev Iteration</button>
              <span id="iter-label" class="iter-counter">Iteration 0</span>
              <button id="btn-next-iter" class="btn btn-secondary" onclick="nextIteration()" title="Next iteration">Next Iteration &rarr;</button>
              <button id="btn-last-step" class="btn btn-secondary" onclick="lastStep()" title="Last step (Optimal tableau)">Last</button>
            </div>

            <div class="iteration-side-right">
              <span id="phase-badge" class="badge badge-purple" style="display: none;">Phase</span>
              <span id="status-badge" class="badge badge-blue">In progress</span>
            </div>
          </div>

          <!-- Row 2: Sub-step Controls (Directly below iteration in the same box) -->
          <div class="stepper-row stepper-row-substep">
            <div class="substep-controls">
              <button id="btn-prev-substep" class="btn btn-secondary" onclick="prevStep()" title="Previous elementary sub-step">&larr; Prev Sub-step</button>
              <span id="substep-counter-label" class="substep-counter-badge">Sub-step 1 of 1</span>
              <button id="btn-next-substep" class="btn btn-secondary" onclick="nextStep()" title="Next elementary sub-step">Next Sub-step &rarr;</button>
            </div>
          </div>
        </div>

        <!-- Sub-step Elementary Row Operation Banner Container -->
        <div id="substep-banner-container"></div>

        <!-- Prominent Phase Banner Container -->
        <div id="phase-banner-container"></div>

        <div style="margin-bottom: 12px; display: flex; justify-content: space-between; align-items: center;">
          <h2 id="tableau-title" style="font-size: 17px; font-weight: 700;">Initial Tableau</h2>
          <div id="tableau-obj" style="font-size: 14px; font-weight: 700; color: var(--purple);">-z = 0</div>
        </div>

        <!-- Tableau Visual Legend: Variable Basis Status & Active Step Highlights -->
        <div class="tableau-legend">
          <div class="legend-item"><span class="col-tag tag-basic">Basic</span> <span><strong>Basic Variable</strong> (in basis &bull; Red)</span></div>
          <div class="legend-item"><span class="col-tag tag-nonbasic">Non-Basic</span> <span><strong>Non-Basic Variable</strong> (&bull; Blue)</span></div>
          <div class="legend-item"><span class="legend-swatch swatch-pivot"></span> <span><strong>Pivot Element</strong> (Gold)</span></div>
          <div class="legend-item"><span class="legend-swatch swatch-active-row"></span> <span><strong>Active Operation Row</strong></span></div>
          <div class="legend-item"><span class="legend-swatch swatch-elim"></span> <span><strong>Eliminated Entry &rarr; 0</strong> (Green)</span></div>
        </div>

        <!-- Main Tableau Table Container -->
        <div class="table-container">
          <table class="tableau" id="tableau-table">
            <!-- Rendered dynamically: Basis | RHS | var_names... (No z column) -->
          </table>
        </div>

        <!-- Current Basic Solution Banner -->
        <div style="background: #f8fafc; border: 1px solid var(--border); border-radius: 8px; padding: 10px 14px; font-size: 13px; margin-bottom: 16px;">
          <strong>Current Basic Solution:</strong>
          <span id="basic-solution-text" style="color: #334155; margin-left: 8px;"></span>
        </div>

        <!-- Bland's Rule Detailed Step Breakdown -->
        <div class="analysis-grid">
          
          <!-- Card 1: Entering Variable -->
          <div class="info-card">
            <h4>1. Entering Variable Selection (Bland's Rule)</h4>
            <div id="entering-analysis" style="font-size: 13px;"></div>
          </div>

          <!-- Card 2: Minimum Theta Test & Leaving Variable -->
          <div class="info-card">
            <h4>2. Minimum &theta; (Theta) Test & Leaving Variable (Bland)</h4>
            <div id="leaving-analysis" style="font-size: 13px;"></div>
          </div>

          <!-- Card 3: Pivot Element & Row Operations -->
          <div class="info-card">
            <h4>3. Row Operations & Pivot</h4>
            <div id="pivot-analysis" style="font-size: 13px;"></div>
          </div>

        </div>

        <!-- Final Result Card (shown on optimal) -->
        <div id="final-result-card" class="result-banner" style="display: none;">
          <h3>Optimal Solution Reached!</h3>
          <div id="final-result-body" style="font-size: 14px;"></div>
        </div>

      </div>
    </div>

    <!-- TAB 3: CSV MANAGER & EDITOR -->
    <div id="tab-csv" class="tab-content">
      <div class="card">
        <div class="card-title">
          <span>Tableau CSV Manager & Interactive Spreadsheet Editor</span>
          <div class="btn-group">
            <button class="btn btn-secondary" onclick="switchTab('choose')">Change Problem</button>
            <button class="btn" onclick="reloadCsvList()">Refresh List</button>
          </div>
        </div>

        <div class="form-row">
          <label style="font-weight: 600; font-size: 14px;">Detected CSV Files:</label>
          <select id="csv-select" class="select-input" style="max-width: 320px;" onchange="loadSelectedCsv()">
            <option value="tableau.csv" selected>tableau.csv</option>
            <option value="example_standard.csv">example_standard.csv</option>
          </select>
          <button class="btn" onclick="solveCurrentCsv()">Solve This CSV</button>
          <label class="toggle-card active" for="csv-force-twophase" style="margin-left: 6px;">
            <span class="toggle-switch">
              <input type="checkbox" id="csv-force-twophase" checked onchange="this.closest('.toggle-card').classList.toggle('active', this.checked)">
              <span class="toggle-slider"></span>
            </span>
            <span class="toggle-content">
              <span class="toggle-title">Force Two-Phase Simplex</span>
            </span>
          </label>
        </div>

        <div style="margin: 16px 0; border-top: 1px solid var(--border); padding-top: 16px;">
          <h4 style="font-size: 14px; font-weight: 700; margin-bottom: 8px;">Or Upload a CSV File from your computer:</h4>
          <input type="file" id="file-upload" accept=".csv" class="select-input" style="max-width: 400px;" onchange="handleFileUpload(event)">
        </div>

        <!-- Editable CSV Grid -->
        <div style="margin-top: 20px;">
          <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px; flex-wrap: wrap; gap: 8px;">
            <h4 style="font-size: 14px; font-weight: 700;">Tableau Preview & Interactive Editor (Basis, RHS, variables...):</h4>
            <div class="btn-group">
              <button class="btn btn-secondary" onclick="addRow()">+ Add Row</button>
              <button class="btn btn-secondary" onclick="saveCsvToDisk()">Save Changes to CSV</button>
              <button class="btn" onclick="solveEditedTableau()">Solve Modified Tableau</button>
            </div>
          </div>
          <div class="table-container">
            <table id="editor-table" class="editor-table">
              <!-- Rendered dynamically -->
            </table>
          </div>
        </div>

      </div>
    </div>

    <!-- TAB 4: MANUAL LP INPUT -->
    <div id="tab-manual" class="tab-content">
      <div class="card">
        <div class="card-title">
          <span>Manual Linear Programming (LP) Input</span>
          <button class="btn btn-secondary" onclick="switchTab('choose')">Change Problem</button>
        </div>
        <p style="font-size: 13px; color: var(--text-muted); margin-bottom: 16px;">
          Enter the problem specification. The program will automatically construct the initial tableau
          (including Two-Phase Phase 1 with artificial variables y1, y2, ... in yellow if &ge; or = constraints are present) and solve it with Bland's Rule.
        </p>

        <div class="form-row" style="align-items: flex-start;">
          <div class="form-group" style="min-width: 190px;">
            <label style="font-weight: 600; font-size: 13px;">Optimization Goal:</label>
            <select id="manual-obj-type" class="select-input" disabled>
              <option value="min" selected>Minimization (MIN Z)</option>
            </select>
          </div>

          <div class="form-group" style="width: 120px;">
            <label style="font-weight: 600; font-size: 13px;">Variables (n):</label>
            <input type="number" id="manual-n-vars" class="select-input" value="2" min="1" max="10" onchange="buildManualForm()">
          </div>

          <div class="form-group" style="width: 120px;">
            <label style="font-weight: 600; font-size: 13px;">Constraints (m):</label>
            <input type="number" id="manual-n-constraints" class="select-input" value="3" min="1" max="15" onchange="buildManualForm()">
          </div>

          <div class="form-group" style="flex: 1; min-width: 290px;">
            <label style="font-weight: 600; font-size: 13px;">Algorithm Mode:</label>
            <label class="toggle-card" for="manual-force-twophase">
              <span class="toggle-switch">
                <input type="checkbox" id="manual-force-twophase" onchange="this.closest('.toggle-card').classList.toggle('active', this.checked)">
                <span class="toggle-slider"></span>
              </span>
              <span class="toggle-content">
                <span class="toggle-title">Force Two-Phase Simplex</span>
                <span class="toggle-desc">Introduce artificial variables (y) for all constraints</span>
              </span>
            </label>
          </div>
        </div>

        <div id="manual-form-container" style="margin-top: 16px;">
          <!-- Generated inputs -->
        </div>

        <button class="btn" style="margin-top: 20px; padding: 10px 20px; font-size: 14px;" onclick="solveManualLP()">
          Solve with Simplex Method (Bland's Rule)
        </button>
      </div>
    </div>

    <!-- TAB 5: DEMO PROBLEMS -->
    <div id="tab-demos" class="tab-content">
      <div class="card">
        <div class="card-title">
          <span>Benchmark & Demonstration Problems (Bland's Rule)</span>
          <button class="btn btn-secondary" onclick="switchTab('choose')">Change Problem</button>
        </div>
        <div class="analysis-grid">
          
          <div class="info-card">
            <h4>Beale's Cycling Problem (1955)</h4>
            <p style="font-size: 13px; color: var(--text-muted); margin-bottom: 10px;">
              The quintessential benchmark problem: under Dantzig's largest-coefficient rule, it cycles indefinitely with period 6.
              With <strong>Bland's Rule</strong>, ties are resolved and it converges to the optimum in 7 iterations!
            </p>
            <button class="btn" onclick="loadDemo('beale')">Load & Solve Beale</button>
          </div>

          <div class="info-card">
            <h4>Standard 2D Problem</h4>
            <p style="font-size: 13px; color: var(--text-muted); margin-bottom: 10px;">
              Min -3x1 - 5x2 subject to standard inequality constraints &le;. Clean pedagogical step-by-step example.
            </p>
            <button class="btn" onclick="loadDemo('standard2d')">Load & Solve 2D</button>
          </div>

          <div class="info-card">
            <h4>Two-Phase Problem (Constraints =, &ge;, &le;)</h4>
            <p style="font-size: 13px; color: var(--text-muted); margin-bottom: 10px;">
              Min 4x1 + x2 with mixed constraints forcing artificial variables <strong>y1, y2, y3</strong> in yellow.
              Demonstrates Phase 1 (-w = rhs &rarr; 0) and Phase 2 (-z = rhs &rarr; -17/5, Min z = 17/5 = 3.4).
            </p>
            <button class="btn" onclick="loadDemo('twophase')">Load & Solve Two-Phase</button>
          </div>

        </div>
      </div>
    </div>

  </div>

  <script>
    let stepsData = [];
    let currentStepIdx = 0;
    let playInterval = null;
    let editorHeaders = [];
    let editorRows = [];

    function switchTab(tabId) {
      document.querySelectorAll('.tab-btn').forEach(btn => btn.classList.remove('active'));
      document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
      const activeBtn = document.getElementById(`tab-btn-${tabId}`);
      if (activeBtn) activeBtn.classList.add('active');
      const activeContent = document.getElementById(`tab-${tabId}`);
      if (activeContent) activeContent.classList.add('active');
    }

    function syncLandingCsv(val) {
      const select = document.getElementById('csv-select');
      if (select) select.value = val;
    }

    // Load CSV list on init
    async function reloadCsvList() {
      try {
        const res = await fetch('/api/csv-files');
        const files = await res.json();
        
        const selects = [document.getElementById('csv-select'), document.getElementById('landing-csv-select')];
        selects.forEach(select => {
          if (!select) return;
          select.innerHTML = '';
          files.forEach(f => {
            const opt = document.createElement('option');
            opt.value = f;
            opt.textContent = f;
            select.appendChild(opt);
          });
          if (files.includes('tableau.csv')) {
            select.value = 'tableau.csv';
          } else if (files.length > 0) {
            select.value = files[0];
          }
        });

        const activeSelect = document.getElementById('csv-select');
        if (activeSelect && activeSelect.value) {
          await loadSelectedCsv();
        }
      } catch (err) {
        console.error('Error fetching CSV list:', err);
      }
    }

    async function solveLandingCsv() {
      const select = document.getElementById('landing-csv-select');
      const filename = (select && select.value) ? select.value : 'tableau.csv';
      await solveCsvFilename(filename);
    }

    async function previewLandingCsv() {
      const select = document.getElementById('landing-csv-select');
      const filename = (select && select.value) ? select.value : 'tableau.csv';
      const csvSelect = document.getElementById('csv-select');
      if (csvSelect) csvSelect.value = filename;
      await loadSelectedCsv();
      switchTab('csv');
    }

    function openManualInputTab() {
      buildManualForm();
      switchTab('manual');
    }

    async function loadSelectedCsv() {
      const select = document.getElementById('csv-select');
      if (!select.value) return;
      try {
        const res = await fetch('/api/load-csv?file=' + encodeURIComponent(select.value));
        const data = await res.json();
        if (data.error) {
          alert('Error loading CSV: ' + data.error);
          return;
        }
        editorHeaders = data.headers;
        editorRows = data.rows;
        renderEditorTable();
      } catch (err) {
        alert('Error loading CSV: ' + err);
      }
    }

    function renderEditorTable() {
      const table = document.getElementById('editor-table');
      table.innerHTML = '';
      
      // Header row
      const thead = document.createElement('thead');
      const trH = document.createElement('tr');
      const thRow = document.createElement('th');
      thRow.style.width = '60px';
      thRow.style.textAlign = 'center';
      thRow.textContent = 'Row';
      trH.appendChild(thRow);
      editorHeaders.forEach((h, cIdx) => {
        const th = document.createElement('th');
        if (cIdx === 1) th.className = 'col-rhs';
        th.innerHTML = `<input type="text" value="${h}" onchange="editorHeaders[${cIdx}] = this.value">`;
        trH.appendChild(th);
      });
      thead.appendChild(trH);
      table.appendChild(thead);

      // Body rows
      const tbody = document.createElement('tbody');
      editorRows.forEach((row, rIdx) => {
        const tr = document.createElement('tr');
        if (rIdx === 0) tr.className = 'cost-row';
        const tdRow = document.createElement('td');
        tdRow.style.textAlign = 'center';
        tdRow.innerHTML = `<span class="row-badge">R${rIdx}</span>`;
        tr.appendChild(tdRow);
        row.forEach((cell, cIdx) => {
          const td = document.createElement('td');
          if (cIdx === 1) td.className = 'cell-rhs';
          td.innerHTML = `<input type="text" value="${cell}" onchange="editorRows[${rIdx}][${cIdx}] = this.value">`;
          tr.appendChild(td);
        });
        tbody.appendChild(tr);
      });
      table.appendChild(tbody);
    }

    function addRow() {
      const newRow = new Array(editorHeaders.length).fill('0');
      let maxX = 0;
      editorHeaders.forEach(h => {
        if (h && h.startsWith('x') && !isNaN(h.slice(1))) {
          maxX = Math.max(maxX, parseInt(h.slice(1), 10));
        }
      });
      newRow[0] = 'x' + (maxX + 1);
      newRow[1] = '0'; // RHS column
      editorRows.push(newRow);
      renderEditorTable();
    }

    async function saveCsvToDisk() {
      const filename = document.getElementById('csv-select').value || 'tableau.csv';
      const nl = String.fromCharCode(10);
      const content = [editorHeaders.join(',')].concat(editorRows.map(r => r.join(','))).join(nl);
      try {
        const res = await fetch('/api/save-csv', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ filename, content })
        });
        const data = await res.json();
        alert(data.message || ('File ' + filename + ' saved successfully!'));
        await reloadCsvList();
      } catch (e) {
        alert('Error saving file: ' + e);
      }
    }

    async function solveCsvFilename(filename) {
      try {
        const forceTp = document.getElementById('csv-force-twophase') ? document.getElementById('csv-force-twophase').checked : true;
        const res = await fetch('/api/solve-file?file=' + encodeURIComponent(filename) + '&force_two_phase=' + forceTp);
        const text = await res.text();
        let data;
        try {
          data = JSON.parse(text);
        } catch (e) {
          console.error("Non-JSON response:", text);
          alert("Server response error: " + text.slice(0, 150));
          return;
        }
        if (data.error) {
          alert("Error solving: " + data.error);
          return;
        }
        setSolutionData(data);
        switchTab('visualizer');
      } catch (err) {
        alert('Error solving: ' + err);
      }
    }

    async function solveCurrentCsv() {
      const select = document.getElementById('csv-select');
      if (!select.value) return;
      await solveCsvFilename(select.value);
    }

    async function solveEditedTableau() {
      try {
        const forceTp = document.getElementById('csv-force-twophase') ? document.getElementById('csv-force-twophase').checked : true;
        const res = await fetch('/api/solve-tableau', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ headers: editorHeaders, rows: editorRows, force_two_phase: forceTp })
        });
        const text = await res.text();
        let data;
        try {
          data = JSON.parse(text);
        } catch (e) {
          console.error("Non-JSON response:", text);
          alert("Server response error: " + text.slice(0, 150));
          return;
        }
        if (data.error) {
          alert("Error solving: " + data.error);
          return;
        }
        setSolutionData(data);
        switchTab('visualizer');
      } catch (err) {
        alert('Error solving: ' + err);
      }
    }

    function setSolutionData(data) {
      if (!data || !data.steps || data.steps.length === 0) {
        alert('No solution steps were returned.');
        return;
      }
      stepsData = data.steps;
      currentStepIdx = 0;
      renderCurrentStep();
    }

    function renderCurrentStep() {
      if (!stepsData || stepsData.length === 0) return;
      const step = stepsData[currentStepIdx];
      const total = stepsData.length;

      const iterLabelText = step.iteration > 0 ? `Iteration ${step.iteration}` : 'Initial Tableau';
      const iterEl = document.getElementById('iter-label');
      if (iterEl) iterEl.textContent = iterLabelText;
      document.getElementById('tableau-title').textContent = step.title;
      const isMin = true;
      const cleanName = step.obj_name.replace('-', '').toUpperCase();

      document.getElementById('tableau-obj').innerHTML =
        `<span style="font-size: 13px; color: #64748b; font-weight: normal; margin-right: 8px;">${step.obj_name} = ${step.raw_rhs} &nbsp;&rarr;</span><strong>MIN ${cleanName} = ${step.obj_val}</strong>`;

      // Stepper Bar Phase Badges (Only displayed for Two-Phase Active Phases)
      const phaseBadge = document.getElementById('phase-badge');
      if (step.is_phase1) {
        phaseBadge.style.display = 'inline-block';
        phaseBadge.textContent = 'Phase 1 (-w = rhs)';
        phaseBadge.className = 'badge badge-amber';
      } else if (step.is_phase2) {
        phaseBadge.style.display = 'inline-block';
        phaseBadge.textContent = `Phase 2 (${step.obj_name} = rhs)`;
        phaseBadge.className = 'badge badge-purple';
      } else {
        phaseBadge.style.display = 'none';
      }

      // Status badge
      const badge = document.getElementById('status-badge');
      if (step.status === 'OPTIMAL') {
        badge.textContent = 'Optimal Reached';
        badge.className = 'badge badge-green';
      } else if (step.status === 'UNBOUNDED') {
        badge.textContent = 'Unbounded Problem';
        badge.className = 'badge badge-red';
      } else {
        badge.textContent = `Iteration ${step.iteration}`;
        badge.className = 'badge badge-blue';
      }

      // Navigation button disabled states in consolidated stepper box
      const btnFirst = document.getElementById('btn-first-step');
      if (btnFirst) btnFirst.disabled = (currentStepIdx === 0);
      const btnPrevIter = document.getElementById('btn-prev-iter');
      if (btnPrevIter) btnPrevIter.disabled = (currentStepIdx === 0);
      const btnNextIter = document.getElementById('btn-next-iter');
      if (btnNextIter) btnNextIter.disabled = (currentStepIdx === total - 1);
      const btnLast = document.getElementById('btn-last-step');
      if (btnLast) btnLast.disabled = (currentStepIdx === total - 1);

      // Sub-step navigation in consolidated stepper box (Row 2)
      const substepCounter = document.getElementById('substep-counter-label');
      if (substepCounter) {
        substepCounter.textContent = `Sub-step ${step.substep_idx} of ${step.substep_total} • (Step ${currentStepIdx + 1} of ${total})`;
      }
      const prevSubBtn = document.getElementById('btn-prev-substep');
      if (prevSubBtn) prevSubBtn.disabled = (currentStepIdx === 0);
      const nextSubBtn = document.getElementById('btn-next-substep');
      if (nextSubBtn) nextSubBtn.disabled = (currentStepIdx === total - 1);

      // Prominent Sub-step Operation Banner (Below Navigation Box)
      const substepContainer = document.getElementById('substep-banner-container');
      if (substepContainer) {
        let pillCls = 'substep-pill-pivot';
        let pillText = `Step ${step.substep_idx} of ${step.substep_total}`;
        if (step.substep_type === 'normalize_pivot') {
          pillCls = 'substep-pill-norm';
          pillText = `Normalize Pivot Row`;
        } else if (step.substep_type === 'eliminate_row') {
          pillCls = 'substep-pill-elim';
          pillText = `Eliminate Row R${step.target_row}`;
        } else if (step.substep_type === 'iteration_complete') {
          pillCls = 'substep-pill-done';
          pillText = `Iteration ${step.iteration} Complete`;
        } else if (step.substep_type === 'optimal') {
          pillCls = 'substep-pill-opt';
          pillText = `Optimal Solution Reached`;
        } else if (step.substep_type === 'unbounded') {
          pillCls = 'substep-pill-unb';
          pillText = `Problem Unbounded`;
        } else {
          pillCls = 'substep-pill-pivot';
          pillText = `Pivot Selection (Bland & θ)`;
        }

        const formulaHtml = step.formula ? `<div class="substep-formula">${step.formula}</div>` : '';
        const descHtml = step.substep_desc ? `<div class="substep-desc">${step.substep_desc}</div>` : '';

        substepContainer.innerHTML = `
          <div class="substep-banner">
            <div class="substep-banner-header">
              <span class="substep-pill ${pillCls}">${pillText}</span>
            </div>
            ${formulaHtml}
            ${descHtml}
          </div>
        `;
      }

      // Phase Banner rendering (Only displayed for Two-Phase Active Phases)
      const bannerContainer = document.getElementById('phase-banner-container');
      if (step.is_phase1) {
        bannerContainer.innerHTML = `
          <div class="phase-banner phase-banner-p1">
            <span class="phase-badge-p1">PHASE 1 ACTIVE</span>
            <div>
              <strong>Phase 1: Finding an Initial Feasible Basis</strong> &mdash;
              Minimizing the sum of artificial variables (<code>y1, y2, ...</code>).
              Tableau row 0 tracks <code>-w = rhs</code>. Goal: Drive <code>-w &rarr; 0</code> so all artificial variables leave the basis.
            </div>
          </div>`;
      } else if (step.is_phase2) {
        bannerContainer.innerHTML = `
          <div class="phase-banner phase-banner-p2">
            <span class="phase-badge-p2">PHASE 2 ACTIVE</span>
            <div>
              <strong>Phase 2: Optimizing Original Objective</strong> &mdash;
              Feasible basis achieved (<code>-w = 0</code>)! All artificial variables eliminated.
              Now optimizing original objective with tableau row 0 tracking <code>${step.obj_name} = rhs</code>.
            </div>
          </div>`;
      } else {
        bannerContainer.innerHTML = '';
      }

      // Basic solution text
      const solParts = [];
      for (const [k, v] of Object.entries(step.basic_solution)) {
        if (k.startsWith("-")) {
          const posK = k.replace('-', '');
          solParts.push(`<strong>${k}</strong> = ${v} &nbsp;(<strong style="color: #7c3aed;">${posK} = ${step.obj_val}</strong>)`);
        } else {
          solParts.push(`<strong>${k}</strong> = ${v}`);
        }
      }
      document.getElementById('basic-solution-text').innerHTML = solParts.join(', ');

      // Build tableau HTML
      const table = document.getElementById('tableau-table');
      table.innerHTML = '';

      // Determine pivot row and col indices
      const pivotRow = step.pivot ? step.pivot.row : -1;
      const pivotCol = step.pivot ? step.pivot.col : -1;
      const isOpSubstep = (step.substep_type === 'normalize_pivot' || step.substep_type === 'eliminate_row');
      const isPivotSubstep = (step.substep_type === 'pivot_selection');

      // Table Header: Row | Basis | RHS | var_names... (NO z column)
      const thead = document.createElement('thead');
      const trH = document.createElement('tr');
      trH.innerHTML = `
        <th class="col-row-idx">Row</th>
        <th class="col-basis">Basis</th>
        <th class="col-rhs">
          <div class="col-header-box">
            <span>RHS</span>
            <span class="col-tag tag-rhs">RHS</span>
          </div>
        </th>`;

      step.var_names.forEach((name, cIdx) => {
        const isPivotCol = (cIdx === pivotCol);
        const isBasic = step.basic_vars && step.basic_vars.includes(name);

        const nameCls = isBasic ? 'col-name-basic' : 'col-name-nonbasic';
        const tagCls = isBasic ? 'tag-basic' : 'tag-nonbasic';
        const tagText = isBasic ? 'Basic' : 'Non-Basic';

        const isPivotColActive = (isPivotSubstep || isOpSubstep) && isPivotCol;
        const pivotHeaderCls = isPivotColActive ? ' col-pivot-active' : '';
        const pivotTagHtml = isPivotColActive ? '<span class="col-tag tag-pivot">Pivot Col</span>' : '';

        trH.innerHTML += `
          <th class="${pivotHeaderCls}">
            <div class="col-header-box">
              <span class="${nameCls}">${name}</span>
              <span class="col-tag ${tagCls}">${tagText}</span>
              ${pivotTagHtml}
            </div>
          </th>`;
      });
      thead.appendChild(trH);
      table.appendChild(thead);

      // Table Body
      const tbody = document.createElement('tbody');
      step.rows.forEach((row, rIdx) => {
        const tr = document.createElement('tr');
        const isPivotRow = (rIdx === pivotRow);
        const isTargetRow = (rIdx === step.target_row);
        const isActiveRow = step.active_rows && step.active_rows.includes(rIdx);
        const isRowDimmed = isOpSubstep && !isActiveRow;

        if (rIdx === 0) {
          tr.className = 'cost-row ' + (step.is_phase1 ? 'cost-row-p1' : 'cost-row-p2');
        }
        if (isPivotRow) tr.classList.add('pivot-row-active');

        // 1. Row Index Cell
        let rowBadgeCls = 'row-badge';
        let rowBadgeTitle = `Row ${rIdx}`;
        let rowBadgeText = `R${rIdx}`;

        if (isOpSubstep) {
          if (step.substep_type === 'normalize_pivot' && isPivotRow) {
            rowBadgeCls = 'row-badge row-badge-pivot';
            rowBadgeTitle = 'Pivot Row being normalized';
            rowBadgeText = `R${rIdx} (Norm)`;
          } else if (step.substep_type === 'eliminate_row') {
            if (isTargetRow) {
              rowBadgeCls = 'row-badge';
              rowBadgeTitle = `Target Row R${rIdx} being eliminated`;
              rowBadgeText = `R${rIdx} (Target)`;
            } else if (isPivotRow) {
              rowBadgeCls = 'row-badge row-badge-pivot';
              rowBadgeTitle = `Pivot Row R${rIdx} reference`;
              rowBadgeText = `R${rIdx} (Ref)`;
            }
          }
        } else if (isPivotRow) {
          rowBadgeCls = 'row-badge row-badge-pivot';
          rowBadgeTitle = 'Pivot Row (Leaving Variable)';
          rowBadgeText = `R${rIdx} (Pivot)`;
        }

        const rowIdxTd = `<td class="cell-row-idx ${isRowDimmed ? 'cell-dimmed' : ''}"><span class="${rowBadgeCls}" title="${rowBadgeTitle}">${rowBadgeText}</span></td>`;

        // 2. Basis Cell
        const isYBasic = row.basic.toLowerCase().startsWith('y');
        let basisBadge = rIdx === 0
          ? `<span class="badge ${step.is_phase1 ? 'badge-amber' : 'badge-purple'}">${row.basic}</span>`
          : (isYBasic ? `<span class="badge badge-amber">${row.basic}</span>` : `<span class="badge badge-red">${row.basic}</span>`);
        if (isRowDimmed) {
          basisBadge = `<span style="color: #94a3b8; font-weight: 600;">${row.basic}</span>`;
        }
        const basisTd = `<td class="cell-basis ${isRowDimmed ? 'cell-dimmed' : ''}">${basisBadge}</td>`;

        // 3. RHS Cell
        let rhsCls = 'cell-rhs';
        if (isRowDimmed) {
          rhsCls += ' cell-dimmed';
        } else if (step.substep_type === 'normalize_pivot' && isPivotRow) {
          rhsCls += ' cell-active-norm';
        } else if (step.substep_type === 'eliminate_row') {
          if (isTargetRow) rhsCls += ' cell-target-row';
          else if (isPivotRow) rhsCls += ' cell-pivot-row';
        }
        const rhsTd = `<td class="${rhsCls}">${row.rhs}</td>`;

        tr.innerHTML = rowIdxTd + basisTd + rhsTd;

        // 4. Variable Cells (Selective focus during operations; clean 1-color-per-cell otherwise)
        row.cells.forEach((cellVal, cIdx) => {
          const varName = step.var_names[cIdx];
          const isPivotCell = (rIdx === pivotRow && cIdx === pivotCol);
          const isElimCell = (step.eliminated_cell && step.eliminated_cell.row === rIdx && step.eliminated_cell.col === cIdx);

          let cellCls = '';
          let titleText = '';

          if (isOpSubstep) {
            if (isRowDimmed) {
              cellCls = 'cell-dimmed';
              titleText = `Row R${rIdx} (inactive during this row operation)`;
            } else if (step.substep_type === 'normalize_pivot') {
              if (isPivotCell) {
                cellCls = 'cell-pivot';
                titleText = `Pivot Element = 1 (Normalized leading 1)`;
              } else {
                cellCls = 'cell-active-norm';
                titleText = `Normalized: R${rIdx} / pivot = ${cellVal}`;
              }
            } else if (step.substep_type === 'eliminate_row') {
              if (isElimCell) {
                cellCls = 'cell-eliminated';
                titleText = `Eliminated: a(R${rIdx}, ${varName}) driven to 0`;
              } else if (isTargetRow) {
                cellCls = 'cell-target-row';
                titleText = `Updated Row R${rIdx} entry: ${cellVal}`;
              } else if (isPivotRow) {
                if (isPivotCell) {
                  cellCls = 'cell-pivot';
                  titleText = `Pivot Reference Element = 1`;
                } else {
                  cellCls = 'cell-pivot-row';
                  titleText = `Pivot Reference Row R${rIdx} entry: ${cellVal}`;
                }
              }
            }
          } else if (isPivotSubstep) {
            // Pivot Selection Step: highlight ONLY the cells/rows/columns used in this step
            if (isPivotCell) {
              cellCls = 'cell-pivot';
              titleText = `Pivot Element: a(R${rIdx}, ${varName}) = ${cellVal}`;
            } else if (rIdx === pivotRow) {
              cellCls = 'cell-pivot-row';
              titleText = `Leaving Variable Row R${rIdx}`;
            } else if (cIdx === pivotCol) {
              cellCls = 'cell-target-row';
              titleText = `Entering Variable Column ${varName}`;
            }
            // All other cells: clean neutral (no color)
          }
          // In iteration_complete, optimal, unbounded, or initial: cellCls remains '' (NO COLORS)!

          tr.innerHTML += `<td class="${cellCls}" title="${titleText}">${cellVal}</td>`;
        });

        tbody.appendChild(tr);
      });
      table.appendChild(tbody);

      // Analysis Cards
      renderAnalysis(step, isMin, cleanName);
    }

    function renderAnalysis(step, isMin, cleanName) {
      const enteringDiv = document.getElementById('entering-analysis');
      const leavingDiv = document.getElementById('leaving-analysis');
      const pivotDiv = document.getElementById('pivot-analysis');
      const finalCard = document.getElementById('final-result-card');

      if (step.status === 'OPTIMAL') {
        enteringDiv.innerHTML = '<span class="badge badge-green">No negative reduced cost</span><p style="margin-top: 8px;">All coefficients in Row 0 are &ge; 0. The tableau is optimal!</p>';
        leavingDiv.innerHTML = 'No variable leaves the basis.';
        pivotDiv.innerHTML = 'Optimal solution reached &mdash; no row operations needed.';
        
        // Show final result card
        finalCard.style.display = 'block';
        const body = document.getElementById('final-result-body');
        
        let reportHtml = `<div style="font-size: 15px; margin-bottom: 8px;">Tableau Row 0 Equation: <strong>${step.obj_name}* = ${step.raw_rhs}</strong> &nbsp;&rarr;&nbsp; Optimal Minimum Value: <strong style="color: #14532d; font-size: 16px;">MIN ${cleanName}* = -(${step.obj_name}*) = ${step.obj_val}</strong></div>`;

        const basicPairs = [];
        for (const [k, v] of Object.entries(step.basic_solution)) {
          if (!k.startsWith("-") && k !== "z" && k !== "w") {
            basicPairs.push(`<code>${k} = ${v}</code>`);
          }
        }
        reportHtml += `<div>Optimal Basic Variables: ${basicPairs.join(', ') || 'None'}</div>`;
        body.innerHTML = reportHtml;

      } else if (step.status === 'UNBOUNDED') {
        enteringDiv.innerHTML = `<span class="badge badge-red">Unbounded Problem</span><p style="margin-top: 8px;">${step.message}</p>`;
        leavingDiv.innerHTML = 'No valid &theta; ratio found (all column entries &le; 0).';
        pivotDiv.innerHTML = 'None';
        finalCard.style.display = 'none';

      } else {
        finalCard.style.display = 'none';

        // 1. Entering analysis
        let entHtml = '';
        if (step.candidates && step.candidates.length > 0) {
          entHtml += `<div style="margin-bottom: 6px; color: var(--text-muted);">Improving candidates (Row 0 coeff &lt; 0):</div>`;
          entHtml += '<ul style="margin-bottom: 8px;">';
          step.candidates.forEach(c => {
            const isChosen = step.entering_var && c.name === step.entering_var.name;
            const chosenBadge = isChosen ? ' <span class="badge badge-blue">Selected by Bland</span>' : '';
            entHtml += `<li>&bull; <strong>${c.name}</strong>: coeff = <code>${c.coeff}</code>, Bland index = <strong>${c.index}</strong>${chosenBadge}</li>`;
          });
          entHtml += '</ul>';
          if (step.entering_var) {
            entHtml += `<div><strong>Bland's Rule:</strong> Chooses <strong>${step.entering_var.name}</strong> (smallest Bland index ${step.entering_var.index}).</div>`;
          }
        } else {
          entHtml = 'No candidates with negative reduced cost.';
        }
        enteringDiv.innerHTML = entHtml;

        // 2. Leaving analysis (Theta minimum ratio test with [RHS] / [Coeff])
        let leavHtml = '';
        if (step.ratios && step.ratios.length > 0) {
          leavHtml += `<table style="width: 100%; font-size: 12px; border-collapse: collapse; margin-bottom: 8px;">`;
          leavHtml += `<tr style="border-bottom: 1.5px solid var(--border); text-align: left;">
            <th style="padding: 5px 4px;">Row</th>
            <th style="padding: 5px 4px;">Basis</th>
            <th style="padding: 5px 4px;">[RHS] / [Coeff]</th>
            <th style="padding: 5px 4px;">&theta; (Theta)</th>
          </tr>`;
          step.ratios.forEach(r => {
            const isLeaving = step.leaving_var && r.row === step.leaving_var.row;
            const bg = isLeaving ? 'style="background: #fef3c7; font-weight: bold;"' : '';
            if (r.valid) {
              const thetaVal = r.theta || r.ratio;
              leavHtml += `<tr ${bg}>
                <td style="padding: 5px 4px;">R${r.row}</td>
                <td style="padding: 5px 4px;">${r.basic} (${r.index})</td>
                <td style="padding: 5px 4px;"><strong>[${r.rhs}]</strong> / <strong>[${r.coeff}]</strong></td>
                <td style="padding: 5px 4px;">&theta; = <strong>${thetaVal}</strong></td>
              </tr>`;
            } else {
              leavHtml += `<tr>
                <td style="padding: 5px 4px; color: #94a3b8;">R${r.row}</td>
                <td style="padding: 5px 4px; color: #94a3b8;">${r.basic}</td>
                <td style="padding: 5px 4px; color: #94a3b8;">[${r.rhs}] / [${r.coeff}]</td>
                <td style="padding: 5px 4px; color: #94a3b8;">Skipped (${r.reason})</td>
              </tr>`;
            }
          });
          leavHtml += `</table>`;
          if (step.leaving_var) {
            leavHtml += `<div>Minimum ratio <strong>&theta;* = ${step.min_ratio}</strong> &rarr; Leaving variable: <strong>${step.leaving_var.basic}</strong> (Row ${step.leaving_var.row}).</div>`;
            if (step.bland_tie_note) {
              const tieNoteClean = step.bland_tie_note.replace(/ratio/gi, '&theta;');
              leavHtml += `<div style="margin-top: 6px; color: var(--warning); font-size: 12px; background: #fffbeb; padding: 6px 8px; border-radius: 4px; border: 1px solid #fde68a;"><strong>Bland's Tie-Break:</strong> ${tieNoteClean}</div>`;
            }
          }
        }
        leavingDiv.innerHTML = leavHtml;

        // 3. Pivot analysis & Row operations
        let pivHtml = '';
        if (step.pivot) {
          pivHtml += `<div style="margin-bottom: 8px;">Pivot Element at <strong>Row R${step.pivot.row}</strong>, <strong>Col ${step.entering_var ? step.entering_var.name : step.pivot.col}</strong>: <code>${step.pivot.val}</code></div>`;
        }
        if (step.formula) {
          pivHtml += `<div style="margin-bottom: 8px;"><strong>Current Operation:</strong><div class="substep-formula" style="margin-top: 4px;">${step.formula}</div></div>`;
        }
        if (step.substep_desc) {
          pivHtml += `<div style="font-size: 12px; color: #475569;">${step.substep_desc}</div>`;
        }
        pivotDiv.innerHTML = pivHtml;
      }
    }

    function firstStep() { currentStepIdx = 0; renderCurrentStep(); }
    function prevStep() { if (currentStepIdx > 0) { currentStepIdx--; renderCurrentStep(); } }
    function nextStep() { if (currentStepIdx < stepsData.length - 1) { currentStepIdx++; renderCurrentStep(); } }
    function lastStep() { currentStepIdx = stepsData.length - 1; renderCurrentStep(); }

    function nextIteration() {
      if (!stepsData || currentStepIdx >= stepsData.length - 1) return;
      const currentIter = stepsData[currentStepIdx].iteration;
      let targetIdx = -1;
      for (let i = currentStepIdx + 1; i < stepsData.length; i++) {
        if (stepsData[i].iteration > currentIter) {
          targetIdx = i;
          break;
        }
      }
      currentStepIdx = (targetIdx !== -1) ? targetIdx : (stepsData.length - 1);
      renderCurrentStep();
    }

    function prevIteration() {
      if (!stepsData || currentStepIdx <= 0) return;
      const currentIter = stepsData[currentStepIdx].iteration;
      let startOfCurrent = currentStepIdx;
      while (startOfCurrent > 0 && stepsData[startOfCurrent - 1].iteration === currentIter) {
        startOfCurrent--;
      }
      if (currentStepIdx > startOfCurrent) {
        currentStepIdx = startOfCurrent;
      } else {
        let targetIdx = currentStepIdx - 1;
        if (targetIdx >= 0) {
          const prevIter = stepsData[targetIdx].iteration;
          while (targetIdx > 0 && stepsData[targetIdx - 1].iteration === prevIter) {
            targetIdx--;
          }
          currentStepIdx = targetIdx;
        }
      }
      renderCurrentStep();
    }

    function togglePlay() {
      const btn = document.getElementById('play-btn');
      if (playInterval) {
        clearInterval(playInterval);
        playInterval = null;
        btn.textContent = 'Auto-play';
      } else {
        btn.textContent = 'Pause';
        playInterval = setInterval(() => {
          if (currentStepIdx < stepsData.length - 1) {
            nextStep();
          } else {
            togglePlay();
          }
        }, 1200);
      }
    }

    async function loadDemo(demoName) {
      try {
        const res = await fetch('/api/demo?name=' + encodeURIComponent(demoName));
        const text = await res.text();
        let data;
        try {
          data = JSON.parse(text);
        } catch (e) {
          alert("Server error: " + text.slice(0, 150));
          return;
        }
        if (data.error) {
          alert("Demo error: " + data.error);
          return;
        }
        setSolutionData(data);
        switchTab('visualizer');
      } catch (err) {
        alert('Error loading demo: ' + err);
      }
    }

    async function handleLandingFileUpload(e) {
      const file = e.target.files[0];
      if (!file) return;
      const text = await file.text();
      parseCsvTextIntoEditor(file.name, text);
      switchTab('csv');
    }

    async function handleFileUpload(e) {
      const file = e.target.files[0];
      if (!file) return;
      const text = await file.text();
      parseCsvTextIntoEditor(file.name, text);
    }

    function parseCsvTextIntoEditor(filename, text) {
      const nl = String.fromCharCode(10);
      const tabChar = String.fromCharCode(9);
      const lines = text.split(nl).map(l => l.trim()).filter(l => l && !l.startsWith('#'));
      if (lines.length < 2) {
        alert('CSV file is too short.');
        return;
      }
      const delim = lines[0].includes(';') ? ';' : (lines[0].includes(tabChar) ? tabChar : ',');
      editorHeaders = lines[0].split(delim).map(c => c.trim());
      editorRows = lines.slice(1).map(line => line.split(delim).map(c => c.trim()));
      renderEditorTable();
    }

    function buildManualForm() {
      const nVars = parseInt(document.getElementById('manual-n-vars').value) || 2;
      const nConstraints = parseInt(document.getElementById('manual-n-constraints').value) || 3;
      const container = document.getElementById('manual-form-container');
      
      let html = '<div style="margin-bottom: 16px;"><strong>Objective Function:</strong><div class="form-row" style="margin-top: 8px;">';
      for (let j = 1; j <= nVars; j++) {
        html += `<div style="display: flex; align-items: center; gap: 4px;">`;
        html += `<input type="text" id="obj-c-${j}" class="select-input" style="width: 80px;" value="" placeholder="0">`;
        html += `<label style="font-weight: 600;">x${j}</label>`;
        if (j < nVars) html += `<span style="font-weight: 700;">+</span>`;
        html += `</div>`;
      }
      html += '</div></div>';

      html += '<strong>Constraints:</strong>';
      html += '<div style="display: flex; flex-direction: column; gap: 8px; margin-top: 8px;">';
      for (let i = 1; i <= nConstraints; i++) {
        html += `<div class="form-row" style="margin: 0; background: #f8fafc; padding: 8px 12px; border-radius: 6px; border: 1px solid var(--border);">`;
        html += `<span style="font-weight: 700; width: 30px;">(${i})</span>`;
        for (let j = 1; j <= nVars; j++) {
          html += `<div style="display: flex; align-items: center; gap: 4px;">`;
          html += `<input type="text" id="c-${i}-${j}" class="select-input" style="width: 70px;" value="" placeholder="0">`;
          html += `<label style="font-weight: 600;">x${j}</label>`;
          if (j < nVars) html += `<span style="font-weight: 700;">+</span>`;
          html += `</div>`;
        }
        
        html += `<select id="rel-${i}" class="select-input" style="width: 75px;">
                  <option value="<=" selected>&le;</option>
                  <option value=">=">&ge;</option>
                  <option value="=">=</option>
                 </select>`;
        
        html += `<input type="text" id="rhs-${i}" class="select-input" style="width: 80px;" value="" placeholder="0">`;
        html += `</div>`;
      }
      html += '</div>';

      container.innerHTML = html;
    }

    async function solveManualLP() {
      const objType = 'min';
      const nVars = parseInt(document.getElementById('manual-n-vars').value) || 2;
      const nConstraints = parseInt(document.getElementById('manual-n-constraints').value) || 3;
      const forceTwoPhase = document.getElementById('manual-force-twophase').checked;

      const objCoeffs = [];
      for (let j = 1; j <= nVars; j++) {
        objCoeffs.push(document.getElementById(`obj-c-${j}`).value.trim() || '0');
      }

      const constraints = [];
      for (let i = 1; i <= nConstraints; i++) {
        const coeffs = [];
        for (let j = 1; j <= nVars; j++) {
          coeffs.push(document.getElementById(`c-${i}-${j}`).value.trim() || '0');
        }
        const rel = document.getElementById(`rel-${i}`).value;
        const rhs = document.getElementById(`rhs-${i}`).value.trim() || '0';
        constraints.push({ coeffs, rel, rhs });
      }

      const hasAnyNonZero = objCoeffs.some(c => c !== '0') || constraints.some(c => c.coeffs.some(x => x !== '0') || c.rhs !== '0');
      if (!hasAnyNonZero) {
        alert("Please enter the coefficients for your LP model before solving.");
        return;
      }

      try {
        const res = await fetch('/api/solve-lp', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ obj_type: objType, obj_coeffs: objCoeffs, constraints, force_two_phase: forceTwoPhase })
        });
        const text = await res.text();
        let data;
        try {
          data = JSON.parse(text);
        } catch (e) {
          alert("Server response error: " + text.slice(0, 150));
          return;
        }
        if (data.error) {
          alert("Error: " + data.error);
          return;
        }
        setSolutionData(data);
        switchTab('visualizer');
      } catch (err) {
        alert("Error solving LP: " + err);
      }
    }

    // Initialize on load
    window.addEventListener('DOMContentLoaded', async () => {
      await reloadCsvList();
      buildManualForm();
      switchTab('choose'); // Landing screen is active on first open!
    });
  </script>
</body>
</html>
"""

# =====================================================================
# HTTP REQUEST HANDLER & REST API
# =====================================================================

class SimplexGUIRequestHandler(BaseHTTPRequestHandler):
    def log_message(self, format, *args):
        pass

    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        query = urllib.parse.parse_qs(parsed.query)

        if path in ("", "/"):
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.end_headers()
            self.wfile.write(HTML_PAGE.encode("utf-8"))
            return

        elif path == "/api/csv-files":
            csv_files = get_available_csv_files()
            self.send_json(csv_files)
            return

        elif path == "/api/load-csv":
            filename = query.get("file", ["tableau.csv"])[0]
            filepath = resolve_csv_path(filename)
            if not os.path.exists(filepath):
                self.send_json({"error": f"File '{filename}' not found"}, status_code=404)
                return
            try:
                with open(filepath, "r", encoding="utf-8-sig") as f:
                    lines = [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]
                if not lines:
                    self.send_json({"error": "Empty CSV file"}, status_code=400)
                    return
                delim = ';' if ';' in lines[0] else ('\t' if '\t' in lines[0] else ',')
                headers = [c.strip() for c in lines[0].split(delim)]
                rows = [[c.strip() for c in line.split(delim)] for line in lines[1:]]
                # Convert any slack 's' columns to 'x'
                s_cols = {}
                existing_x = [int(h[1:]) for h in headers if h.startswith('x') and h[1:].isdigit()]
                curr_x = max(existing_x, default=0) + 1
                for idx, h in enumerate(headers):
                    if h.lower().startswith('s') and h[1:].isdigit():
                        new_name = f"x{curr_x}"
                        curr_x += 1
                        s_cols[h] = new_name
                        headers[idx] = new_name
                if s_cols:
                    for row in rows:
                        if row and row[0] in s_cols:
                            row[0] = s_cols[row[0]]
                self.send_json({"headers": headers, "rows": rows})
            except Exception as e:
                self.send_json({"error": str(e)}, status_code=500)
            return

        elif path == "/api/solve-file":
            filename = query.get("file", ["tableau.csv"])[0]
            force_tp = query.get("force_two_phase", ["true"])[0].lower() != "false"
            filepath = resolve_csv_path(filename)
            if not os.path.exists(filepath):
                self.send_json({"error": f"File '{filename}' not found"}, status_code=404)
                return
            try:
                tab, detected_obj_type = load_tableau_from_csv(filepath)
                if force_tp:
                    result = record_tableau_steps_twophase(tab, obj_type=detected_obj_type)
                else:
                    result = record_tableau_steps(tab, obj_type=detected_obj_type)
                self.send_json(result)
            except Exception as e:
                import traceback
                traceback.print_exc()
                self.send_json({"error": str(e)}, status_code=500)
            return

        elif path == "/api/demo":
            demo_name = query.get("name", ["beale"])[0]
            try:
                if demo_name == "beale":
                    lp = get_beale_cycling_example()
                elif demo_name == "twophase":
                    lp = get_twophase_example()
                elif demo_name == "standard2d":
                    lp = get_standard_2d_example()
                else:
                    self.send_json({"error": f"Unknown demo '{demo_name}'"}, status_code=400)
                    return
                result = record_lp_steps(lp)
                self.send_json(result)
            except Exception as e:
                import traceback
                traceback.print_exc()
                self.send_json({"error": str(e)}, status_code=500)
            return

        self.send_json({"error": "Endpoint not found."}, status_code=404)

    def do_POST(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        content_length = int(self.headers.get("Content-Length", 0))
        post_body = self.rfile.read(content_length)

        try:
            data = json.loads(post_body.decode("utf-8")) if post_body else {}

            if path == "/api/save-csv":
                filename = data.get("filename", "tableau.csv")
                content = data.get("content", "")
                filepath = resolve_csv_path(filename, for_writing=True)
                with open(filepath, "w", encoding="utf-8") as f:
                    f.write(content)
                self.send_json({"success": True, "message": f"Saved '{os.path.basename(filepath)}'"})
                return

            elif path == "/api/solve-tableau":
                headers = data.get("headers", [])
                rows = data.get("rows", [])
                force_tp = data.get("force_two_phase", True)
                if not headers or not rows:
                    self.send_json({"error": "Empty tableau submitted."}, status_code=400)
                    return

                row0_label = rows[0][0].strip() if rows[0] else "-z"
                obj_name = row0_label if row0_label else "-z"
                obj_type = "min" if obj_name.startswith("-") else "max"

                var_names = [h.strip() for h in headers[2:]]
                var_indices = {name: idx + 1 for idx, name in enumerate(var_names)}
                matrix = []
                for r in rows:
                    rhs_val = parse_fraction(r[1].strip())
                    coeffs = [parse_fraction(cell.strip()) for cell in r[2:]]
                    matrix.append(coeffs + [rhs_val])

                basic_vars = [r[0].strip() for r in rows[1:]]
                tab = SimplexTableau(var_names, var_indices, matrix, basic_vars, obj_name=obj_name)
                if force_tp:
                    result = record_tableau_steps_twophase(tab, obj_type=obj_type)
                else:
                    result = record_tableau_steps(tab, obj_type=obj_type)
                self.send_json(result)
                return

            elif path == "/api/solve-lp":
                obj_type = data.get("obj_type", "min")
                obj_coeffs = [parse_fraction(c) for c in data.get("obj_coeffs", [])]
                raw_constraints = data.get("constraints", [])
                force_two_phase = data.get("force_two_phase", False)
                constraints = []
                for c in raw_constraints:
                    coeffs = [parse_fraction(x) for x in c["coeffs"]]
                    rel = c["rel"]
                    rhs = parse_fraction(c["rhs"])
                    constraints.append((coeffs, rel, rhs))
                lp = LinearProgram(obj_type, obj_coeffs, constraints, force_two_phase=force_two_phase)
                result = record_lp_steps(lp)
                self.send_json(result)
                return

            self.send_json({"error": "Endpoint not found."}, status_code=404)
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.send_json({"error": str(e)}, status_code=500)

    def send_json(self, data, status_code=200):
        self.send_response(status_code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.end_headers()
        def fallback(obj):
            if isinstance(obj, Fraction):
                return format_fraction(obj)
            return str(obj)
        self.wfile.write(json.dumps(data, default=fallback).encode("utf-8"))


def find_free_port(start_port: int = 5050) -> int:
    port = start_port
    while port < 65535:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            if s.connect_ex(('127.0.0.1', port)) != 0:
                return port
        port += 1
    return start_port


def run_gui(port: int = 5050):
    port = find_free_port(port)
    server_address = ('127.0.0.1', port)
    httpd = HTTPServer(server_address, SimplexGUIRequestHandler)
    url = f"http://127.0.0.1:{port}"

    print("=" * 82)
    print("      GRAPHICAL USER INTERFACE (GUI) - SIMPLEX METHOD WITH BLAND'S RULE     ")
    print("=" * 82)
    print(f"Server successfully started at: {url}")
    print("Opening browser automatically...")
    print("Press Ctrl+C in terminal to stop server.")
    print("=" * 82 + "\n")

    try:
        webbrowser.open(url)
    except Exception as e:
        print(f"Could not open browser automatically: {e}")
        print(f"Please open this URL manually: {url}")

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down server...")
        httpd.server_close()
        print("Server stopped.")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1].endswith(".csv") and "--cli" in sys.argv:
        csv_file = sys.argv[1]
        resolved = resolve_csv_path(csv_file)
        display_path = resolved if os.path.exists(resolved) else csv_file
        print(f"Loading tableau from '{display_path}'...")
        try:
            tab, obj_type = load_tableau_from_csv(resolved)
            solve_imported_tableau(tab, obj_type)
        except Exception as e:
            print(f"[!] Error loading CSV file '{csv_file}': {e}")
    else:
        run_gui()
