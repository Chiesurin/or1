#!/usr/bin/env python3
"""
Simplex Method with Bland's Rule
================================
A complete implementation of the Simplex Method (including Two-Phase Simplex
for general constraints) using Bland's Smallest-Subscript Rule to prevent cycling.

Features:
- Exact arithmetic using Python's `fractions.Fraction` (no float rounding errors).
- Step-by-step display of every simplex tableau with clean ASCII formatting.
- Detailed Bland's rule step-by-step explanations:
  * Entering variable: smallest index among variables with negative reduced costs.
  * Leaving variable: smallest basic variable index among rows tied for minimum ratio.
- Handles Maximization and Minimization problems.
- Handles all constraint types: <=, >=, and = via Two-Phase Simplex.
- Unboundedness and Infeasibility detection.
- Interactive manual input (guided step-by-step or quick line format),
  plus built-in classic benchmarks (e.g., Beale's cycling problem).
"""

import sys
import os
import csv
import re
import glob
from fractions import Fraction
from typing import List, Tuple, Dict, Optional


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


def get_display_obj_val(raw_val: Fraction, obj_type: str = "min", obj_name: str = "z") -> Fraction:
    """Returns the true objective value for display (ensuring non-negative cost when minimizing)."""
    if obj_name == "w" or obj_type == "min":
        return -raw_val if raw_val < 0 else raw_val
    return raw_val


class SimplexTableau:
    """
    Represents a Simplex Tableau.
    Variables are uniquely indexed with 1-based indexing for Bland's Rule:
    x1, x2, ..., xn, s1, s2, ..., e1, ..., y1, ...
    """
    def __init__(self, var_names: List[str], var_indices: Dict[str, int],
                 tableau: List[List[Fraction]], basic_vars: List[str], obj_name: str = "z"):
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
        display_cost = get_display_obj_val(self.tableau[0][-1], obj_type, self.obj_name)
        row0_cells = [self.obj_name, format_fraction(display_cost)]
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

        # Calculate optimal column widths
        num_columns = len(col_headers)
        col_widths = []
        for c in range(num_columns):
            max_len = max(len(table_cells[r][c]) for r in range(len(table_cells)))
            col_widths.append(max(max_len, 5))

        sep_line = "+" + "+".join("-" * (w + 2) for w in col_widths) + "+"

        print(sep_line)
        header_str = "| " + " | ".join(table_cells[0][c].center(col_widths[c]) for c in range(num_columns)) + " |"
        print(header_str)
        print(sep_line)

        for r in range(1, len(table_cells)):
            row_str = "| " + " | ".join(
                table_cells[r][c].center(col_widths[c]) if c == 0
                else table_cells[r][c].rjust(col_widths[c])
                for c in range(num_columns)
            ) + " |"
            print(row_str)
            if r == 1:
                # Separator after objective row
                print(sep_line)
        print(sep_line)

        # Print current basic solution
        print("\n  Current Basic Solution:")
        print(f"    {self.obj_name} = {format_approx(self.tableau[0][-1])}")
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


def record_simplex_steps(initial_tableau: SimplexTableau, phase_name: str = "",
                         forbidden_vars: Optional[set] = None,
                         obj_type: str = "min") -> Tuple[str, List[Dict], SimplexTableau]:
    """
    Executes the Simplex Method with Bland's Rule, recording every step into a structured list
    for GUI visualization.
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
        title = f"{phase_name} - Tableau {iteration}" if phase_name else f"Tableau {iteration}"
        if iteration == 0:
            title += " (Initial Tableau)"

        display_cost = get_display_obj_val(tab.tableau[0][-1], obj_type, tab.obj_name)
        step_data: Dict = {
            "iteration": iteration,
            "title": title,
            "phase_name": phase_name,
            "is_phase1": (phase_name == "Phase 1" or tab.obj_name == "w"),
            "is_phase2": (phase_name == "Phase 2"),
            "var_names": list(tab.var_names),
            "var_indices": dict(tab.var_indices),
            "basic_vars": list(tab.basic_vars),
            "obj_name": tab.obj_name,
            "obj_type": obj_type,
            "rows": [],
            "basic_solution": {},
            "obj_val": format_approx(display_cost),
            "raw_obj_val": format_fraction(display_cost),
            "status": "CONTINUE",
        }

        # Build rows for the UI table (RHS on left next to basic, NO z column)
        for r in range(tab.num_rows()):
            b_var = tab.obj_name if r == 0 else tab.basic_vars[r - 1]
            cells = [format_fraction(val) for val in tab.tableau[r][:-1]]
            if r == 0:
                rhs_val = display_cost
            else:
                rhs_val = tab.tableau[r][-1]
            rhs_str = format_fraction(rhs_val)
            step_data["rows"].append({
                "basic": b_var,
                "rhs": rhs_str,
                "cells": cells
            })
            if r > 0:
                step_data["basic_solution"][b_var] = rhs_str

        # Check candidate entering variables
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

        step_data["candidates"] = candidates

        if not candidates:
            step_data["status"] = "OPTIMAL"
            step_data["message"] = "All coefficients in Row 0 are >= 0. Optimality criterion satisfied!"
            steps.append(step_data)
            return "OPTIMAL", steps, tab

        # Bland's Rule: Choose candidate with smallest index
        candidates.sort(key=lambda c: c["index"])
        entering = candidates[0]
        step_data["entering_var"] = entering
        entering_col = entering["col"]

        # Ratio test
        ratios = []
        for row_idx in range(1, tab.num_rows()):
            b_var = tab.basic_vars[row_idx - 1]
            b_idx = tab.var_indices[b_var]
            a_ij = tab.tableau[row_idx][entering_col]
            rhs = tab.tableau[row_idx][-1]

            if a_ij > 0:
                ratio = rhs / a_ij
                ratios.append({
                    "row": row_idx,
                    "basic": b_var,
                    "index": b_idx,
                    "ratio": format_fraction(ratio),
                    "ratio_val": ratio,
                    "valid": True,
                    "rhs": format_fraction(rhs),
                    "coeff": format_fraction(a_ij)
                })
            else:
                reason = "denominator = 0" if a_ij == 0 else f"denominator < 0 ({format_fraction(a_ij)})"
                ratios.append({
                    "row": row_idx,
                    "basic": b_var,
                    "index": b_idx,
                    "valid": False,
                    "reason": reason,
                    "rhs": format_fraction(rhs),
                    "coeff": format_fraction(a_ij)
                })

        step_data["ratios"] = ratios
        valid_ratios = [r for r in ratios if r["valid"]]

        if not valid_ratios:
            step_data["status"] = "UNBOUNDED"
            step_data["message"] = f"All coefficients in column '{entering['name']}' are <= 0. The problem is UNBOUNDED!"
            steps.append(step_data)
            return "UNBOUNDED", steps, tab

        min_ratio_val = min(r["ratio_val"] for r in valid_ratios)
        tied_candidates = [r for r in valid_ratios if r["ratio_val"] == min_ratio_val]

        step_data["min_ratio"] = format_fraction(min_ratio_val)
        step_data["tied"] = len(tied_candidates) > 1

        if len(tied_candidates) > 1:
            tied_candidates.sort(key=lambda r: r["index"])
            leaving_item = tied_candidates[0]
            step_data["bland_tie_note"] = (
                f"Tie detected between {len(tied_candidates)} rows with ratio {leaving_item['ratio']}. "
                f"Bland's Rule selects the basic variable with the smallest index: '{leaving_item['basic']}' (Index {leaving_item['index']})."
            )
        else:
            leaving_item = tied_candidates[0]
            step_data["bland_tie_note"] = None

        pivot_row = leaving_item["row"]

        # Clean up Fraction ratio_val from dictionaries for JSON serialization
        for r in ratios:
            r.pop("ratio_val", None)
        leaving = dict(leaving_item)
        leaving.pop("ratio_val", None)
        step_data["leaving_var"] = leaving
        pivot_element = tab.tableau[pivot_row][entering_col]
        step_data["pivot"] = {
            "row": pivot_row,
            "col": entering_col,
            "val": format_fraction(pivot_element)
        }

        # Row operations
        ops = []
        ops.append(f"R{pivot_row} -> R{pivot_row} / ({format_fraction(pivot_element)})")
        for r in range(tab.num_rows()):
            if r != pivot_row:
                factor = tab.tableau[r][entering_col]
                if factor != 0:
                    r_name = tab.obj_name if r == 0 else f"R{r}"
                    sign = "-" if factor > 0 else "+"
                    ops.append(f"{r_name} -> {r_name} {sign} ({format_fraction(abs(factor))}) * R{pivot_row}")
        step_data["row_operations"] = ops

        steps.append(step_data)

        # Elementary row operations
        tab.tableau[pivot_row] = [x / pivot_element for x in tab.tableau[pivot_row]]
        for r in range(tab.num_rows()):
            if r != pivot_row:
                factor = tab.tableau[r][entering_col]
                if factor != 0:
                    tab.tableau[r] = [
                        tab.tableau[r][c] - factor * tab.tableau[pivot_row][c]
                        for c in range(len(tab.tableau[r]))
                    ]
        tab.basic_vars[pivot_row - 1] = entering["name"]

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

    slack_counter = 1
    surplus_counter = 1
    artificial_counter = 1

    constraint_aux: List[Dict[str, Optional[str]]] = []

    for coeffs, rel, rhs in norm_constraints:
        aux: Dict[str, Optional[str]] = {"slack": None, "surplus": None, "artificial": None}
        if rel == "<=":
            s_name = f"s{slack_counter}"
            slack_counter += 1
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
        phase1_tab = SimplexTableau(full_var_names, var_indices, phase1_tableau_data, basic_vars, obj_name="w")

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
    phase2_tab = SimplexTableau(phase2_var_names, phase2_var_indices, phase2_tableau_data,
                                phase2_basic_vars, obj_name="z")

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
        opt_z = get_display_obj_val(raw_z, lp.obj_type, "z")
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
            "solution_vector": sol_dict,
            "aux_vector": aux_dict,
            "basic_vars": list(final_tab.basic_vars),
            "alternative_optima": alt_optima
        }

    return result


class LinearProgram:
    """Represents a Linear Programming problem."""
    def __init__(self, obj_type: str, obj_coeffs: List[Fraction],
                 constraints: List[Tuple[List[Fraction], str, Fraction]],
                 var_names: Optional[List[str]] = None,
                 force_two_phase: bool = False):
        self.obj_type = obj_type.strip().lower()  # 'max' or 'min'
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
        print(f"  {self.obj_type.upper()} {obj_str}")
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

        slack_counter = 1
        surplus_counter = 1
        artificial_counter = 1

        constraint_aux: List[Dict[str, Optional[str]]] = []

        for coeffs, rel, rhs in norm_constraints:
            aux: Dict[str, Optional[str]] = {"slack": None, "surplus": None, "artificial": None}
            if rel == "<=":
                s_name = f"s{slack_counter}"
                slack_counter += 1
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
            phase1_tab = SimplexTableau(full_var_names, var_indices, phase1_tableau_data, basic_vars, obj_name="w")

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
        phase2_tab = SimplexTableau(phase2_var_names, phase2_var_indices, phase2_tableau_data,
                                    phase2_basic_vars, obj_name="z")

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
            opt_z = get_display_obj_val(raw_z, self.obj_type, "z")

            print(" [STATUS] OPTIMAL SOLUTION FOUND!")
            print(f" • Optimal Objective Value ({self.obj_type.upper()} Z*): {format_approx(opt_z)}")
            if self.obj_type == "min":
                print(f"   (Calculated from tableau value {format_fraction(raw_z)} -> Minimum Z* = {format_approx(opt_z)})")

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

    # 1. Objective Type
    while True:
        ot = input("\nProblem type (max / min) [default: max]: ").strip().lower()
        if ot in ("", "max", "maximize"):
            obj_type = "max"
            break
        elif ot in ("min", "minimize"):
            obj_type = "min"
            break
        print("Invalid input. Please enter 'max' or 'min'.")

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
    
    Maximize z = 10*x1 - 57*x2 - 9*x3 - 24*x4
    s.t.
      0.5*x1 - 5.5*x2 - 2.5*x3 + 9*x4 <= 0
      0.5*x1 - 1.5*x2 - 0.5*x3 + 1*x4 <= 0
      1*x1                            <= 1
      x1, x2, x3, x4 >= 0
    """
    obj_type = "max"
    obj_coeffs = [Fraction(10), Fraction(-57), Fraction(-9), Fraction(-24)]
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
    Standard textbook maximization problem:
    Maximize z = 3*x1 + 5*x2
    s.t.
      x1        <= 4
      2*x2      <= 12
      3*x1 + 2*x2 <= 18
      x1, x2 >= 0
    Optimal: x1 = 2, x2 = 6, z* = 36.
    """
    obj_type = "max"
    obj_coeffs = [Fraction(3), Fraction(5)]
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
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File '{file_path}' not found in current folder.")

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
    obj_name = "z"

    if has_basis_col:
        row0_label = data_rows[0][0].strip().lower()
        if row0_label and row0_label not in ("z", "obj", "cost"):
            obj_name = row0_label
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
                basic_vars.append(f"s{r}")

    var_indices = {name: idx + 1 for idx, name in enumerate(var_names)}
    tab = SimplexTableau(var_names, var_indices, matrix, basic_vars, obj_name=obj_name)
    # Default to minimization as requested
    return tab, "min"


def solve_imported_tableau(tab: SimplexTableau, obj_type: str = "min"):
    """Solves an imported tableau using Bland's rule and prints the final report."""
    print("\n" + "=" * 82)
    print(" SIMPLEX METHOD ON IMPORTED CSV TABLEAU")
    print("=" * 82)
    print(f" Target Direction: {obj_type.upper()}")
    print(" Variables and Bland's Rule Index Ordering:")
    for v_name in tab.var_names:
        print(f"   Index {tab.var_indices[v_name]:<2}: {v_name}")
    print(f" Initial Basis: {tab.basic_vars}")
    print("=" * 82)

    status, final_tab = solve_simplex(tab, obj_type=obj_type)

    print("\n" + "=" * 82)
    print(" FINAL RESULTS SUMMARY")
    print("=" * 82)

    if status == "OPTIMAL":
        final_tab.print_tableau(title="Final Optimal Tableau", obj_type=obj_type)
        raw_z = final_tab.tableau[0][-1]
        opt_z = get_display_obj_val(raw_z, obj_type, tab.obj_name)

        print(" [STATUS] OPTIMAL SOLUTION FOUND!")
        print(f" • Optimal Objective Value ({obj_type.upper()} {tab.obj_name}*): {format_approx(opt_z)}")
        if obj_type == "min":
            print(f"   (Calculated from tableau value {format_fraction(raw_z)} -> Minimum {tab.obj_name}* = {format_approx(opt_z)})")

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


# ---------------------------------------------------------------------
# MAIN MENU
# ---------------------------------------------------------------------

def main():
    # If a filename argument was provided on the command line, run it directly!
    if len(sys.argv) > 1 and sys.argv[1].endswith(".csv"):
        csv_file = sys.argv[1]
        print(f"Loading tableau from '{csv_file}'...")
        try:
            tab, obj_type = load_tableau_from_csv(csv_file)
            solve_imported_tableau(tab, obj_type)
        except Exception as e:
            print(f"[!] Error loading CSV file '{csv_file}': {e}")
        return

    # Look for CSV files in current folder
    csv_files = sorted(glob.glob("*.csv"))
    default_csv = "tableau.csv" if "tableau.csv" in csv_files else (csv_files[0] if csv_files else "")

    print("=" * 82)
    print("       SIMPLEX METHOD WITH BLAND'S RULE (PYTHON IMPLEMENTATION)       ")
    print("=" * 82)
    print("This program implements the Simplex Algorithm with Bland's Smallest-Subscript")
    print("Rule to guarantee anti-cycling and finite termination.")
    print("All arithmetic is computed exactly using rational fractions.\n")
    if csv_files:
        print(f"CSV files detected in current folder: {', '.join(csv_files)}")
    print("\nChoose an option:")
    default_opt_label = f" [default: {default_csv}]" if default_csv else ""
    print(f"  [1] Import starting tableau from CSV file in current folder{default_opt_label}")
    print("  [2] Enter a problem manually (Interactive step-by-step wizard)")
    print("  [3] Enter a problem manually (Paste multi-line text)")
    print("  [4] Run Demo: Beale's Cycling Problem (Proof of Bland's Rule anti-cycling)")
    print("  [5] Run Demo: Two-Phase Simplex (<=, >=, and = constraints)")
    print("  [6] Run Demo: Standard 2D Maximization (Textbook example)")
    print("  [g] Launch Graphical User Interface (Web Browser GUI)")
    print("  [q] Quit")
    print("-" * 82)

    choice = input("Enter choice (1-6, g, or q) [default: 1]: ").strip().lower()

    if choice in ("", "1"):
        if csv_files:
            print("\nAvailable CSV files in current folder:")
            for idx, fname in enumerate(csv_files, 1):
                marker = " (default)" if fname == default_csv else ""
                print(f"  [{idx}] {fname}{marker}")
            file_choice = input(f"\nEnter file number or filename [default: {default_csv}]: ").strip()
            if not file_choice:
                chosen_file = default_csv
            elif file_choice.isdigit() and 1 <= int(file_choice) <= len(csv_files):
                chosen_file = csv_files[int(file_choice) - 1]
            else:
                chosen_file = file_choice
        else:
            chosen_file = input("Enter path to CSV file [default: tableau.csv]: ").strip() or "tableau.csv"

        try:
            tab, obj_type = load_tableau_from_csv(chosen_file)
            solve_imported_tableau(tab, obj_type)
        except Exception as e:
            print(f"\n[!] Error loading CSV file '{chosen_file}': {e}")
            return
        return

    elif choice == "2":
        lp = prompt_guided_input()
    elif choice == "3":
        lp = prompt_multiline_input()
    elif choice == "4":
        print("\nLoading Beale's Cycling Problem (1955)...")
        print("Note: Dantzig's rule would cycle indefinitely with period 6 on this problem.")
        print("Bland's rule breaks ties using variable indices and terminates in 6 pivots!")
        lp = get_beale_cycling_example()
    elif choice == "5":
        print("\nLoading Two-Phase Simplex Problem (Min 4x1 + x2)...")
        lp = get_twophase_example()
    elif choice == "6":
        print("\nLoading Standard 2D Textbook Problem (Max 3x1 + 5x2)...")
        lp = get_standard_2d_example()
    elif choice in ("g", "gui"):
        from gui import run_gui
        run_gui()
        return
    elif choice in ("q", "quit", "exit"):
        print("Exiting.")
        return
    else:
        print("Invalid choice, defaulting to CSV import.")
        try:
            tab, obj_type = load_tableau_from_csv(default_csv or "tableau.csv")
            solve_imported_tableau(tab, obj_type)
        except Exception as e:
            print(f"\n[!] Error: {e}")
        return

    lp.solve()


if __name__ == "__main__":
    main()

