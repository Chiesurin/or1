#!/usr/bin/env python3
"""
Simplex Method with Bland's Rule - Native Desktop GUI (Tkinter)
==============================================================
A desktop window application built with Python's standard `tkinter`.

Usage:
  python3 gui_tk.py
"""

import os
import sys
import glob
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from fractions import Fraction

from simplex import (
    SimplexTableau,
    load_tableau_from_csv,
    record_tableau_steps,
    format_fraction,
    format_approx,
)


class SimplexTkinterApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("Simplex Method with Bland's Rule")
        self.root.geometry("980x720")
        self.root.minsize(800, 600)

        self.steps_data = []
        self.current_step_idx = 0

        self.setup_ui()
        self.load_default_csv()

    def setup_ui(self):
        # Top toolbar
        toolbar = ttk.Frame(self.root, padding="10")
        toolbar.pack(fill=tk.X)

        ttk.Label(toolbar, text="CSV File:", font=("Helvetica", 11, "bold")).pack(side=tk.LEFT, padx=(0, 6))

        self.csv_combo = ttk.Combobox(toolbar, state="readonly", width=25)
        self.refresh_csv_list()
        self.csv_combo.pack(side=tk.LEFT, padx=4)
        self.csv_combo.bind("<<ComboboxSelected>>", self.on_csv_selected)

        ttk.Button(toolbar, text="Browse...", command=self.browse_csv).pack(side=tk.LEFT, padx=4)
        ttk.Button(toolbar, text="Solve CSV", command=self.solve_selected_csv).pack(side=tk.LEFT, padx=8)

        # Status and title
        title_frame = ttk.Frame(self.root, padding="10 0 10 10")
        title_frame.pack(fill=tk.X)

        self.title_label = ttk.Label(title_frame, text="Initial Tableau", font=("Helvetica", 13, "bold"))
        self.title_label.pack(side=tk.LEFT)

        self.obj_label = ttk.Label(title_frame, text="Z = 0", font=("Helvetica", 12, "bold"), foreground="#7c3aed")
        self.obj_label.pack(side=tk.RIGHT)

        # Treeview (Tableau Table)
        table_frame = ttk.Frame(self.root, padding="10 0 10 0")
        table_frame.pack(fill=tk.BOTH, expand=True)

        self.tree = ttk.Treeview(table_frame, show="headings", height=8)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        scrollbar_y = ttk.Scrollbar(table_frame, orient=tk.VERTICAL, command=self.tree.yview)
        scrollbar_y.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree.configure(yscrollcommand=scrollbar_y.set)

        # Stepper Controls
        control_frame = ttk.Frame(self.root, padding="10")
        control_frame.pack(fill=tk.X)

        self.btn_first = ttk.Button(control_frame, text="⏮ First", command=self.first_step)
        self.btn_first.pack(side=tk.LEFT, padx=2)

        self.btn_prev = ttk.Button(control_frame, text="◀ Prev", command=self.prev_step)
        self.btn_prev.pack(side=tk.LEFT, padx=2)

        self.step_label = ttk.Label(control_frame, text="Step 0 / 0", font=("Helvetica", 11, "bold"))
        self.step_label.pack(side=tk.LEFT, padx=12)

        self.btn_next = ttk.Button(control_frame, text="Next ▶", command=self.next_step)
        self.btn_next.pack(side=tk.LEFT, padx=2)

        self.btn_last = ttk.Button(control_frame, text="Last ⏭", command=self.last_step)
        self.btn_last.pack(side=tk.LEFT, padx=2)

        self.status_label = ttk.Label(control_frame, text="In progress", font=("Helvetica", 10, "italic"))
        self.status_label.pack(side=tk.RIGHT, padx=10)

        # Bottom Notebook: Analysis details & Solution
        notebook = ttk.Notebook(self.root)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Tab 1: Bland's Rule Analysis
        self.text_analysis = tk.Text(notebook, wrap=tk.WORD, font=("Courier", 11), height=8)
        notebook.add(self.text_analysis, text="🔍 Step Analysis (Bland's Rule)")

        # Tab 2: Solution Summary
        self.text_solution = tk.Text(notebook, wrap=tk.WORD, font=("Courier", 11), height=8)
        notebook.add(self.text_solution, text="🎉 Optimal Solution")

    def refresh_csv_list(self):
        files = sorted(glob.glob("*.csv"))
        self.csv_combo["values"] = files
        if "tableau.csv" in files:
            self.csv_combo.set("tableau.csv")
        elif files:
            self.csv_combo.set(files[0])

    def browse_csv(self):
        filename = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if filename:
            self.solve_csv_file(filename)

    def on_csv_selected(self, event):
        self.solve_selected_csv()

    def load_default_csv(self):
        default = self.csv_combo.get()
        if default and os.path.exists(default):
            self.solve_csv_file(default)

    def solve_selected_csv(self):
        selected = self.csv_combo.get()
        if selected:
            self.solve_csv_file(selected)

    def solve_csv_file(self, filepath: str):
        try:
            tab, obj_type = load_tableau_from_csv(filepath)
            result = record_tableau_steps(tab, obj_type)
            self.steps_data = result["steps"]
            self.current_step_idx = 0
            self.render_step()
        except Exception as e:
            messagebox.showerror("Error", f"Could not solve CSV file:\n{e}")

    def render_step(self):
        if not self.steps_data:
            return

        step = self.steps_data[self.current_step_idx]
        total = len(self.steps_data)

        self.title_label.config(text=step["title"])
        self.obj_label.config(text=f"{step['obj_name'].upper()} = {step['obj_val']}")
        self.step_label.config(text=f"Step {self.current_step_idx + 1} of {total}")

        if step["status"] == "OPTIMAL":
            self.status_label.config(text="OPTIMAL REACHED", foreground="#16a34a")
        elif step["status"] == "UNBOUNDED":
            self.status_label.config(text="UNBOUNDED", foreground="#dc2626")
        else:
            self.status_label.config(text=f"Iteration {step['iteration']}", foreground="#2563eb")

        # Setup Treeview Columns
        cols = ["Basis", step["obj_name"]] + step["var_names"] + ["RHS"]
        self.tree["columns"] = cols
        for c in cols:
            self.tree.heading(c, text=c)
            self.tree.column(c, width=70, anchor=tk.CENTER)

        # Clear existing rows
        for item in self.tree.get_children():
            self.tree.delete(item)

        # Populate rows
        for r_idx, r in enumerate(step["rows"]):
            vals = [r["basic"], r["z_col"]] + r["cells"] + [r["rhs"]]
            self.tree.insert("", tk.END, values=vals)

        # Populate Analysis Text
        self.text_analysis.delete("1.0", tk.END)
        lines = []

        if step["status"] == "OPTIMAL":
            lines.append("=" * 60)
            lines.append("OPTIMALITY CRITERION SATISFIED!")
            lines.append("All Row 0 coefficients are >= 0.")
            lines.append("=" * 60)
        else:
            lines.append("--- 1. ENTERING VARIABLE SELECTION (BLAND'S RULE) ---")
            if step.get("candidates"):
                for c in step["candidates"]:
                    chosen = " [CHOSEN BY BLAND]" if c["name"] == step["entering_var"]["name"] else ""
                    lines.append(f"  • {c['name']:<5} (Index {c['index']:<2}) | Reduced cost: {c['coeff']}{chosen}")
                lines.append(f"Bland's Rule: selected '{step['entering_var']['name']}' (minimum index among candidates)\n")

            lines.append("--- 2. MINIMUM RATIO TEST (LEAVING VARIABLE) ---")
            if step.get("ratios"):
                for r in step["ratios"]:
                    if r["valid"]:
                        leaving = " [LEAVES BASIS]" if r["row"] == step["leaving_var"]["row"] else ""
                        lines.append(f"  • Row {r['row']} ({r['basic']:<5}, Index {r['index']:<2}): {r['rhs']} / {r['coeff']} = {r['ratio']}{leaving}")
                    else:
                        lines.append(f"  • Row {r['row']} ({r['basic']:<5}): No ratio ({r['reason']})")
                if step.get("bland_tie_note"):
                    lines.append(f"\n[TIE DETECTED]: {step['bland_tie_note']}")
                lines.append(f"Leaving variable: '{step['leaving_var']['basic']}' (Row {step['leaving_var']['row']})\n")

            lines.append("--- 3. PIVOT & ELEMENTARY ROW OPERATIONS ---")
            if step.get("pivot"):
                lines.append(f"Pivot on Row {step['pivot']['row']}, Column '{step['entering_var']['name']}' (Value: {step['pivot']['val']})")
                for op in step.get("row_operations", []):
                    lines.append(f"  • {op}")

        self.text_analysis.insert(tk.END, "\n".join(lines))

        # Populate Solution Tab
        self.text_solution.delete("1.0", tk.END)
        last_step = self.steps_data[-1]
        s_lines = []
        if last_step["status"] == "OPTIMAL":
            s_lines.append(f"OPTIMAL OBJECTIVE VALUE: {last_step['obj_name'].upper()}* = {last_step['obj_val']}\n")
            s_lines.append("OPTIMAL BASIC SOLUTION:")
            for k, v in last_step["basic_solution"].items():
                s_lines.append(f"  • {k:<6} = {v} (Basic)")
            non_basics = [v for v in last_step["var_names"] if v not in last_step["basic_solution"]]
            if non_basics:
                s_lines.append(f"\nNon-basic variables (= 0): {', '.join(non_basics)}")
        elif last_step["status"] == "UNBOUNDED":
            s_lines.append("THE LINEAR PROGRAM IS UNBOUNDED.")
        self.text_solution.insert(tk.END, "\n".join(s_lines))

    def first_step(self):
        self.current_step_idx = 0
        self.render_step()

    def prev_step(self):
        if self.current_step_idx > 0:
            self.current_step_idx -= 1
            self.render_step()

    def next_step(self):
        if self.current_step_idx < len(self.steps_data) - 1:
            self.current_step_idx += 1
            self.render_step()

    def last_step(self):
        self.current_step_idx = len(self.steps_data) - 1
        self.render_step()


def run_tkinter():
    root = tk.Tk()
    app = SimplexTkinterApp(root)
    root.mainloop()


if __name__ == "__main__":
    run_tkinter()
