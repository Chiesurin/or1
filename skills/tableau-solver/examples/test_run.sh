#!/usr/bin/env bash
# Test all execution modes for the Tableau Solver Skill

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOLVER="${SCRIPT_DIR}/../scripts/solve_tableau.py"

echo "=== 1. Test Solve from CSV File ==="
python3 "$SOLVER" --csv sample_tableau.csv

echo -e "\n=== 2. Test Solve from Raw Inline String ==="
python3 "$SOLVER" --raw "Basis,RHS,x1,x2,x3,x4
-z,0,-3,-5,0,0
x3,4,1,0,1,0
x4,12,0,2,0,1"

echo -e "\n=== 3. Test Two-Phase LP Formulation with Markdown Report Export ==="
python3 "$SOLVER" --lp "min 4*x1 + x2
s.t.
  3*x1 + x2 = 3
  4*x1 + 3*x2 >= 6
  x1 + 2*x2 <= 4
  x1, x2 >= 0" --output-md "${SCRIPT_DIR}/test_twophase_solution.md"

echo -e "\n=== 4. Test JSON Output ==="
python3 "$SOLVER" --csv sample_tableau.csv --format json

echo -e "\nAll tests completed successfully!"
