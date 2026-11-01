# Supported Tableau & Problem Input Formats

The Python solver engine supports multiple flexible input formats.

---

## 1. CSV File or Raw CSV Text (`--csv` or `--raw`)

### Format A: RHS on the Left (Recommended Standard)
```csv
Basis,RHS,x1,x2,x3,x4,x5
-z,0,-3,-5,0,0,0
x3,4,1,0,1,0,0
x4,12,0,2,0,1,0
x5,18,3,2,0,0,1
```

### Format B: RHS on the Right
```csv
Basis,x1,x2,x3,x4,x5,RHS
-z,-3,-5,0,0,0,0
x3,1,0,1,0,0,4
x4,0,2,0,1,0,12
x5,3,2,0,0,1,18
```

### Format C: With an explicit `z` Column (Stripped automatically)
```csv
Basis,z,x1,x2,x3,x4,RHS
-z,1,-3,-5,0,0,0
x3,0,1,0,1,0,4
x4,0,0,2,0,1,12
```

### Delimiters Supported:
- Comma `,`
- Semicolon `;`
- Tab `\t`
- Whitespace (space-separated matrix)

### Number Formats Supported:
- Integers: `0`, `4`, `-12`
- Irreducible / Reducible Fractions: `1/2`, `-11/2`, `3/5`
- Decimals: `0.5`, `-2.5`, `3.14`
- European Comma Decimals: `0,5`, `-2,5`

---

## 2. Inline Prompt Text (`--raw` or `--stdin`)

When a user pastes a tableau matrix directly in the chat, pass it via `--raw`:

```bash
python3 skills/tableau-solver/scripts/solve_tableau.py --raw "Basis,RHS,x1,x2,x3
-z,0,-10,57,9
x4,0,1/2,-11/2,-5/2
x5,1,1,0,0"
```

Or piped via stdin:
```bash
cat << 'EOF' | python3 skills/tableau-solver/scripts/solve_tableau.py --stdin
Basis,RHS,x1,x2,x3
-z,0,-10,57,9
x4,0,1/2,-11/2,-5/2
x5,1,1,0,0
EOF
```

---

## 3. Linear Program Formulation Text (`--lp`)

Pass standard algebraic linear programs directly via `--lp`:

```bash
python3 skills/tableau-solver/scripts/solve_tableau.py --lp "min -3*x1 - 5*x2
s.t.
  x1 <= 4
  2*x2 <= 12
  3*x1 + 2*x2 <= 18
  x1, x2 >= 0"
```

### Two-Phase LP Example:
```bash
python3 skills/tableau-solver/scripts/solve_tableau.py --lp "min 4*x1 + x2
s.t.
  3*x1 + x2 = 3
  4*x1 + 3*x2 >= 6
  x1 + 2*x2 <= 4
  x1, x2 >= 0"
```
The engine automatically detects the equality and greater-than-or-equal constraints and executes Two-Phase Simplex.
