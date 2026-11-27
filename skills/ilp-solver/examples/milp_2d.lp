# 2D Mixed-Integer Linear Program (x1 integer, x2 continuous)
# Root continuous relaxation: x1 = 7/2 (3.5), x2 = 2, z = 37/2 (18.5)
# Optimal MILP solution: x1 = 3, x2 = 7/3 (~2.333), z = 55/3 (~18.333)

max 3*x1 + 4*x2
s.t.
2*x1 + 3*x2 <= 13
2*x1 + x2 <= 9
x1 >= 0, integer
x2 >= 0, continuous
