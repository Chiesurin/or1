# Classical 2D Integer Linear Program for Gomory Cutting Planes
# Continuous relaxation optimum is fractional: x1 = 103/26 (~3.962), x2 = 21/13 (~1.615), z = 166/13 (~12.77)
# Optimal integer solution: x1 = 3, x2 = 2, z = 12

max 2*x1 + 3*x2
s.t.
2*x1 + 5*x2 <= 16
6*x1 + 2*x2 <= 27
x1, x2 >= 0
