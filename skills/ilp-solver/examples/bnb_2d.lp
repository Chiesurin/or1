# Classic 2D Integer Linear Program for Branch and Bound / Branch and Cut
# Root continuous relaxation: x1 = 9/4 (2.25), x2 = 15/4 (3.75), z = 165/4 (41.25)
# Optimal integer solution: x1 = 0, x2 = 5, z = 40

max 5*x1 + 8*x2
s.t.
x1 + x2 <= 6
5*x1 + 9*x2 <= 45
x1, x2 >= 0, integer
