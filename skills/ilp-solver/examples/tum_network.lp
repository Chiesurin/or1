# Totally Unimodular (TUM) Integer Linear Program
# Node-arc / interval structure ensures all square submatrix determinants in {-1, 0, 1}
# With integer RHS b, Hoffman-Kruskal guarantees an Ideal Polyhedron (P = P_I)
# Continuous relaxation directly yields an integer optimal solution without cutting planes

max x1 + x2
s.t.
x1 <= 3
x2 <= 4
x1 + x2 <= 5
x1, x2 >= 0
