# Copyright (c) 2025 Felipe Matheus Oliveira Silva
# This code is licensed under the MIT License.
# See the LICENSE file in the project root for more information.

import mpmath as mp
import fel_utils as futls
import revised_simplex as Simplex

def f(x):
    return mp.exp(str(x)) - mp.mpf(2)

setproblem = futls.SimplexProblemSetup(100)

points = setproblem.f_linspace(0, 3, 4)

polynomial_degree = 3


'''A = setproblem.construct_matrix(polynomial_degree, points)
b = setproblem.construct_vector_b(f, points)
c = setproblem.construct_vector_c(polynomial_degree)'''

'''
min -3x1 + 2x2 - x3
s.t. 2x1 -  x2 + x3 <= 10
     -x1 - 2x2 + x3 <= 4
'''

A = mp.matrix( [ [2, -1, 1], [-1, -2, 1] ] )
b = mp.matrix( [10, 4] )
c = mp.matrix( [-3, 2, -1] )


solver = Simplex.RevisedSimplex(A, b, c)

solution, objective_value, status = solver.solve()
'''
print(f"Status: {status}")
print(f"Solution: \n{solution}")
print(f"Objective value: {objective_value}")'''

