#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 17:12:37 2025

@author: felipe
"""
import mpmath as mp
import fel_utils as futls
import revised_simplex as Simplex

def f(x):
    return mp.exp(str(x)) - mp.mpf(2)

setproblem = futls.SimplexProblemSetup(100)

points = setproblem.f_linspace(0, 3, 4)
print(points)

polynomial_degree = 3


'''A = setproblem.construct_matrix(polynomial_degree, points)
b = setproblem.construct_vector_b(f, points)
c = setproblem.construct_vector_c(polynomial_degree)'''

A = mp.matrix([[1, -1], [1, 1]])
b = mp.matrix([2, 6])
c = mp.matrix([-2, -1])


solver = Simplex.RevisedSimplex(A, b, c)

solution, objective_value, status = solver.solve()

print(f"Status: {status}")
print(f"Solution: {solution}")
print(f"Objective value: {objective_value}")

