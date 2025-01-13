#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 16:11:04 2025

@author: felipe
"""

import mpmath as mp


class SimplexProblemSetup:
    
    def __init__(self, precision = 50):
        mp.mp.dps = int(precision);
        self.prec =  mp.mp.dps
        
    def f_linspace(self, start, end, num=50):
        '''
            x_i = start + i * (end - start) / (num - 1)
            i = 0, ... , num -1 
        '''
        points = [mp.mpf(str(start)) + mp.mpf(i) * (mp.mpf(str(end)) - mp.mpf(str(start))) / (mp.mpf(num) - mp.mpf(1)) for i in range(num)]
        
        return points
    
    def construct_matrix(self, n, points, omega = lambda x: mp.mpf(1)):
        m = len(points)
        matrix = mp.matrix(2*m , n+2) 
        
        for i in range(m):
            matrix[i, 0] = omega(points[i])
            matrix[i, 1] = mp.mpf(1)
            
        for i in range(m, 2*m):
            matrix[i, 0] = omega(points[i-m])
            matrix[i, 1] = mp.mpf(-1)
            
        for M in range(m):
            for N in range(1, n+1):
                matrix[M, N+1] = mp.power(points[M], N)
                matrix[M + m, N+1] = -mp.power(points[M], N)
            
        return matrix

    def construct_vector_b(self, f, points):
        m = len(points)
        b = mp.matrix(2 * m, 1)
        
        for i in range(m):
            b[i] = f(points[i])
            b[i + m] = -f(points[i])
        
        return b

    def construct_vector_c(self, n):
        c = mp.matrix(n + 2, 1)
        c[0] = mp.mpf(1)
        
        return c


setproblem = SimplexProblemSetup(100)



