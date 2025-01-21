# Copyright (c) 2025 Felipe Matheus Oliveira Silva
# This code is licensed under the MIT License.
# See the LICENSE file in the project root for more information.

import mpmath as mp

class SimplexProblemSetup:
    """
    A class to set up linear programming problems for function approximation.
    
    This class provides methods to construct the necessary matrices and vectors
    for solving approximation problems using the simplex method with arbitrary precision
    arithmetic through the mpmath library.
    
    Attributes:
        prec (int): The number of decimal places used for precision in calculations
    """
    
    def __init__(self, precision=50):
        """
        Initialize the SimplexProblemSetup with specified precision.
        
        Args:
            precision (int, optional): Number of decimal places for arithmetic precision.
                                     Defaults to 50.
        """
        mp.mp.dps = int(precision)
        self.prec = mp.mp.dps
        
    def f_linspace(self, start, end, num=50):
        """
        Create a linearly spaced sequence of points with high precision.
        
        Implements the formula: x_i = start + i * (end - start) / (num - 1)
        for i = 0, ..., num-1
        
        Args:
            start (float/str): The starting point of the sequence
            end (float/str): The ending point of the sequence
            num (int, optional): Number of points to generate. Defaults to 50.
            
        Returns:
            list: A list of mpf (arbitrary-precision floating-point) numbers
                 representing the linear space.
        """
        points = [mp.mpf(str(start)) + mp.mpf(i) * (mp.mpf(str(end)) - mp.mpf(str(start))) / 
                 (mp.mpf(num) - mp.mpf(1)) for i in range(num)]
        return points
    
    def construct_matrix(self, n, points, omega_sup, omega_inf):
        """
        Construct the constraint matrix for the linear programming problem.
        
        Creates a matrix of size (2m × (n+2)) where m is the number of points,
        incorporating both upper and lower bound constraints.
        
        Args:
            n (int): Degree of polynomial approximation
            points (list): List of points where constraints are evaluated
            omega_sup (callable): Upper weight function
            omega_inf (callable): Lower weight function
            
        Returns:
            matrix: An mpmath matrix representing the constraints
                   Size: (2m × (n+2))  
        """
        m = len(points)
        matrix = mp.matrix(2*m, n+2)
        
        # Construct upper bound constraints
        for i in range(m):
            matrix[i, 0] = omega_sup(points[i])
            matrix[i, 1] = mp.mpf(1)
            
        # Construct lower bound constraints
        for i in range(m, 2*m):
            matrix[i, 0] = omega_inf(points[i-m])
            matrix[i, 1] = mp.mpf(-1)
            
        # Fill polynomial terms
        for M in range(m):
            for N in range(1, n+1):
                matrix[M, N+1] = mp.power(points[M], N)
                matrix[M + m, N+1] = -mp.power(points[M], N)
            
        return matrix

    def construct_vector_b(self, f, points):
        """
        Construct the right-hand side vector for the linear programming problem.
        
        Args:
            f (callable): The function to be approximated
            points (list): List of points where the function is evaluated
            
        Returns:
            matrix: An mpmath column vector of size (2m × 1) where m = len(points)
                   containing function values at the specified points
        """
        m = len(points)
        b = mp.matrix(2 * m, 1)
        
        for i in range(m):
            b[i] = f(points[i])
            b[i + m] = -f(points[i])
        
        return b

    def construct_vector_c(self, n):
        """
        Construct the objective function coefficient vector.
        
        Creates a vector that represents the objective function coefficients
        for minimizing the approximation error.
        
        Args:
            n (int): Degree of polynomial approximation
            
        Returns:
            matrix: An mpmath column vector of size ((n+2) × 1) with 1 in the
                   first position and 0s elsewhere
        """
        c = mp.matrix(n + 2, 1)
        c[0] = mp.mpf(1)
        return c

    def getPrec(self):
        """
        Get the current precision setting.
        
        Returns:
            int: The number of decimal places used for precision
        """
        return self.prec