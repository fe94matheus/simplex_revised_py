# Copyright (c) 2025 Felipe Matheus Oliveira Silva
# This code is licensed under the MIT License.
# See the LICENSE file in the project root for more information.

import fel_utils as futls
import revised_simplex as Simplex
import mpmath as mp
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,           
    "font.family": "serif",       
    "text.latex.preamble": r"\usepackage{amsmath}"   
})

class OptimalPolynomial:
    """
       A class for computing and visualizing optimal polynomial approximations of functions.
       
       This class uses linear programming techniques to find the best polynomial
       approximation of a given function under specified bounds. It includes
       functionality for:
       - Computing optimal polynomial coefficients
       - Visualizing the approximation alongside the original function
       - Handling arbitrary-precision arithmetic       
   """
    def  __init__(self, precision):
        """
        Initialize the OptimalPolynomial solver.
        
        Args:
            precision (int): Number of decimal places for arithmetic precision
        """
        self.prec = int(precision)
    
    def get_coefs(self, f, degree, a, b, num, omega_sup = lambda x: mp.mpf(1),  omega_inf = lambda x: mp.mpf(1)):
        """
        Compute optimal polynomial coefficients for approximating a function.
        
        This method finds the coefficients of a polynomial that best approximates
        the given function f(x) on the interval [a,b], subject to the constraints
        defined by omega_sup and omega_inf.
        
        Args:
            f (callable): Function to be approximated
            degree (int): Degree of the approximating polynomial
            a (float/str): Left endpoint of the approximation interval
            b (float/str): Right endpoint of the approximation interval
            num (int): Number of discretization points
            omega_sup (callable, optional): Upper weight function. Defaults to constant 1
            omega_inf (callable, optional): Lower weight function. Defaults to constant 1
            
        Returns:
            mpmath.matrix: Vector of polynomial coefficients in descending order
                         (highest degree first)
        
        Notes:
            The optimization problem is formulated as a linear programming problem
            and solved using the revised simplex method. The resulting polynomial
            minimizes the maximum absolute error between f(x) and the approximation.
        """ 
        setproblem = futls.SimplexProblemSetup(self.prec)
    
        self.points = setproblem.f_linspace(a, b, num)
        polynomial_degree = degree
    
        A = -1 * setproblem.construct_matrix(polynomial_degree, self.points, omega_sup, omega_inf)
        b = -1 * setproblem.construct_vector_b(f, self.points)
        c = setproblem.construct_vector_c(polynomial_degree)
    
        solver = Simplex.RevisedSimplex(A, b, c, setproblem.getPrec())
    
        self.solution, self.objective_value, self.it, self.status = solver.solve()
        
        x = self.solution[1 : ]
        self.x = x[::-1]
        
        return self.x
    
    def print_status(self):
        """
        Print the status and results of the polynomial approximation.
        
        Displays:
        - Solution status from the simplex solver
        - Optimal polynomial coefficients
        - Number of iterations taken
        - Achieved approximation error
        """
        
        print(f"Status: {self.status}")
        print(f"Coefficients: \n{self.x}")
        print(f"Iterations: {self.it}")
        print(f"Error: {self.objective_value}")
        
    def plot_fig(self, f, coefs):
        """
        Create a visualization of the original function and its polynomial approximation.
        
        Generates a publication-quality plot showing both the original function
        and its polynomial approximation using LaTeX rendering for labels and titles.
        
        Args:
            f (callable): Original function being approximated
            coefs (mpmath.matrix): Polynomial coefficients in descending order
            
        Notes:
            - The original function is shown in blue
            - The approximating polynomial is shown in dashed red
            
      
        """
        
        # Plotting
        f_vals = [f(mp.mpf(x)) for x in self.points]  # Evaluate the function f(x)
        p_vals = [mp.polyval(coefs, xx) for xx in self.points]  # Evaluate the polynomial

        # Convert results to floats for plotting
        f_vals_float = list(map(float, f_vals))
        p_vals_float = list(map(float, p_vals))
        
        # Create the plot
        plt.figure(figsize=(10, 6), dpi=300)
        plt.plot(self.points, f_vals_float, label="$f(x)$", color="blue", linewidth=2)
        plt.plot(self.points, p_vals_float, label="a(x)", color="red", linestyle="--", linewidth=2)
        plt.xlabel("$x$")
        plt.ylabel("$y$")
        plt.title("Function $f(x)$ and its Polynomial Approximation")
        plt.legend()
        plt.grid(True)
        plt.show()


