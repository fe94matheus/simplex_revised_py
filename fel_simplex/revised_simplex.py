# Copyright (c) 2025 Felipe Oliveira
# This code is licensed under the MIT License.
# See the LICENSE file in the project root for more information.

import mpmath as mp


class RevisedSimplex:
    """
    A class implementing the Revised Simplex Method for solving linear programming problems.
    
    This implementation includes:
    - Handling of unrestricted variables through variable splitting (x = x⁺ - x⁻)
    - Automatic slack variable addition
    - High-precision arithmetic using mpmath

  
    """
    def __init__(self, A, b, c, precision=50):
        """
        Initialize the Revised Simplex Method solver.
        
        Args:
            A (mpmath.matrix): Constraint coefficients matrix
            b (mpmath.matrix): Right-hand side constraints
            c (mpmath.matrix): Objective function coefficients
            precision (int, optional): Number of decimal places for arithmetic precision.
                                     Defaults to 50.
        """
        self.original_A = A
        self.b = b
        self.original_c = c
        
        self.m = self.original_A.rows
        self.n = self.original_A.cols
        
        mp.mp.dps = int(precision);
        
        self.A, self.c = self.transform_unrestricted()
                
        # Add slack variables
        self.augmented_matrix, self.new_c = self.add_slack_variables()
        
        # Initialize basis
        self.basis = list(range(self.A.cols, self.A.cols + self.A.rows))
        self.nonbasis = list(range(self.A.cols))
        
        self.len_basis = len(self.basis)
        self.len_nonbasis= len(self.nonbasis)      
       
    def transform_unrestricted(self):
        """
        Transform unrestricted variables using the splitting method x = x⁺ - x⁻.
        
        This method doubles the number of variables, representing each original
        variable x as the difference of two non-negative variables (x⁺ and x⁻).
        
        Returns:
            tuple: (new_A, new_c) where:
                  - new_A is the transformed constraint matrix
                  - new_c is the transformed objective coefficients
        """
        new_A = mp.matrix(self.m, 2 * self.n)
        new_c = mp.matrix(2 * self.n, 1)
        
        for j in range(self.n):
            new_A[ : , j ] = self.original_A[ : , j ]
            new_c[ j , 0 ] = self.original_c[ j , 0 ]
            
            new_A[ : , j + self.n ] = -self.original_A[ : , j ]
            new_c[ j + self.n , 0 ] = -self.original_c[ j , 0 ]
            
        return new_A, new_c
        
    def add_slack_variables(self):
        """
        Add slack variables to create an initial basic feasible solution.
        
        Adds an identity matrix to the constraint matrix to create slack variables,
        which form the initial basis.
        
        Returns:
            tuple: (augmented_matrix, new_c) where:
                  - augmented_matrix includes the slack variables
                  - new_c is extended with zeros for slack variables
        """
        identity = mp.eye(self.m)
        augmented_matrix = mp.matrix(self.A.rows, self.A.cols + identity.cols)  
        
        augmented_matrix[:, :self.A.cols ] = self.A  
        augmented_matrix[:, self.A.cols :] = identity 
        
        new_c = mp.matrix(self.A.rows + self.A.cols, 1)
        
        new_c[:self.A.cols, 0] = self.c

        return augmented_matrix, new_c
        
        
    def get_basis_inverse(self):
        """
        Calculate the inverse of the current basis matrix.
        
        Extracts the columns corresponding to basic variables and computes
        the inverse using mpmath's matrix inversion.
        
        Returns:
            mpmath.matrix: Inverse of the current basis matrix
        """
        B = mp.matrix(self.len_basis, self.len_basis)  

        for i in range(self.len_basis):
            for j, basis_col in enumerate(self.basis):
                B[i, j] = self.augmented_matrix[i, basis_col]
        
        return B**-1
        
    def compute_reduced_costs(self, B_inv):
        """
        Compute reduced costs for all non-basic variables.
        
        The reduced cost for variable j is: cⱼ - yᵀaⱼ
        where y = (cᵦᵀB⁻¹)ᵀ and aⱼ is the column of variable j.
        
        Args:
            B_inv (mpmath.matrix): Inverse of the current basis matrix
            
        Returns:
            list: Pairs of (variable_index, reduced_cost) for non-basic variables
        """
        c_B = mp.matrix(self.len_basis, 1)
        
       
        for j, basis_col in enumerate(self.basis):
            c_B[j, 0] = self.new_c[basis_col, 0]
        
        y_T = c_B.T * B_inv
        
        reduced_costs = []
        
        for j, n_basis_col  in enumerate(self.nonbasis):
            r_cost = (self.new_c[n_basis_col, 0] - (y_T * self.augmented_matrix[ : , n_basis_col]))[0]
            reduced_costs.append((n_basis_col, r_cost))
        
        return reduced_costs
    
    def get_entering_variable(self, reduced_costs):
        """
        Select the entering variable using minimum reduced cost criterion.
        
        Args:
            reduced_costs (list): List of (variable_index, reduced_cost) pairs
            
        Returns:
            tuple: (entering_variable_index, minimum_reduced_cost)
        """
        min_reduced_cost = mp.inf
        entering_var = None
        
        for j, rc in reduced_costs:
            if rc < min_reduced_cost:
                min_reduced_cost = rc
                entering_var = j
                
        return entering_var, min_reduced_cost
    
    def get_leaving_variable(self, B_inv, entering_var):
        """
        Select the leaving variable using the minimum ratio test.
        
        Args:
            B_inv (mpmath.matrix): Inverse of the current basis matrix
            entering_var (int): Index of the entering variable
            
        Returns:
            tuple: (leaving_variable_index, maximum_ratio) or (None, None) if unbounded
        """
        a_j = mp.matrix(self.m, 1)
        a_j = self.augmented_matrix[:, entering_var]
       
        d = B_inv * a_j
        
        '''if all(d <= 0):
            return None, None  # Unbounded solution'''
        
        b_bar = B_inv * self.b
        
        ratios = []
        
        for i in range(self.len_basis):
            if d[i, 0] > 0:
                ratio =   b_bar[i, 0] / d[i, 0]
                ratios.append((self.basis[i], ratio))
                
        if not ratios:
            return None, None
                  
        leaving_var, max_ratio = min(ratios, key=lambda x: x[1])
            
        return leaving_var, max_ratio
        
        
    
    def update_basis(self, entering_var, leaving_var):
        """
        Update the basis and non-basis lists.
        
        Args:
            entering_var (int): Index of the entering variable
            leaving_var (int): Index of the leaving variable
        """
        self.basis[self.basis.index(leaving_var)] = entering_var
        self.nonbasis.remove(entering_var)
        self.nonbasis.append(leaving_var)
        
    def get_solution(self, B_inv):
        """
        Calculate the current basic solution.
        
        Args:
            B_inv (mpmath.matrix): Inverse of the current basis matrix
            
        Returns:
            mpmath.matrix: Current solution vector
        """
        b_star = B_inv * self.b
        x = mp.matrix(len(self.new_c), 1)
        for i, b in enumerate(self.basis):
            x[b] = b_star[i]
        return x
    
    def get_original_solution(self, transformed_solution):
        """
        Convert the solution from split variables back to original variables.
        
        Converts each pair of split variables (x⁺, x⁻) back to the original
        unrestricted variable x = x⁺ - x⁻.
        
        Args:
            transformed_solution (mpmath.matrix): Solution in terms of split variables
            
        Returns:
            tuple: (original_solution, original_objective_coefficients)
        """
        original_solution = mp.matrix(self.original_A.cols, 1)
        original_c = mp.matrix(self.original_A.cols, 1)
                
        # For each original variable, compute x = x⁺ - x⁻
        for i in range(self.original_A.cols):
            pos_part = transformed_solution[i]
            neg_part = transformed_solution[self.original_A.cols + i]
            original_solution[i] = pos_part - neg_part
        
        original_c[: , 0] = self.new_c[: self.original_A.cols , 0]
                
        return original_solution, original_c
    
    def solve(self, max_iterations=1000):
        """
        Solve the linear programming problem using the revised simplex method.
        
        Args:
            max_iterations (int, optional): Maximum number of iterations. Defaults to 1000.
            
        Returns:
            tuple: (x, obj_val, iterations, status) where:
                  - x is the optimal solution vector (or None if not found)
                  - obj_val is the optimal objective value (or None if not found)
                  - iterations is the number of iterations performed
                  - status is one of: "Optimal", "Unbounded", "Max iterations reached"
                  
        Notes:
            The algorithm terminates when:
            - An optimal solution is found (all reduced costs are non-negative)
            - The problem is determined to be unbounded
            - The maximum number of iterations is reached
        """
        iteration = 0
        
        while iteration < max_iterations:
            # Get basis inverse
            B_inv = self.get_basis_inverse()

            # Compute reduced costs
            reduced_costs = self.compute_reduced_costs(B_inv)
            
            # Get entering variable
            entering_var, min_reduced_cost = self.get_entering_variable(reduced_costs)
            
            # Check optimality  
            tolerance = mp.mpf(f'1e-{int(mp.mp.dps * 0.5)}')  
            
            if min_reduced_cost >= -tolerance:
                transformed_solution = self.get_solution(B_inv)
                x, c = self.get_original_solution(transformed_solution)
                
                return x, (c.T * x)[0], iteration, "Optimal"
            
            # Get leaving variable
            leaving_var, max_ratio = self.get_leaving_variable(B_inv, entering_var)

            # Check if unbounded
            if leaving_var is None:
                return None, None, "Unbounded"
            
            # Update basis
            self.update_basis(entering_var, leaving_var)
            
            iteration += 1
            
        return None, None, "Max iterations reached"



