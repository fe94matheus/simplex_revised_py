# Copyright (c) 2025 Felipe Oliveira
# This code is licensed under the MIT License.
# See the LICENSE file in the project root for more information.

import mpmath as mp
import optimal_poly as op


def f(x):
    return mp.exp(str(x)) - mp.mpf(2)

def omega_sup(x):
    return 0

def omega_inf(x):
    return 0
        
p = op.OptimalPolynomial(15)

coefs = p.get_coefs(f, 1, 0, 3, 20)

p.plot_fig(f, coefs)



'''
Para polinômio de grau 15 e 25 pontos: Com precisão de 10 casas decimais, 
temos um problema de matriz numericamente singular.
'''

