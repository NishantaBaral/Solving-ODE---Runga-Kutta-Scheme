
import numpy as np
from typing import Callable
from PeriodComputation import rk23_find_period

def gauss_legendre(f: Callable[[float], float], a: float, b: float, n_nodes: int) -> float:
    """Gauss-Legendre quadrature on [a,b] with n_nodes.
    """
    nodes, weights = np.polynomial.legendre.leggauss(n_nodes)

    nodes = (b-a)/2*nodes+(a+b)/2
    sum = 0
    for i in range(n_nodes):
        sum = sum + weights[i]*f(nodes[i])  
    sum = sum*(b-a)/2
    return sum


f = lambda x: 1/np.sqrt(1+np.sin(x)**2)
a=0
b=np.pi/2
n_nodes = 20
integral = gauss_legendre(f, a, b, n_nodes)
C_quad = 4 * np.sqrt(2) * integral

a_values = [0.5, 1, 2, 4, 8] 

# Print the header (3 columns as requested, plus the 'a' column)
print(f"{'a':<5} | {'T_RK23':<15} | {'T_quad':<15} | {'Relative Error':<15}")
print("-" * 57)
for a in a_values:
    # Calculate the exact period using the scaling law
    T_quad = C_quad / a
    
    # Calculate the numerical period using your RK23 
    T_rk23, accepted, rejected = rk23_find_period(a, tol=1e-8)
    
    # Relative Error: |approx - exact| / exact
    rel_error = abs(T_rk23 - T_quad) / T_quad
    
    # Print the formatted row
    print(f"{a:<5.1f} | {T_rk23:<15.8f} | {T_quad:<15.8f} | {rel_error:<15.2e}")