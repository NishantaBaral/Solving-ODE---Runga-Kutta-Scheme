import numpy as np
import matplotlib.pyplot as plt

def rk23_solver(f, u0, t0, tf, tol):
    '''
    Runge-Kutta 2(3) method for solving ODEs. Parameters are:
    f: function representing the ODE system (du/dt = f(t, u))
    u0: initial condition
    t0: initial time, tf = final time
    tol = tolerance for error control
    '''
    t = t0
    u = np.array(u0, dtype=float)
    
    # Initial step size guess
    h = 0.1
    
    # Lists to store the accepted time steps and state vectors
    t_out = [t]
    u_out = [u]
    
    while t < tf:
        #Embedded RK23 Tableau 
        k1 = f(t, u)
        k2 = f(t + h/2, u + (h/2)*k1)
        k3 = f(t + 3*h/4, u + (3*h/4)*k2)
        
        # The 3rd-order update
        u_next = u + h * ( (2/9)*k1 + (1/3)*k2 + (4/9)*k3 )
        
        # Evaluate k4 at the new 3rd-order point 
        k4 = f(t + h, u_next)
        
        # The 2nd-order update
        u_next_hat = u + h * ( (7/24)*k1 + (1/4)*k2 + (1/3)*k3 + (1/8)*k4 )
        
        # Local Error Estimate maximum absolute difference between the 2nd and 3rd order estimates
        error = np.max(np.abs(u_next - u_next_hat))

        # Step Rejection and Adaptation
        if error <= tol:
            # Step Accepted!
            t = t + h
            u = u_next
            
            # Store the data
            t_out.append(t)
            u_out.append(u)
            
            # A larger step size for the next loop
            h = h * 0.9 * (tol / error)**(1/3)
            
        else:
            # Step Rejected! Try again with a smaller step size.
            h = h * 0.9 * (tol / error)**(1/3)

    return np.array(t_out), np.array(u_out)

def oscillator(t, u):
    # u[0] is y, u[1] is v
    return np.array([u[1], -u[0]**3])

# Pick an arbitrary amplitude (a = 2.0) and run it for 15 seconds
a = 2.0
t0 = 0.0
tf = 15.0
initial_conditions = [a, 0.0]

# Calling the solver
t_results, u_results = rk23_solver(oscillator, initial_conditions,t0,tf, tol=1e-8)

# Plot the position (y) over time
plt.plot(t_results, u_results[:, 0], 'b.-', label='Position (y)')
plt.title("RK23 Test: Nonlinear Oscillator (a=2.0)")
plt.xlabel("Time")
plt.legend()
plt.grid()
plt.show()