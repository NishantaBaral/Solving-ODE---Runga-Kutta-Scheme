import numpy as np

def oscillator(t, u):
    return np.array([u[1], -u[0]**3])

def rk23_find_period(a, tol=1e-8):
    f = oscillator
    # Initial setup
    t = 0.0
    u = np.array([a, 0.0]) 
    h = 0.1
    
    # Counters for the final table
    accepted_steps = 0
    rejected_steps = 0
    
    
    while True:
        #RK23 Tableau.
        k1 = f(t,u)
        k2 = f(t + h/2, u + (h/2)*k1)
        k3 = f(t + 3*h/4, u + (3*h/4)*k2)
        u_next = u + h * ( (2/9)*k1 + (1/3)*k2 + (4/9)*k3 )
        k4 = f(t + h, u_next)
        u_next_hat = u + h * ( (7/24)*k1 + (1/4)*k2 + (1/3)*k3 + (1/8)*k4 )
        
        # Error Estimate
        error = np.max(np.abs(u_next - u_next_hat))
        
        # Step Acceptance and Event Detection 
        if error <= tol:
            accepted_steps += 1
            
            y_old = u[0]
            y_new = u_next[0]
            
            # EVENT DETECTION
            if y_old > 0 and y_new <= 0:
                # Calculate exact crossing time using linear interpolation
                t_star = t + h * ( (0 - y_old) / (y_new - y_old) )
                
                # The period is 4 times this crossing time
                period = 4 * t_star

                return period, accepted_steps, rejected_steps
            
            # If we didn't cross zero, we update normally
            t = t + h
            u = u_next
            
            # Increase h
            h = h * 0.9 * (tol / error)**(1/3)
            
        else:
            # Step rejected
            rejected_steps += 1
            
            # Decrease h
            h = h * 0.9 * (tol / error)**(1/3)

# The list of amplitudes
a_values = [0.5, 1.0, 2.0, 4.0, 8.0]
tolerance = 1e-8

print(f"{'a':<5} | {'T(a)':<15} | {'Accepted Steps':<15} | {'Rejected Steps':<15}")
print("-" * 55)

for a in a_values:
    T_a, acc, rej = rk23_find_period(a, tol=tolerance)
    
    # Print a nicely formatted row
    print(f"{a:<5.1f} | {T_a:<15.6f} | {acc:<15} | {rej:<15}")