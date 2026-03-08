import numpy as np
import matplotlib.pyplot as plt

def F(t,X):
    x1 = X[0]
    x2 = X[1]

    dx1dt = -2*x1 - (x1**2)/(1+x2**4)
    dx2dt = -30*x2 - (x2**2)/(1+x1**4)

    return np.array([dx1dt, dx2dt])

#settimg up the time and initial conditions
h= 0.01 # step size
t = np.arange(0, 10+h, h) # time vector
X_results = np.zeros((len(t), 2))
X_results[0] = np.array([1.0, 1.0])

#Heun's method
for n in range(len(t) - 1):
    t_n = t[n]
    X_n = X_results[n]
    # Predictor step
    k1 = F(t_n, X_n)

#Guess the next value using the predictor step
    X_predictor = X_n + h * k1

    # Corrector: Find the slope at the predicted point (k2)
    t_next = t[n+1]
    k2 = F(t_next, X_predictor)
    # Final Step: Average the slopes to step forward
    X_results[n+1] = X_n + (h / 2) * (k1 + k2)

# 4. Plot the Results!
# X_results[:, 0] grabs all the x1 values, X_results[:, 1] grabs all the x2 values
plt.plot(t, X_results[:, 0], label='x1(t)', color='blue', linewidth=2)
plt.plot(t, X_results[:, 1], label='x2(t)', color='red', linestyle='--', linewidth=2)

plt.title("Heun's Method Solution for the Non-Linear System")
plt.xlabel("Time (t)")
plt.ylabel("State Variables")
plt.legend()
plt.grid(True)
plt.show()
    