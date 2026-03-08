We start with a non-linear oscilaltor $y''(t) + y(t)^3 = 0,$ with initial conditions $y(0) = a > 0,  y'(0) = 0.$
This equation corresponds to a particle moving in the potential $V(y) = \frac{1}{4}y^4.$ After rewriting the equation as a first-order system for
$$ u(t) = \begin{pmatrix} y(t) \\ v(t) \end{pmatrix}, \quad v = y',$$ we solve the ODE. RK23.py implements the Runge–Kutta method of order 2(3) (RK23) for
systems of the form $u'(t)=f(t,u).$ The implementation includes the embedded RK23 tableau, a local error estimate obtained from the difference of the order 2 and order 3 updates, step rejection and step-size adaptation, and a user-defined tolerance parameter tol.
As the first time $t_* > 0$ such that $y(t_*) = 0$ satisfies $T(a) = 4t_*,$ where $T(a)$ is the period, PeriodComputation.py detects the zero-crossing time inside the final acceptance step, 
i.e., computes $t_*$ by detecting the first sign change of $y(t)$.  For tolerance $\texttt{tol} = 10^{-8}$, and tabulates $T(a)$ for $a \in \{0.5,1,2,4,8\}$. 
Using conservation of energy, 
$$ T(a) = 4 \int_{0}^{a} \frac{dy}{\sqrt{2 \left(E - \frac{1}{4}y^4\right)}}, \quad E = \frac{1}{4}a^4. $$
So, GaussLegendre Comaprison,.py computes numerically
$$ C = 4\sqrt{2} \int_{0}^{1} \frac{du}{\sqrt{1 - u^4}} $$
using Gauss--Legendre quadrature (firstly by change of variables $\u = \sin\theta$. Outputs table comparing $T_{\text{RK23}}(a)$ and $T_{\text{quad}}(a)$ and their relative error (3 columns), for $a \in \{0.5, 1, 2, 4, 8\}$.

stiff syste,.py integrates the stiff system, using the Heun integration scheme,
    \[
    \dot x_1 = -2x_1-\dfrac{x_1^2}{1+x_2^4}, \quad \dot x_2 = -30x_2-\dfrac{x_2^2}{1+x_1^4}, \quad x(0)=(1,1).
    \]

