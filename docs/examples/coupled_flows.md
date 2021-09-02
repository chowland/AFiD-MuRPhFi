# Flows coupled to phase change

## 2-D Rayleigh-Bénard with a melting boundary
This example is taken from appendix A.3 of [Favier et al. (2019)](https://doi.org/10.1017/jfm.2018.773).
We consider the convection flow between a heated lower plate and a moving phase boundary.
There is no analytic form for the solution of this flow, so we can only perform a relative convergence study and compare to the results of [Favier et al. (2019)](https://doi.org/10.1017/jfm.2018.773).

The initial condition consists of a linear temperature profile subject to a perturbation in the liquid phase.
The melting temperature is set to $T_m=0.5$ and the intial height of the interface is at $x=0.5$.
Precisely, the initial fields are

$$
T(x,y) = \begin{cases}
            1 - x, & x\geq 0.5, \\
            1 - x + a \sin^2(2\pi x) \sin(4\pi y), & x < 0.5,
         \end{cases}
$$

$$
\phi(x,y) = \frac{1 + \tanh((x - 0.5)/2\varepsilon)}{2},
$$

and zero velocity.
Recall that $x$ is the wall-normal direction in our formulation.
The perturbation amplitude $a$ is not specified in [Favier et al. (2019)](https://doi.org/10.1017/jfm.2018.773), but we set it to 0.1.