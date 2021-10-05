# Phase-field model validation

## 1-D melting from a hot boundary: phase-field parameters

We begin by investigating the model parameters $A$ and $C$ introduced by the phase-field method.
As a reminder, the governing equation for the phase field $\phi$ takes the form

$$
\partial_t \phi = \frac{A}{Pe_T} \nabla^2 \phi - \frac{A}{\varepsilon^2 Pe_T} \phi (1 - \phi) (1 - 2\phi + C(T - T_m)) ,
$$

so $A$ represents the diffusion coefficient of the phase field, and $C$ captures the strength of the coupling between the temperature field and the phase field.
As described in the [model description](equations.md#phase-field-method), the product $AC$ should be set by the Stefan number.

Since the thermal coupling term in the above equation is premultiplied by $1/\varepsilon^2$, where $\varepsilon$ is equal to the refined grid spacing, its magnitude can become very large.
This can cause numerical instability of the phase-field equation if the time step is too large.
From testing the [1-D melting example](examples/stefan.md) for a range of values of $A$ and $C$, we find the following restrictions for the time step:

$A$  | $C$ | ${\Delta t}_\mathrm{max}$
:--: | :-: | :------------------------
1    | 1   | $10^{-3}$
1    | 10  | $2\times 10^{-4}$
1    | 100 | $2\times 10^{-5}$
0.1  | 1   | $10^{-2}$
0.1  | 10  | $2\times 10^{-3}$
0.1  | 100 | $2\times 10^{-4}$
0.01 | 1   | $10^{-1}$
0.01 | 10  | $2\times 10^{-2}$
0.01 | 100 | $2\times 10^{-3}$

These tests were all performed with $Pe_T=10^3$, with a base resolution of $N_x = 512$ and a refined grid resolution of $N_x^r = 1024$.
It is clear from these results that the prefactor $AC$ constrains the maximum possible time step.
Since the model requires that $AC=1.2 \mathcal{S}^{-1}$, we anticipate the time step restriction to be less severe at high Stefan number.

With that in mind, let us now consider how the value of $C$ affects the results produced by the phase-field method.
As a reminder, $C=\varepsilon H\Delta T/\Gamma$, where $\Gamma$ is related to the surface energy of the interface.
Larger values of $C$ therefore reduce the impact of the Gibbs-Thomson effect, where high local curvature increases the local melting temperature.
In our 1-D example, the value of $C$ should therefore play no significant effect.

First, we consider the position of the phase boundary over time.
The analytic solution for the phase boundary height is

$$
h(t) = 2\Lambda \sqrt{(t+t_0)/Pe_T} ,
$$

where $\Lambda = 0.62$ for $\mathcal{S}=1$, and $t_0$ is a value such that $h(0)=0.1$.
In the below figure, we plot the analytic solution against the position of $\phi=0.5$ for the phase-field model with $1\leq C \leq 100$.


{% include "melting_validation_interface.html" %}

As expected, the value of $C$ appears to have little impact on the phase-field solution.
Small differences between the simulations develop, but this is likely due to the initial condition of the phase-field, which is taken to be the equilibrium solution $\phi=0.5(1 - \tanh ((x-0.1)/2\varepsilon))$.
This initial condition may not perfectly represent the initial motion of the interface at early times.

Now we can take a closer look at the temperature profiles as we vary the parameters.
Below we plot the temperature profile $T(x,t)$ for a range of discrete times between $t=0$ and $t=100$.

{% include "melting_validation_temperature.html" %}

The results at first seem almost indistinguishable from the analytic solution for each value of $C$.
However on close inspection, we can see that the temperature in the solid phase near the phase boundary is slightly above the melting temperature $T=0$ for the cases $C>1$.
This does not appear to play a significant role in the development of the temperature field in the liquid, or the interface position, but this result should be kept in mind when proceeding with the other validation cases.

## 1-D Freezing from a cold boundary: resolution study

Now we have a good idea of how the model parameters behave, we perform a convergence study to investigate how the accuracy of the simulation depends on the resolution, both of the base grid and the refined grid.
We consider the inverse problem to that above, where a liquid phase freezes from a cold boundary and a solidification front moves across the domain.
Thanks to the symmetry with the previous problem, the position of the phase boundary should be unchanged from that above.
Following the above results, we perform these tests with phase-field parameters fixed at $A=0.01$, $C=100$.
We vary the size of the base grid `nxm` between 64 and 512, and vary the size of the refined grid `nxmr` between 64 and 1024.

{% include "freezing_convergence_interface.html" %}

The results highlight the importance of high resolution for the phase-field variable.
Even for the lowest base resolution of $n_x=64$, we find that the phase-field model tracks the interface reasonably accurately if the refined grid is sufficiently large, say $n_x^r=512$.
Such a disparity between the resolutions does however result in some oscillations in the interface position due to the coupling between $\phi$ and the temperature field.
Such oscillations only appear significant when the refined grid is more than twice as fine as the base grid.

Inspection of the temperature profiles shows a similar trend, where under-resolution of the phase-field can lead to excessive diffusion at the phase boundary, and significant deviations from the melting temperature.

{% include "freezing_convergence_temperature.html" %}

## 1-D Freezing of a supercooled melt

As we shall see, the stability of the phase-field model is more sensitive to the parameter choices for the case of supercooling.
In this example, the solid phase grows out of a liquid that is cooled to below the melting temperature.
An analytic solution can be derived for a semi-infinite domain, which we can compare to our finite domain when the interface is sufficiently far from the upper boundary.

For $C=1$ the one-dimensional growth of the solid phase matches the analytic solution well, even as the effect of the upper boundary becomes more significant.
As seen below, both the temperature profiles and the interface position are consistent with the analytic solution.

{% include "supercooling_profile_validation.html" %}

{% include "supercooling_interface_validation.html" %}

However, for other values of $C$ the model breaks down.
The phase field $\phi$ takes begins to take values away from $0$ and $1$, and a solid phase begins to grow from the upper boundary.
The heat fluxes also appear inconsistent at the phase boundary, with relatively large temperatures spreading throughout the domain.

{% include "supercooling_errors.html" %}

The reason for the instability of the system is currently unclear, but given these results, it seems sensible to fix the model parameter $C=1$ in simulations where the effect of supercooling may occur.

## 2-D axisymmetric melting

<video width="100%" controls>
  <source src="../assets/AxisymMelting.mp4" type="video/mp4">
</video>

## 2-D Rayleigh-Bénard convection with a melting boundary

Now, we introduce an example where the evolution of the temperature field and the phase boundary are also coupled to the fluid flow.
This again follows a validation example used by [Favier et al. (2019)](), and details of the initial condition can be found on the [examples page](examples/coupled_flows.md).
The system consists of a fluid heated from below that underlies a solid phase.
The evolution of the system is shown in the following video:

<video width="100%" controls>
  <source src="../assets/MeltingRBC.mp4" type="video/mp4">
</video>

We follow Favier et al. in performing a convergence study, taking into account the resolutions of the two grids separately.
The following plot presents the height of the phase boundary $\phi=0.5$ at the time $t=0.1 Pe$ for a range of resolutions.

{% include "melting_RBC_interface_convergence.html" %}

As in the previous cases, we see that the resolution of the phase field is most important for the accurate evolution of the phase boundary.
This highlights a key advantage of our multiple-resolution technique, where only the phase-field variable needs evolving on the largest grid.
The above figure presents results for the model parameters $A=0.01$, $C=100$.

As we described above, the product $AC$ should be set by the physical parameters.
This leaves the parameter $C$ as the sole free model parameter.
Below, we vary the value of $C$ while keeping $AC=1$ fixed, with a resolution of $n_x=512$, $n_x^r=1024$.
We observe monotonic dependence of the interface position with $A$, but only small variations relative to those seen in the resolution study.

{% include "melting_RBC_interface_parameters.html" %}
