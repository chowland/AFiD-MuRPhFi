# Quantities saved to `means.h5`

The file `means.h5` contains a number of variables in the form $f(x,t)$, where $x$ is the wall-normal coordinate.
These variables are typically quantities averaged in the periodic ($y,z$) directions such that further averaging the profiles in $x$ then gives a time series for a volume average.
All profiles are evaluated at the cell centres.
Quantities involving the variables $S$ or $\phi$ are evaluated on the refined grid `xmr` whereas all other quantities are evaluated on the base grid `xm`.

## Mean profiles
The variables named `vybar`, `vzbar`, `Tbar`, `Sbar`, and `phibar` are simply $yz$-averages of the flow variables, i.e.

$$
\overline{f}(x,t) = \frac{1}{L_y L_z} \int_0^{L_z} \int_0^{L_y} f(x,y,z,t) \,\mathrm{d}y\,\mathrm{d}z ,
$$

for $f\in\{v,w,T,S,\phi\}$.
We will use the overbar as notation for the $yz$-average throughout this page.

*Note: due to the divergence-free velocity field, there is no point writing out $\overline{u}(x,t)$ since it must be identically zero.*

## Root mean square profiles
The variables `vxrms`, `vyrms`, `vzrms`, `Trms`, `Srms`, and `phirms` are defined as the square root of the $yz$-average of the squared flow variables, that is

$$
f_\mathrm{rms}(x,t) = \sqrt{\overline{f^2}} .
$$

Often it is useful to quantify the rms value for the turbulent fluctuations rather than the total.
We do not write this out from the simulation explicitly, since it can already be calculated from the output as

$$
f^\prime_\mathrm{rms}(x,t) = \sqrt{\overline{\left(f - \overline{f}\right)^2}} = \sqrt{ {f_\mathrm{rms}}^2 - \overline{f}^2} .
$$

## Fluxes

For the scalars $T$ and $S$, we also output the fluxes `vxT`, `vyT`, `vzT`, `vxS`, `vyS`, `vzS` defined as the components of $\overline{\boldsymbol{u}T}$ and $\overline{\boldsymbol{u}S}$.
Fluxes for salinity are computed using the interpolated velocity fields on the refined grid.
As above, it is often useful to separate the fluxes into a component due to the mean profiles and that due to turbulent fluctuations.
Thanks to the nice properties of the $yz$-average, this can be done easily as shown below for the turbulent flux of salinity in the $y$-direction:

$$
\overline{v'S'} = \overline{(v - \overline{v})(S - \overline{S})} = \overline{vS} - \overline{v}\overline{S} .
$$

Momentum fluxes (or Reynolds stresses) are also output to `means.h5`.
Since terms such as $\overline{uu}$ are already attainable from the rms quantities, we only need to save the cross-terms `vxvy`, `vxvz`, and `vyvz`, defined as $\overline{uv}$, $\overline{uw}$, and $\overline{vw}$ respectively.

Note that since $\overline{u}\equiv 0$, any flux of the form $\overline{uf}$ is equal to its turbulent counterpart $\overline{u'f'}$.

## Dissipation rates

Finally, we also output the dissipation rates of kinetic energy and scalar variance as `epsilon`, `chiT`, and `chiS`.
These are defined as

$$
\varepsilon(x,t) = \frac{1}{Re}\overline{\frac{\partial u_i}{\partial x_j} \frac{\partial u_i}{\partial x_j}} , \qquad \chi_T(x,t) = \frac{1}{RePr} \overline{|\nabla T|^2}, \qquad \chi_S(x,t) = \frac{1}{ReSc} \overline{|\nabla S|^2} ,
$$

where $Re$ is the bulk Reynolds number, equal to $\sqrt{Pr/Ra_T}$ if the temperature free-fall scale is used (`FFscaleS=0`) or equal to $\sqrt{Sc/Ra_S}$ if the salinity free-fall scale is used (`FFscaleS=1`).

As with the rms values and fluxes, it can be useful to determine the turbulent contributions to these dissipation rates.
This can be done by computing the wall-normal derivatives of the mean profiles (for which there is a helper function `ddx` in the `afidtools` module).
For example, the turbulent kinetic energy dissipation rate can be computed as

$$
\varepsilon'(x,t) = \varepsilon - \frac{1}{Re} \left|\frac{\partial \boldsymbol{\overline{u}}}{\partial x} \right|^2
$$

## Deriving other quantities from the output

The outputs can also be useful for determining key quantities such as energy budget terms and the rms displacement of the phase boundary.

### Energy budgets
As an example, we consider the energy budget of the volume-averaged turbulent kinetic energy (TKE) $E_K' = \frac{1}{2}\langle|\boldsymbol{u}'|^2\rangle$.
Here and on the rest of this page, we use $\langle \cdot\rangle$ to denote a volume average.
A volume average can be calculated from the mean profile data using the helper function `xmean` in the `afidtools` module.

Derived from the momentum equation, the evolution equation for the TKE reads

$$
\frac{dE_K'}{dt} = \mathcal{P}_S + \mathcal{J}' - \varepsilon' ,
$$

where the shear production $\mathcal{P}_S$ represents the transfer of energy between the mean flow and the turbulence, the buoyancy flux $\mathcal{J}'$ represents the exchange between kinetic and potential energy, and the dissipation rate $\varepsilon'$ has been discussed above.
The two new terms are defined as

$$
\mathcal{P}_S = -\left\langle \overline{u'\boldsymbol{u}'} \cdot \frac{\partial \overline{\boldsymbol{u}}}{\partial x} \right\rangle , \qquad
\mathcal{J}' = \langle w'b'\rangle = \langle w'(R_\rho T' - S') \rangle
$$

Example code for computing the shear production term can be found below
```python
folder = "/path/to/simulation/"
import afidtools as afid
grid = afid.Grid(folder)
vxvy, vxvz = afid.read_mean(folder, "vxvy"), afid.read_mean(folder, "vxvz")
vybar, vzbar = afid.read_mean(folder, "vybar"), afid.read_mean(folder, "vzbar")
dvydx, dvzdx = afid.ddx(vybar, grid.xm), afid.ddx(vzbar, grid.xm)
shear_prod = -vxvy*dvydx - vxvz*dvzdx
P_S = afid.xmean(shear_prod, grid.xc)
```

A function `energy_budgets` is also provided that collects all the relevant energy budget terms as time series into a Pandas DataFrame.