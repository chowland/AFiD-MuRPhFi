# Power spectra and time-averaging

## Spectrum calculation (`afid_spectra`)

The module [`afid_spectra`](https://github.com/chowland/AFiD-MuRPhFi/blob/main/src/spectra.F90) is used to compute and store power spectra from the simulation fields.
Since the computation of two-dimensional Fourier transforms involves the expensive operation of transposing the global 3-D arrays, spectra are not calculated by default.
This option is activated by including a file `spectra.in` in the simulation directory that specifies the time from which spectra should be calculated.

Fourier transforms are performed in the `y` and `z` directions in `CalcFourierCoef` to produce the discretized equivalent of the three-dimensional array

$$
\hat{f}(x,k_y, k_z) = \iint f(x,y,z) e^{-i (k_y y + k_z z)} \,\mathrm{d}y \,\mathrm{d}z
$$

Currently, the power spectra $\Phi_{ff}=|\hat{f}|^2$ of the three velocity components `vx`, `vy`, `vz` and temperature `temp` are computed and stored along with the co-spectrum $\Phi_{u\theta}=\hat{S}^* \hat{u}$ of the wall-normal concentration flux.
Each spectrum is computed at each time step and integrated in time to provide a final time-averaged three-dimensional spectrum that is written to `outputdir` by `WriteSpectra`.

## Time-averaged fields (`afid_averaging`)

A similar procedure is followed by the module `afid_averaging` to provide three-dimensional fields of time-averaged variables.
The same input file `spectra.in` is needed to activate this operation, and the averaging occurs from the time specified in the input file.
The subroutine `UpdateTemporalAverages` integrates the four basic variables in time at each time step to construct the temporal averages.
The time-averaged fields are eventually written to `outputdir/fields/xx_mean.h5`.