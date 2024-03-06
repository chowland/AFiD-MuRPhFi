# Output file overview

```
outputdir/
-- cordin_info.h5
-- means.h5
-- (spectra.h5)
-- flowmov/
  -- movie_xcut.h5
  -- movie_ycut.h5
  -- movie_zcut.h5
-- fields/
  -- (00000_phi.h5)
  -- (00000_sal.h5)
  -- 00000_temp.h5
  -- 00000_vx.h5
  -- 00000_vy.h5
  -- 00000_vz.h5
  -- 00001...
  -- ...
```

## `cordin_info.h5`
This file contains 1-D arrays of the grid coordinates in the following format:
```
cordin_info.h5:
- /xc(r)
- /xm(r)
- /ym(r)
- /zm(r)
```
`xc` contains the wall-normal coordinates of the grid cell faces and is of length `nxm+1`.
`xm`, `ym`, and `zm` are the coordinates of the grid cell centres, and are arrays of length `nxm`, `nym`, and `nzm`.
The post-processing module `afidtools` reads this file by calling the class `afid.Grid(fld)`.

## `means.h5`
```
means.h5:
- /vybar
  -- /0000
  -- /0001
  -- ...
- /vxrms
- /epsilon
- /chi
```
Profiles as a function of $x$ and time, output at the regular interval defined by `t_out` in `bou.in`.
The helper function `read_mean()` is provided in the `afidtools` module to read this data.

## `flowmov`
```
flowmov/
  -- movie_xcut.h5:
    -- /temp
      -- /00000
      -- /00001
      -- ...
    -- (/sal)
      -- ...
    -- (/phi)
      -- ...
    -- /vx
      -- ...
    -- /vy
      -- ...
    -- /vz
      -- ...
  -- movie_ycut.h5:
    -- ...
  -- movie_zcut.h5:
    -- ...
```
The files in `flowmov` contain 2D slices of the state variables ($\boldsymbol{u}$, $T$, $S$, $\phi$).
The three plane slices produced are at $x\approx 0$ (taken at the first grid point from the wall), $y=L_y/2$ and $z=L_z/2$.
The slices are output at the regular time interval defined by `t_frame` in `bou.in`.

## `fields`

Each file in `fields` contains one 3D array from a single time snapshot, under the group name `/var`.