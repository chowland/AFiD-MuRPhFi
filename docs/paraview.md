# 3D visualization with ParaView

## Connecting remotely to a ParaView HPC session

Many HPC systems are set up with ParaView installed, and allow you to make use of their superior computational power to load and render 3-D fields.
This can simply be done by running ParaView on the HPC and using a VNC to access the GUI on the remote system.
But ParaView also allows you to connect a local ParaView installation to a server running on the HPC.
This way, the remote system only focuses on rendering the data and does not need to generate the GUI and stream it via VNC.

The Snellius [visualisation documentation](https://servicedesk.surf.nl/wiki/display/WIKI/ParaView+client-sever+mode+for+parallel+visualization) by SURFsara provides a walkthrough for setting up this ParaView connection.
Optionally, you can request access to the GPU visualisation nodes, which will give you better rendering performance.

We will show how to connect to a CPU node here, but for sizeable visualisation tasks, we recommend using GPU nodes.
The script [`paraview_submit.sh`]() in the `tools` directory can be used to launch a ParaView server on a CPU node, which calls 
```
module load 2022
module load ParaView-server-osmesa/5.10.1-foss-2022a-mpi
srun pvserver --force-offscreen-rendering
```
As set up, the script will request access to exclusive use of a CPU node for one hour.
We can use `squeue` to check when the job starts running.
Once it is running, look at `NODELIST` to find out which node the server is running on (let's say we saw it is running on `tcn237`).
On your *local* machine, now run the following ssh command in a new terminal
```
ssh -L 11111:tcn237:11111 username@snellius.surf.nl
```
Here, the 5 digit port number is somewhat arbitrary, but the node must match that the server is running on.

ParaView must be installed on your machine before you can connect.
Specifically, you need the **_same version_** of ParaView as is used on the HPC.
For Snellius, the software stack `2022` currently uses version `5.10`.

Finally, with the correct version of ParaView installed, you can connect to the remote server with the following steps:

1. *File* -> *Connect*
2. Optionally add the server with Add Server (host `localhost`, port `11111`, server type `Client / Server`, Startup Type `Manual`)
3. Select `cs://localhost:11111` in the Choose Server Configuration dialog
4. Click *Connect*

## Using the AFiDTools module to make the simulation output readable from ParaView

We can load the data into ParaView by using XDMF.
This data format uses XML to inform ParaView about the structure of the HDF5 files that we want to load, so we need to produce these XML files.
Python functions are provided to produce these XML files for the 2D slices stored in `outputdir/flowmov` and the 3D fields stored in `outputdir/fields`.
A (somewhat poorly maintained) overview of XDMF and the structure of these XML files can be found on the [XDMF website](https://www.xdmf.org/index.php/XDMF_Model_and_Format).

### Viewing 3D fields
Performing volume rendering in ParaView is most efficient when the data is presented to ParaView on a uniform grid.
We therefore provide post-processing functions in the `afidtools` Python module to interpolate the 3D fields stored in `outputdir/fields` and produce new 3D fields which get stored in `outputdir/viz`.

For visualisation purposes, these fields do not need to be as high resolution as in the simulations.
We therefore downscale the fields as part of this interpolation.
The downscaling also helps prevent issues with memory during visualisation.

The following Python script can be used to generate the downscaled 3D fields, and the corresponding XDMF files.
Here, the fields are downscaled by a factor of 4, whereas the default is 2 (if `scale` is not specified).

```python
import afidtools as afid
folder = "/path/to/my/simulation"
vars = ["vx", "vy", "vz", "temp", "sal"]
for var in vars:
    afid.interpolate_field_to_uniform(folder, var, scale=4)
    print("Interpolated "+var+" to uniform")
    afid.generate_uniform_xmf(folder, var, scale=4)
    print("Created xmf for "+var)
```
For each variable, this creates a `var_fields.xmf` file which is what should be opened in ParaView.

### 2D slices from `flowmov`
Although viewing the 2D slices can be easily done in [Jupyter notebooks](jupyter.md), it is also possible to load the planes in ParaView.
This also has the benefit of plotting the various slices in 3-D space.
The Python function `afid.generate_cut_xmf(folder, plane)` produces an XML file for all the simulation variables in each plane slice.
These are stored in `outputdir/flowmov` and can be opened in ParaView using the XDMF Reader.

<!-- ### 3D fields from the `fields` directory
The grids used by AFiD can be represented in three ways in XDMF through the `TopologyType` specified in the XML:

- 3DSMesh: This is the most generic form of a structured grid. It requires three 3D arrays specifying the coordinates $x$, $y$, and $z$ for each point of the arrays describing the flow field. These grids are not output by AFiD, but can be generated by the Python function `create_3D_grids` from the AFiDTools module. Volume rendering is *possible* in ParaView when this type of grid is used, but it is very slow.
- 3DRectMesh: This type describes a structured grid with perpendicular axes, as is used by AFiD. The grids from `outputdir/cordin_info.h5` can be used to describe the coordinates of the arrays. However, volume rendering is not possible in ParaView using this type of mesh.
- 3DCoRectMesh: This describes a grid with uniform grid spacing in each direction. This is the best type of grid for high-performance volume rendering. However, if grid stretching is used by the simulation, this will distort the visualisation of the fields. We aim to soon produce a small Fortran program using the interpolation routines of MuRPhFi to convert 3D fields to uniform grids for visualisation. -->