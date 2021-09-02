# 3D visualization with Paraview

## Connecting remotely to a ParaView HPC session

Many HPC systems are set up with ParaView installed, and allow you to make use of their superior computational power to load and render 3-D fields.
This can simply be done by running ParaView on the HPC and using a VNC to access the GUI on the remote system.
But ParaView also allows you to connect a local ParaView installation to a server running on the HPC.
This way, the remote system only focuses on rendering the data and does not need to generate the GUI and stream it via VNC.

The Cartesius [visualisation documentation](https://userinfo.surfsara.nl/systems/cartesius/paraview_client_server) by SURFsara provides an excellent walkthrough for setting up this ParaView connection.
On the HPC side, you need access to the GPU visualisation nodes, from which you can run the following commands
```
module load 2020
module load remotevis
gcn_paraview -w 1h -n 2 2
```
This will submit a job requesting 2 GPU nodes with 2 processes each for 1 hour.
Once the job is ready, information on how to connect will appear.
On your *local* machine, now run a command such as
```
ssh -L 11111:gcn4:11111 username@vis.cartesius.surfsara.nl
```
The exact command you should run, specifying the host node and ports, will be shown by the output of `gcn_paraview` on the HPC.

ParaView must be installed on your machine before you can connect.
Specifically, you need the **_same version_** of ParaView as is used by `gcn_paraview` on the HPC.
For Cartesius, the software stack `2020` currently uses version `5.8`.

**Note**: the `2021` software stack that will be used on the new Snellius HPC does not currently have a `remotevis` module but does have version `5.9` of ParaView available, so we should expect the new `remotevis` (or equivalent) to use this version.

Finally, with the correct version of ParaView installed, you can connect to the remote server with the following steps:

1. *File* -> *Connect*
2. Optionally add the server with Add Server (host `localhost`, port `11111`, server type `Client / Server`, Startup Type `Manual`)
3. Select `cs://localhost:11111` in the Choose Server Configuration dialog
4. Click *Connect*

## Using the AFiDTools Python module to make the simulation output readable from ParaView