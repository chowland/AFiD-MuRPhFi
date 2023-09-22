"""
This script generates a movie from the slices output by a simulation.
Plots are generated in parallel and saved as individual images to a temporary
directory before being converted into an MP4.
ffmpeg must be installed for the conversion to operate.
If running this script on a HPC, do NOT run on a login node, submit it as a job.
"""

# Load the necessary modules
import multiprocessing as mp
import matplotlib.pyplot as plt
import h5py
import afidtools as afid
import os

# Locate the simulation directory
simdir = '/scratch-shared/howland/afid_RBC/gcc'
# Choose a variable to plot
var = 'sal'

# Load the grid used in the simulation
grid = afid.Grid(simdir)

# Load an initial field
T = afid.read_cut(simdir, var, 0, 'z')

# Create a figure and axis
fig, ax = plt.subplots(layout='constrained')

# Make a pcolor plot of the 2D field
pc = ax.pcolormesh(
    grid.ymr, grid.xmr, T.T,
    shading='nearest', rasterized=True,
    cmap='inferno',
    vmin=-0.5, vmax=0.5
)

# Add a colorbar
cb = fig.colorbar(pc, ax=[ax], label="$\\theta/\Delta$")

# Add axis labels and fix aspect ratio
ax.set(
    aspect='equal',
    xlabel='$x/H$',
    ylabel='$z/H$'
)

# Rescale the output image so that the resolution of the image matches the simulation
bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
nymr = grid.ymr.size
dpi = nymr/bbox.width

# Create a temporary directory to store the image output
tmpdir = './tmp'
os.mkdir(tmpdir)

# Define a function to load data from a new time and save the plot with that data
def update_and_save_plot(idx):
    T = afid.read_cut(simdir, var, idx, 'z')
    pc.set_array(T.T)
    fig.savefig(tmpdir+'/%03i.png' % idx, dpi=dpi)
    return

# Count the number of snapshots we need to save
with h5py.File(simdir+'/outputdir/flowmov/movie_zcut.h5', 'r') as f:
    nmov = len(list(f[var].keys()))

fig.canvas.draw()
fig.set_layout_engine('none')

# Create a multiprocessing pool to split up the work on 4 processors
pool = mp.Pool(4)
# Save plots for all of the 
pool.map(update_and_save_plot, range(nmov))

# Convert the images to a mp4 file (assumes ffmpeg is installed)
os.system('ffmpeg -y -i '+tmpdir+'/%03d.png -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" RBC.mp4')

# Remove the temporary directory
for filename in os.listdir(tmpdir):
    os.remove(tmpdir+'/'+filename)
os.rmdir(tmpdir)