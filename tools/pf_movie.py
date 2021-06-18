import h5py
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("ggplot")

folder = "test_runs/new_continua"

class Grid:
    def __init__(self, xm, xmr, xc, xcr, ym, ymr, yc, ycr, zm, zmr):
        self.xm = xm
        self.xmr= xmr
        self.xc = xc
        self.xcr= xcr
        self.ym = ym
        self.ymr= ymr
        self.yc = yc
        self.ycr= ycr
        self.zm = zm
        self.zmr= zmr

def read_grid(folder):
    with h5py.File(folder+"/outputdir/cordin_info.h5","r") as f:
        xm = f["xm"][()]
        xmr = f["xmr"][()]
        xc = f["xc"][()]
        xcr = f["xcr"][()]
        ym = f["ym"][()]
        ymr = f["ymr"][()]
        zm = f["zm"][()]
        zmr = f["zmr"][()]
    dy, dyr = ym[1] - ym[0], ymr[1] - ymr[0]
    yc, ycr = np.arange(0, dy*(ym.size+1), dy), np.arange(0, dyr*(ymr.size+1), dyr)
    return Grid(xm, xmr, xc, xcr, ym, ymr, yc, ycr, zm, zmr)

def read_zcut(folder, var, idx):
    varname = var+"/"+"%05d" % idx
    with h5py.File(folder+"/outputdir/flowmov/movie_zcut.h5","r") as f:
        A = f[varname][()]
    return A

def count_samples(folder):
    with h5py.File(folder+"/outputdir/flowmov/movie_zcut.h5","r") as f:
        Nsamp = len(f["temp"].keys())
    return Nsamp

grid = read_grid(folder)

print("Finished reading grid.")

T = read_zcut(folder, "temp", 0)
phi = read_zcut(folder, "phi", 0)

fig, ax = plt.subplots(figsize=(5.0, 4.0), dpi=128)
pc = ax.pcolormesh(grid.yc, grid.xc, T.T, cmap="RdBu_r", vmin=-0.6, vmax=1)
ct = ax.contour(grid.ymr, grid.xmr, phi.T, levels=[0.5], colors="k")
ax.set_aspect("equal")

Nsamp = count_samples(folder)

print("Constructed figure, starting loop.")

for i in range(Nsamp):
    ct.set_alpha(0)
    T = read_zcut(folder, "temp", i)
    phi = read_zcut(folder, "phi", i)
    pc.set_array(T.T.flatten())
    ct = ax.contour(grid.ymr, grid.xmr, phi.T, levels=[0.5], colors="k")
    fig.savefig("tmp/%05d.png" % i)
    if i % 10 == 0:
        print("Finished saving image",i)

print("Image generation complete.")