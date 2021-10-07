# Tools for pre- and post-processing AFiD data

import h5py
import numpy as np

class Grid:
    def __init__(self, xm, xmr, xc, xcr, ym, ymr, yc, ycr, zm, zmr, zc, zcr):
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
        self.zc = zc
        self.zcr= zcr

def read_grid(folder):
    with h5py.File(folder+"/outputdir/cordin_info.h5","r") as f:
        xm = f["xm"][()]
        xc = f["xc"][()]
        ym = f["ym"][()]
        zm = f["zm"][()]
        if "xmr" in list(f.keys()):
            xmr = f["xmr"][()]
            xcr = f["xcr"][()]
            ymr = f["ymr"][()]
            zmr = f["zmr"][()]
        else:
            xmr, xcr = np.array([]), np.array([])
            ymr, zmr = np.array([]), np.array([])
    if ym.size > 1:
        dy = ym[1] - ym[0]
        yc = np.arange(0, dy*(ym.size+1), dy)
    else: yc = np.array([0, 2*ym[0]])
    if ymr.size > 1:
        dyr = ymr[1] - ymr[0]
        ycr = np.arange(0, dyr*(ymr.size+1), dyr)
    elif ymr.size==1:
        ycr = np.array([0,2*ymr[0]])
    else: ycr = []
    if zm.size > 1:
        dz = zm[1] - zm[0]
        zc = np.arange(0, dz*(zm.size+1), dz)
    else: zc = np.array([0, 2*zm[0]])
    if zmr.size > 1:
        dzr = zmr[1] - zmr[0]
        zcr = np.arange(0, dzr*(zmr.size+1), dzr)
    elif zmr.size==1:
        zcr = np.array([0,2*zmr[0]])
    else: zcr = []
    return Grid(xm, xmr, xc, xcr, ym, ymr, yc, ycr, zm, zmr, zc, zcr)

def read_mean(folder, varname):
    with h5py.File(folder+"/outputdir/means.h5","r") as f:
        samplist = list(f[varname].keys())
        Nsamp = len(samplist)
        nx = f[varname+"/"+samplist[0]].size
        var = np.zeros((nx,Nsamp))
        i = 0
        for num in samplist:
            var[:,i] = f[varname+"/"+num][()]
            i = i + 1
    return var

def read_cut(folder, var, idx, plane):
    varname = var+"/"+"%05d" % idx
    with h5py.File(folder+"/outputdir/flowmov/movie_"+plane+"cut.h5","r") as f:
        A = np.array(f[varname][()])
    return A

def continua_master_from_input(folder):
    with open(folder+"/bou.in", "r") as f:
        for i, line in enumerate(f):
            if i==2:
                nxm, nym, nzm = [float(n) for n in line.split()[0:3]]
            if i==6:
                nxmr, nymr, nzmr = [float(n) for n in line.split()[1:4]]
            if i==26:
                ylen, zlen = [float(n) for n in line.split()[1:3]]
            if i==30:
                istr3, str3, istr3r = [float(n) for n in line.split()]
        time = 0.0
        with h5py.File(folder+"/outputdir/continua_master.h5","w") as f:
            f["nx"], f["ny"], f["nz"] = nxm + 1, nym + 1, nzm + 1
            f["nxr"], f["nyr"], f["nzr"] = nxmr + 1, nymr + 1, nzmr + 1
            f["ylen"], f["zlen"] = ylen, zlen
            f["istr3"], f["str3"], f["istr3r"] = istr3, str3, istr3r
            f["time"] = time
    return

def write_continua(folder, varname, var):
    with h5py.File(folder+"/outputdir/continua_"+varname+".h5","w") as f:
        f["var"] = var
    return

def z_vorticity(folder, idx, grid):
    vx = read_cut(folder, "vx", idx, "z")
    vy = read_cut(folder, "vy", idx, "z")

    dvxdy, dvydx = np.zeros(vx.shape), np.zeros(vy.shape)
    dy = (grid.ym[1] - grid.ym[0])

    dvxdy[1:-1,:] = (vx[2:,:] - vx[:-2,:])/2/dy
    dvxdy[0,:] = (vx[1,:] - vx[-1,:])/2/dy
    dvxdy[-1,:] = (vx[0,:] - vx[-2,:])/2/dy

    dvydx[:, 1:-1] = (vy[:,2:] - vy[:,:-2])/(grid.xm[2:] - grid.xm[:-2])
    dvydx[:, 0] = (vy[:,1] - 0)/(grid.xm[1] - 0)
    dvydx[:, -1] = (0 - vy[:,-2])/(1 - grid.xm[-2])

    dyu = 0.5*(dvxdy[:,1:] + dvxdy[:,:-1])

    dxv = np.zeros(vy.shape)
    dxv[:-1,:] = 0.5*(dvydx[1:,:] + dvydx[:-1,:])
    dxv[-1,:] = 0.5*(dvydx[0,:] + dvydx[-1,:])

    ζ = dxv - dyu

    return ζ