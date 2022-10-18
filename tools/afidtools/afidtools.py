# Tools for pre- and post-processing AFiD data

import h5py
import numpy as np

class Grid:
    """
    Class containing the coordinates of all grids used in the simulation
    """
    def __init__(self, folder):
        with h5py.File(folder+"/outputdir/cordin_info.h5","r") as f:
            self.xm = f["xm"][()]
            self.xc = f["xc"][()]
            self.ym = f["ym"][()]
            self.zm = f["zm"][()]
            if "xmr" in list(f.keys()):
                self.xmr = f["xmr"][()]
                self.xcr = f["xcr"][()]
                self.ymr = f["ymr"][()]
                self.zmr = f["zmr"][()]
            else:
                self.xmr, self.xcr = np.array([]), np.array([])
                self.ymr, self.zmr = np.array([]), np.array([])
        if self.ym.size > 1:
            dy = self.ym[1] - self.ym[0]
            self.yc = np.arange(0, dy*(self.ym.size+1), dy)
        else: self.yc = np.array([0, 2*self.ym[0]])
        if self.ymr.size > 1:
            dyr = self.ymr[1] - self.ymr[0]
            self.ycr = np.arange(0, dyr*(self.ymr.size+1), dyr)
        elif self.ymr.size==1:
            self.ycr = np.array([0,2*self.ymr[0]])
        else: self.ycr = []
        if self.zm.size > 1:
            dz = self.zm[1] - self.zm[0]
            self.zc = np.arange(0, dz*(self.zm.size+1), dz)
        else: self.zc = np.array([0, 2*self.zm[0]])
        if self.zmr.size > 1:
            dzr = self.zmr[1] - self.zmr[0]
            self.zcr = np.arange(0, dzr*(self.zmr.size+1), dzr)
        elif self.zmr.size==1:
            self.zcr = np.array([0,2*self.zmr[0]])
        else: self.zcr = []

class InputParams:
    """
    Class containing the simulation parameters specified by
    the input file `bou.in`.
    """
    def __init__(self, folder):
        fname = folder+"/bou.in"
        with open(fname,"r") as f:
            bou = f.readlines()
        if len(bou)==55:
            # Current format
            self.nxm, self.nym, self.nzm, self.nsst = [
                int(n) for n in bou[2].split()
            ]
            self.multires, self.nxmr, self.nymr, self.nzmr = [
                int(n) for n in bou[6].split()
            ]
            self.flagsal, self.flagpf = [
                n!="0" for n in bou[10].split()
            ]
            self.tout, self.tframe = [
                float(n) for n in bou[22].split()[:2]
            ]
            self.save_3D = float(bou[22].split()[-1])
            self.alx3, self.ylen, self.zlen = [
                float(n) for n in bou[26].split()
            ]
            self.istr3 = int(bou[30].split()[0])
            self.str3 = float(bou[30].split()[1])
            self.istr3r = int(bou[30].split()[2])
            self.RayT, self.PraT, self.RayS, self.PraS = [
                float(n) for n in bou[34].split()[:-1]
            ]
            self.FFscaleS = bou[34].split()[-1]=="1"
            self.inslwS, self.inslwN, self.TfixS, self.TfixN, self.SfixS, self.SfixN = [
                n!="0" for n in bou[42].split()
            ]
            self.active_T, self.active_S = [
                n!="0" for n in bou[46].split()[:-1]
            ]
            self.gAxis = int(bou[46].split()[-1])
            self.xplusU, self.xminusU, self.dPdy, self.dPdz = [
                float(n) for n in bou[50].replace("d","e").split()[:-1]
            ]
            self.pf_D, self.pf_A, self.pf_S, self.pf_Tm = [
                float(n) for n in bou[54].split()[:4]
            ]
            self.IBM = bou[54].split()[4]!="0"
            self.pf_IC = int(bou[54].split()[-1])
        
        elif len(bou)==24:
            self.nxm, self.nym, self.nzm = [
                int(n) for n in bou[1].split()[:3]
            ]
            self.nxmr, self.nymr, self.nzmr = [
                int(n) for n in bou[3].split()
            ]
            self.tout = float(bou[5].split()[2])
            self.FFscaleS = bou[5].split()[-1]=="0"
            self.alx3 = float(bou[7].split()[0])
            self.ylen, self.zlen = [
                float(n) for n in bou[9].split()
            ]
            self.RayT, self.PraT, self.RayS, self.PraS = [
                float(n) for n in bou[11].split()[:4]
            ]
            self.inslwS, self.inslwN, self.TfixS, self.TfixN, self.SfixS, self.SfixN = [
                n!="0" for n in bou[15].split()[:6]
            ]
            self.gAxis = int(bou[15].split()[-1])
            self.xplusU, self.xminusU, self.dPdz = [
                float(n) for n in bou[21].replace("d","e").split()
            ]
            self.tframe = float(bou[23].split()[0])
            self.flagsal = True



def read_mean(folder, varname):
    """
    Returns a 2D array containing the time-series of the wall-normal
    profile of the specified variable `varname`.
    """
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
    """
    Returns a 2D slice of sample number `idx` of the variable `var`.
    `plane` can be 'x', 'y', or 'z' and specifies which dimension
    is sliced to obtain the cut. `folder` specifies the root directory
    of the simulation.
    """
    varname = var+"/"+"%05d" % idx
    with h5py.File(folder+"/outputdir/flowmov/movie_"+plane+"cut.h5","r") as f:
        A = np.array(f[varname][()])
    return A

def read_plane_mean(folder, var, idx, plane):
    """
    Returns a 2D plane of sample number `idx` of the variable `var`.
    `plane` can be 'y' or 'z' and specifies which dimension
    is averaged over to obtain the plane. `folder` specifies the root directory
    of the simulation.
    """
    varname = var+"/"+"%05d" % idx
    with h5py.File(folder+"/outputdir/flowmov/movie_"+plane+"mean.h5","r") as f:
        A = np.array(f[varname][()])
    return A

def continua_master_from_input(folder, time=0.0):
    """
    This function generates the file `continua_master.h5` needed to run
    a simulation from the continua files. Optionally pass the `time` variable
    to set the start time value for the simulation.
    """
    # with open(folder+"/bou.in", "r") as f:
    #     for i, line in enumerate(f):
    #         if i==2:
    #             nxm, nym, nzm = [int(n) for n in line.split()[0:3]]
    #         if i==6:
    #             nxmr, nymr, nzmr = [int(n) for n in line.split()[1:4]]
    #         if i==26:
    #             ylen, zlen = [int(n) for n in line.split()[1:3]]
    #         if i==30:
    #             istr3, str3, istr3r = [float(n) for n in line.split()]
    inputs = InputParams(folder)

    with h5py.File(folder+"/outputdir/continua_master.h5","w") as f:
        f["nx"], f["ny"], f["nz"] = inputs.nxm + 1, inputs.nym + 1, inputs.nzm + 1
        if inputs.multires:
            f["nxr"], f["nyr"], f["nzr"] = inputs.nxmr + 1, inputs.nymr + 1, inputs.nzmr + 1
        f["ylen"], f["zlen"] = inputs.ylen, inputs.zlen
        f["istr3"], f["str3"] = int(inputs.istr3), inputs.str3
        if inputs.multires: f["istr3r"] = int(inputs.istr3r)
        f["time"] = time
    return

def write_continua(folder, varname, var):
    """
    Generates a continua file for the variable `varname` from the
    data in array `var`.
    """
    with h5py.File(folder+"/outputdir/continua_"+varname+".h5","w") as f:
        f["var"] = var
    return

def z_vorticity(folder, idx, grid):
    """
    Returns a cut of the z-component of vorticity, calculated by the
    derivatives of the vx and vy fields. Sample number `idx`, simulation
    directory `folder` and grid data `grid` are needed.
    """
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