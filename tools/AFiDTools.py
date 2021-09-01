# Tools for pre- and post-processing AFiD data

import h5py
import numpy as np
import os

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
    else: yc = []
    if ymr.size > 1:
        dyr = ymr[1] - ymr[0]
        ycr = np.arange(0, dyr*(ymr.size+1), dyr)
    else: ycr = []
    if zm.size > 1:
        dz = zm[1] - zm[0]
        zc = np.arange(0, dz*(zm.size+1), dz)
    else: zc = []
    if zmr.size > 1:
        dzr = zmr[1] - zmr[0]
        zcr = np.arange(0, dzr*(zmr.size+1), dzr)
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

from xml.etree.ElementTree import Element, SubElement, Comment
from xml.dom import minidom
from xml.etree import ElementTree

def generate_cut_xmf(folder, plane):
    """
    Generates an xmf file in the Xdmf format to allow reading of
    the movie slices in ParaView. Specify the `plane` ("x", "y", "z")
    of the slices and the `folder` containing the simulation.
    """

    # Read the grid data from the simulation
    grid = read_grid(folder)
    nxm, nym, nzm = grid.xm.size, grid.ym.size, grid.zm.size
    nxmr, nymr, nzmr = grid.xmr.size, grid.ymr.size, grid.zmr.size

    # Store the dimensions of the slices in dims and dimsr
    if plane=="x":
        dims, dimsr = (nzm, nym), (nzmr, nymr)
    elif plane=="y":
        dims, dimsr = (nzm, nxm), (nzmr, nxmr)
    elif plane=="z":
        dims, dimsr = (nym, nxm), (nymr, nxmr)
    else:
        raise NameError(plane)

    # Count the number of movie samples recorded in the time series
    with h5py.File(folder+"/outputdir/flowmov/movie_"+plane+"cut.h5","r") as f:
        Nsamp = len(list(f["temp"].keys()))
    varlist = ["vx", "vy", "vz", "temp"]

    # Check the time interval used for frame writing from bou.in
    with open(folder+"/bou.in","r") as f:
        for i, line in enumerate(f):
            if i==22:
                t_frame = float(line.split()[1])


    ### BUILD THE XMF STRUCTURE ###
    Xdmf = Element("Xdmf")

    top_domain = SubElement(Xdmf, "Domain")

    # Create Grid element to store the time series
    time_series = SubElement(top_domain, "Grid", attrib={
        "Name":"Time", "GridType":"Collection", "CollectionType":"Temporal"
    })

    # Create empty tuples to build element tree with
    cuts, geom, xdata, ydata, zdata = (), (), (), (), ()
    # Create empty tuples for each variable
    var_att, var_data = [(), (), (), ()], [(), (), (), ()]
    for i in range(Nsamp):
        cuts = cuts + (SubElement(time_series, "Grid", attrib={
            "Name":plane+"cut", "GridType":"Uniform"
        }),)

        SubElement(cuts[i], "Time", attrib={"Value":"%.1f" % float(t_frame*i)})

        SubElement(cuts[i], "Topology", attrib={
            "TopologyType":"3DRectMesh", "Dimensions":"%i %i" % dims
        })
        
        geom = geom + (SubElement(cuts[i], "Geometry", attrib={
            "GeometryType":"VXVYVZ"
        }),)
        zdata = zdata + (SubElement(geom[i], "DataItem", attrib={"Dimensions":"1"}),)
        zdata[i].text = "0.0"
        ydata = ydata + (SubElement(geom[i], "DataItem", attrib={
            "Dimensions":"%i" % dims[0], "Format":"HDF"
        }),)
        if plane=="z":
            ydata[i].text = "../cordin_info.h5:/ym"
        else:
            ydata[i].text = "../cordin_info.h5:/zm"
        xdata = xdata + (SubElement(geom[i], "DataItem", attrib={
            "Dimensions":"%i" % dims[1], "Format":"HDF"
        }),)
        if plane=="x":
            xdata[i].text = "../cordin_info.h5:/ym"
        else:
            xdata[i].text = "../cordin_info.h5:/xm"
        
        for j, var in enumerate(varlist):
            var_att[j] = var_att[j] + (SubElement(cuts[i], "Attribute", attrib={
                "Name":var, "AttributeType":"Scalar", "Center":"Node"
            }),)
            var_data[j] = var_data[j] + (SubElement(var_att[j][i], "DataItem", attrib={
                "Dimensions":"%i %i" % dims, "Format":"HDF"
            }),)
            var_data[j][i].text = "movie_"+plane+"cut.h5:/"+var+"/%05i" % i

    # Convert xmf structure to string
    rough_xmf = ElementTree.tostring(Xdmf, "utf-8")

    # Format string for readability and write to file
    formatted_xmf = minidom.parseString(rough_xmf)
    with open(folder+"/outputdir/flowmov/movie_"+plane+"cut.xmf","w") as f:
        f.write(formatted_xmf.toprettyxml(indent="  "))



def generate_field_xmf(folder, var):
    """
    Generates an xmf file in the Xdmf format to allow reading of
    the 3D fields in ParaView. Specify the variable `var`
    ("vx", "vy", "vz", "temp", "sal", "phi") and the `folder`
    containing the simulation.
    """

    # Read the grid data from the simulation
    grid = read_grid(folder)
    nxm, nym, nzm = grid.xm.size, grid.ym.size, grid.zm.size
    nxmr, nymr, nzmr = grid.xmr.size, grid.ymr.size, grid.zmr.size

    # Store the appropriate grid sizes and names based on the variable
    fulldims = (nzm, nym, nxm+1)
    xx, yy, zz = "xm", "ym", "zm"
    if var=="vx":
        dims = (nzm, nym, nxm+1)
        xx = "xc"
    elif var in "phisal":
        dims = (nzmr, nymr, nxmr)
        fulldims = (nzmr, nymr, nxmr+1)
        xx, yy, zz = "xmr", "ymr", "zmr"
    else:
        dims = (nzm, nym, nxm)
        if var=="vy":
            yy = "yc"
        elif var=="vz":
            zz = "zc"
    
    # Collect indices of saved fields
    samplist = []
    for file in os.listdir(folder+"/outputdir/fields"):
        if var in file:
            samplist.append(int(file[:5]))

    # Check the time interval used for field writing from bou.in
    with open(folder+"/bou.in","r") as f:
        for i, line in enumerate(f):
            if i==22:
                t_field = float(line.split()[2])


    ### BUILD THE XMF STRUCTURE ###
    Xdmf = Element("Xdmf")

    top_domain = SubElement(Xdmf, "Domain")

    # Create Grid element to store the time series
    time_series = SubElement(top_domain, "Grid", attrib={
        "Name":"Time", "GridType":"Collection", "CollectionType":"Temporal"
    })

    fields, geom, xdata, ydata, zdata = (), (), (), (), ()
    var_att, var_slab, slab_data, var_data = (), (), (), ()

    for i in samplist:
        fields = fields + (SubElement(time_series, "Grid", attrib={
            "Name":"field", "GridType":"Uniform"
        }),)

        SubElement(fields[i], "Time", attrib={"Value":"%.1f" % float(t_field*i)})

        SubElement(fields[i], "Topology", attrib={
            "TopologyType":"3DRectMesh", "Dimensions":"%i %i %i" % dims
        })
        
        geom = geom + (SubElement(fields[i], "Geometry", attrib={
            "GeometryType":"VXVYVZ"
        }),)
        zdata = zdata + (SubElement(geom[i], "DataItem", attrib={
            "Dimensions":"%i" % dims[2], "Format":"HDF"
        }),)
        zdata[i].text = "cordin_info.h5:/"+xx
        ydata = ydata + (SubElement(geom[i], "DataItem", attrib={
            "Dimensions":"%i" % dims[1], "Format":"HDF"
        }),)
        ydata[i].text = "cordin_info.h5:/"+yy
        xdata = xdata + (SubElement(geom[i], "DataItem", attrib={
            "Dimensions":"%i" % dims[0], "Format":"HDF"
        }),)
        xdata[i].text = "cordin_info.h5:/"+zz
        
        var_att = var_att + (SubElement(fields[i], "Attribute", attrib={
            "Name":var, "AttributeType":"Scalar", "Center":"Node"
        }),)
        var_slab = var_slab + (SubElement(var_att[i], "DataItem", attrib={
            "ItemType":"HyperSlab", "Dimensions":"%i %i %i" % dims,
            "Type": "HyperSlab"
        }),)
        slab_data = slab_data + (SubElement(var_slab[i], "DataItem", attrib={
            "Dimensions":"3 3", "Format":"XML"
        }),)
        slab_data[i].text = "0 0 0 1 1 1 %i %i %i" % dims
        var_data = var_data + (SubElement(var_slab[i], "DataItem", attrib={
            "Dimensions":"%i %i %i" % fulldims, "Format":"HDF"
        }),)
        var_data[i].text = "fields/%05i_" % i + var +".h5:/var"

    # Convert xmf structure to string
    rough_xmf = ElementTree.tostring(Xdmf, "utf-8")

    # Format string for readability and write to file
    formatted_xmf = minidom.parseString(rough_xmf)
    with open(folder+"/outputdir/"+var+"_fields.xmf","w") as f:
        f.write(formatted_xmf.toprettyxml(indent="  "))
