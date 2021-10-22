"""
A collection of helper functions used to write XML files
enabling the reading of AFiD's HDF5 output data by ParaView

"""

import h5py
import os

from xml.etree.ElementTree import Element, SubElement
from xml.dom import minidom
from xml.etree import ElementTree
from numpy import zeros

from .afidtools import Grid

def generate_cut_xmf(folder, plane):
    """
    Generates an xmf file in the Xdmf format to allow reading of
    the movie slices in ParaView. Specify the `plane` ("x", "y", "z")
    of the slices and the `folder` containing the simulation.
    """

    # Read the grid data from the simulation
    grid = Grid(folder)
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
    grid = Grid(folder)
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
    samplist.sort()

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

    for i, j in enumerate(samplist):
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
        var_data[i].text = "fields/%05i_" % j + var +".h5:/var"

    # Convert xmf structure to string
    rough_xmf = ElementTree.tostring(Xdmf, "utf-8")

    # Format string for readability and write to file
    formatted_xmf = minidom.parseString(rough_xmf)
    with open(folder+"/outputdir/"+var+"_fields.xmf","w") as f:
        f.write(formatted_xmf.toprettyxml(indent="  "))

def generate_uniform_xmf(folder, var):
    """
    Generates an xmf file in the Xdmf format to allow reading of
    the 3D fields in ParaView. This function produces a 3D array 
    on a grid with uniform spacing, allowing for fast volume rendering.
    If grid stretching is used in the simulation, the visualisation
    will not be accurate! Specify the variable `var`
    ("vx", "vy", "vz", "temp", "sal", "phi") and the `folder`
    containing the simulation.
    """

    # Read the grid data from the simulation
    grid = Grid(folder)
    nxm, nym, nzm = grid.xm.size, grid.ym.size, grid.zm.size
    nxmr, nymr, nzmr = grid.xmr.size, grid.ymr.size, grid.zmr.size
    dx, dy, dz = 1/nxm, grid.yc[-1]/nym, grid.zc[-1]/nzm
    if nxmr!=0:
        dxr, dyr, dzr = 1/nxmr, grid.ycr[-1]/nymr, grid.zcr[-1]/nzmr

    # Store the appropriate grid sizes and names based on the variable
    fulldims = (nzm, nym, nxm+1)
    if var=="vx":
        dims = (nzm, nym, nxm+1)
    elif var in "phisal":
        dims = (nzmr, nymr, nxmr)
        fulldims = (nzmr, nymr, nxmr+1)
    else:
        dims = (nzm, nym, nxm)
    
    # Collect indices of saved fields
    samplist = []
    for file in os.listdir(folder+"/outputdir/fields"):
        if var in file:
            samplist.append(int(file[:5]))
    samplist.sort()

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

    fields, geom, origin_data, step_data = (), (), (), ()
    var_att, var_slab, slab_data, var_data = (), (), (), ()

    for i, j in enumerate(samplist):
        fields = fields + (SubElement(time_series, "Grid", attrib={
            "Name":"field", "GridType":"Uniform"
        }),)

        SubElement(fields[i], "Time", attrib={"Value":"%.1f" % float(t_field*i)})

        SubElement(fields[i], "Topology", attrib={
            "TopologyType":"3DCoRectMesh", "Dimensions":"%i %i %i" % dims
        })
        
        geom = geom + (SubElement(fields[i], "Geometry", attrib={
            "GeometryType":"ORIGIN_DXDYDZ"
        }),)
        origin_data = origin_data + (SubElement(geom[i], "DataItem", attrib={"Dimensions":"3"}),)
        if var=="vx":
            origin_data[i].text = 3*"%.5f " % (0, dy/2, dz/2)
        elif var=="vy":
            origin_data[i].text = 3*"%.5f " % (dx/2, 0, dz/2)
        elif var=="vz":
            origin_data[i].text = 3*"%.5f " % (dx/2, dy/2, 0)
        elif var=="temp":
            origin_data[i].text = 3*"%.5f " % (dx/2, dy/2, dz/2)
        else:
            origin_data[i].text = 3*"%.5f " % (dxr/2, dyr/2, dzr/2)
        step_data = step_data + (SubElement(geom[i], "DataItem", attrib={"Dimensions":"3"}),)
        if var in "phisal":
            step_data[i].text = 3*"%.5f " % (dxr, dyr, dzr)
        else:
            step_data[i].text = 3*"%.5f " % (dx, dy, dz)

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
        var_data[i].text = "fields/%05i_" % j + var +".h5:/var"

    # Convert xmf structure to string
    rough_xmf = ElementTree.tostring(Xdmf, "utf-8")

    # Format string for readability and write to file
    formatted_xmf = minidom.parseString(rough_xmf)
    with open(folder+"/outputdir/"+var+"_fields.xmf","w") as f:
        f.write(formatted_xmf.toprettyxml(indent="  "))


def generate_rawfield_xmf(folder, var):
    """
    Generates an xmf file in the Xdmf format to allow reading of
    the 3D fields in ParaView. Creates a 3D structured mesh that requires
    the `grids.h5` file to be produced by the helper function 
    `create_3D_grids`. Volume rendering is possible, but very slow with
    this output. Specify the variable `var`
    ("vx", "vy", "vz", "temp", "sal", "phi") and the `folder`
    containing the simulation.
    """

    # Read the grid data from the simulation
    grid = Grid(folder)
    nxm, nym, nzm = grid.xm.size, grid.ym.size, grid.zm.size
    nxmr, nymr, nzmr = grid.xmr.size, grid.ymr.size, grid.zmr.size

    # Store the appropriate grid sizes and names based on the variable
    if var in "phisal":
        dims = (nzmr, nymr, nxmr+1)
        xx, yy, zz = "xcr", "ymr", "zmr"
    else:
        dims = (nzm, nym, nxm+1)
        xx, yy, zz = "xc", "ym", "zm"
    
    # Collect indices of saved fields
    samplist = []
    for file in os.listdir(folder+"/outputdir/fields"):
        if var in file:
            samplist.append(int(file[:5]))
    samplist.sort()

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
    var_att, var_data = (), ()

    for i, j in enumerate(samplist):
        fields = fields + (SubElement(time_series, "Grid", attrib={
            "Name":"field", "GridType":"Uniform"
        }),)

        SubElement(fields[i], "Time", attrib={"Value":"%.1f" % float(t_field*i)})

        SubElement(fields[i], "Topology", attrib={
            "TopologyType":"3DSMesh", "Dimensions":"%i %i %i" % dims
        })
        
        geom = geom + (SubElement(fields[i], "Geometry", attrib={
            "GeometryType":"X_Y_Z"
        }),)
        zdata = zdata + (SubElement(geom[i], "DataItem", attrib={
            "Dimensions":"%i %i %i" % dims, "Format":"HDF"
        }),)
        zdata[i].text = "fields/grids.h5:/"+xx
        ydata = ydata + (SubElement(geom[i], "DataItem", attrib={
            "Dimensions":"%i %i %i" % dims, "Format":"HDF"
        }),)
        ydata[i].text = "fields/grids.h5:/"+yy
        xdata = xdata + (SubElement(geom[i], "DataItem", attrib={
            "Dimensions":"%i %i %i" % dims, "Format":"HDF"
        }),)
        xdata[i].text = "fields/grids.h5:/"+zz
        
        var_att = var_att + (SubElement(fields[i], "Attribute", attrib={
            "Name":var, "AttributeType":"Scalar", "Center":"Node"
        }),)
        var_data = var_data + (SubElement(var_att[i], "DataItem", attrib={
            "Dimensions":"%i %i %i" % dims, "Format":"HDF"
        }),)
        var_data[i].text = "fields/%05i_" % j + var +".h5:/var"

    # Convert xmf structure to string
    rough_xmf = ElementTree.tostring(Xdmf, "utf-8")

    # Format string for readability and write to file
    formatted_xmf = minidom.parseString(rough_xmf)
    with open(folder+"/outputdir/"+var+"_fields.xmf","w") as f:
        f.write(formatted_xmf.toprettyxml(indent="  "))

def create_3D_grids(folder):
    """
    Creates a grids.h5 file containing 3D arrays with the
    locations of the xc, ym, zm grids at each spatial point.
    """
    grid = Grid(folder)
    dims = (grid.zm.size, grid.ym.size, grid.xc.size)
    tmp = zeros(dims)
    for i in range(dims[0]):
        tmp[i,:,:] = grid.zm[i]
    with h5py.File(folder+"/outputdir/fields/grids.h5","a") as f:
        f["zm"] = tmp
    for j in range(dims[1]):
        tmp[:,j,:] = grid.ym[j]
    with h5py.File(folder+"/outputdir/fields/grids.h5","a") as f:
        f["ym"] = tmp
    for k in range(dims[2]):
        tmp[:,:,k] = grid.xc[k]
    with h5py.File(folder+"/outputdir/fields/grids.h5","a") as f:
        f["xc"] = tmp
