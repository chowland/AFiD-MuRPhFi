module AFiDTools

export
    read_grid,
    read_mean,
    read_cut,
    continua_master_from_input,
    write_continua

using
    DelimitedFiles,
    HDF5

"""
    function read_grid(folder::String)
    
Returns an named tuple containing staggered coordinates `xm`, `ym`, `zm`,
cell-edge coordinates e.g. `xc`, and refined grid coordinates e.g. `xmr`
"""
function read_grid(folder::String)
    filename = folder*"/outputdir/cordin_info.h5"
    fileid = h5open(filename,"r")
    xm, xc = read(fileid["xm"]), read(fileid["xc"])
    ym, zm = read(fileid["ym"]), read(fileid["zm"])
    if "xmr" in keys(fileid)
        xmr, xcr = read(fileid["xmr"]), read(fileid["xcr"])
        ymr, zmr = read(fileid["ymr"]), read(fileid["zmr"])
    else
        xmr, xcr, ymr, zmr = [], [], [], []
    end
    close(fileid)
    if length(ym)>1
        dy = ym[2] - ym[1]
        yc = 0:dy:dy*length(ym)
    else
        yc = []
    end
    if length(ymr)>1
        dyr = ymr[2] - ymr[1]
        ycr = 0:dyr:dyr*length(ymr)
    else
        ycr = []
    end
    if length(zm)>1
        dz = zm[2] - zm[1]
        zc = 0:dz:dz*length(zm)
    else
        zc = []
    end
    if length(zmr)>1
        dzr = zmr[2] - zmr[1]
        zcr = 0:dzr:dzr*length(zmr)
    else
        zcr = []
    end
    return (xm=xm, xmr=xmr, xc=xc, xcr=xcr, ym=ym, ymr=ymr, yc=yc, ycr=ycr, zm=zm, zmr=zmr, zc=zc, zcr=zcr)
end

"""
    function read_mean(folder::String,varname::String)
    
Returns a 2D array containing a time series of mean profiles f̄(x,t) of `varname`
from the simulation output in `folder`.
"""
function read_mean(folder::String,varname::String)
    filename = folder*"/outputdir/means.h5"
    fileid = h5open(filename,"r")
    varlist = keys(filied)
    if varname in varlist
        list = keys(fileid[varname])
        Nsamp = size(list)[1]
        nxx = length(fileid[varname*"/"*list[1]])
        var = zeros((nxx,Nsamp))
        i = 1
        for num in list
            var[:,i] = read(fileid[varname*"/"*num])
            i = i + 1;
        end
    else
        varstring = string()
        for var in varlist
            varstring *= var*", "
        end
        error(varname*" not a variable in means.h5. Valid options are "*varstring)
    end
    close(fileid)
    return var
end

"""
    function read_cut(folder::String,varname::String,idx::Int,plane::String)
    
Returns an array slice of `varname` ("vx", "vy", "vz", "temp", "sal", or "phi")
from a `plane` ("x", "y", or "z") at time index `idx` from the simulation
output in `folder`.
"""
function read_cut(folder::String,varname::String,idx::Int,plane::String)
    filename = folder*"/outputdir/flowmov/movie_"*plane*"cut.h5"
    fileid = h5open(filename,"r")
    var = read(fileid[varname*"/"*string(idx,pad=5)])
    close(fileid)
    return var
end

function continua_master_from_input(folder::String)
    input_data = readdlm(folder*"/bou.in")
    nx, ny, nz = input_data[3,1:3] .+ 1
    nxr, nyr, nzr = input_data[6,2:4] .+ 1
    ylen, zlen = input_data[21,2:3]
    istr3, str3, istr3r = input_data[24,1:3]
    time = 0.0
    fname = folder*"/outputdir/continua_master.h5"
    fid = h5open(fname,"w")
    fid["nx"], fid["ny"], fid["nz"] = nx, ny, nz
    fid["nxr"], fid["nyr"], fid["nzr"] = nxr, nyr, nzr
    fid["ylen"], fid["zlen"] = ylen, zlen
    fid["istr3"], fid["str3"], fid["istr3r"] = istr3, str3, istr3r
    fid["time"] = time
    close(fid)
    return
end

function write_continua(folder::String,varname::String,var::Array{Float64, 3})
    fname = folder*"/outputdir/continua_"*varname*".h5"
    fid = h5open(fname,"w")
    fid["var"] = var
    close(fid)
    return
end

end # module