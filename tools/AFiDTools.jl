module AFiDTools

export
    read_grid,
    read_mean,
    read_cut,
    continua_master_from_input,
    write_continua,
    tmean,
    xmean,
    ddx,
    coarsen,
    refine

using
    DelimitedFiles,
    HDF5,
    Interpolations

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
    varlist = keys(fileid)
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

"""
    function tmean(var::Array{Float64})
    
Computes the time average of data obtained from the `read_mean()` function.
Passing a 2D array produces a time series `Vector`, whereas passing a 1D
array gives the mean and standard deviation from the time series as a named tuple.
"""
function tmean(var::Matrix{Float64})
    if size(var)[2] < 1350
        A = mean(var[:,150:end], dims=2)
    else
        A = mean(var[:,end-1250:end], dims=2)
    end
    return A
end
function tmean(var::Vector{Float64})
    if length(var) < 1350
        tmp = var[150:end]
    else
        tmp = var[end-1250:end]
    end
    return (mean=mean(tmp), std=std(tmp))
end

"""
    function xmean(var::Array{Float64}, grid)

Computes the domain average of the mean profile time series data
obtained from the `read_mean()` function.
Requires an input `grid` that matches the output from `read_grid()`.
"""
function xmean(var::Array{Float64}, grid)
    if size(var)[1]==length(grid.xm)
        dx = grid.xc[2:end] - grid.xc[1:end-1]
    else
        dx = grid.xcr[2:end] - grid.xcr[1:end-1]
    end
    mvar = sum(dx .* var, dims=1)[:]
    return mvar
end

"""
    function ddx(f::Array{Float64}, x::Vector{Float64}; BClow::Float64=0.0, BChigh::Float64=0.0)
Computes the x-derivative of the mean profile quantities generated by
`read_mean()` when the relevant x-coordinates `x` are passed to the function.
Optional arguments for the lower and upper boundary values can be passed as
`BClow` and `BChigh`. By default, these are set to zero.
"""
function ddx(f::Matrix{Float64}, x::Vector{Float64}; BClow::Float64=0.0, BChigh::Float64=0.0)
    df = zeros(size(f))
    df[1,:] = (f[2,:] .- BClow)./(x[2] - 0)
    df[2:end-1,:] = (f[3:end,:] .- f[1:end-2,:])./(x[3:end] .- x[1:end-2])
    df[end,:] = (BChigh .- f[end-1,:])./(1 - x[end-1])
    return df
end
function ddx(f::Vector{Float64}, x::Vector{Float64}; BClow::Float64=0.0, BChigh::Float64=0.0)
    df = zeros(size(f))
    df[1] = (f[2] - BClow)./(x[2] - 0)
    df[2:end-1] = (f[3:end] - f[1:end-2])./(x[3:end] - x[1:end-2])
    df[end] = (BChigh - f[end-1])./(1 - x[end-1])
    return df
end

"""
    function coarsen(var::Array{Float64}, grid)
Interpolates a `read_mean()` quantity `var` from the base grid `grid.xm`
to the refined grid `grid.xmr`.
"""
function coarsen(var::Matrix{Float64}, grid)
    nodes = (grid.xmr,)
    Ns = size(var)[2]
    varc = zeros(length(grid.xm),Ns)
    for i=1:Ns
        itp = interpolate(nodes, var[:,i], Gridded(Linear()))
        varc[:,i] = itp.(grid.xm)
    end
    return varc
end
function coarsen(var::Vector{Float64}, grid)
    nodes = (grid.xmr,)
    itp = interpolate(nodes, var[:], Gridded(Linear()))
    varc = itp.(grid.xm)
    return varc
end
"""
    function refine(var::Array{Float64}, grid)
Interpolates a `read_mean()` quantity `var` from the refined grid `grid.xmr`
to the base grid `grid.xm`.
"""
function refine(var::Matrix{Float64}, grid)
    x = vcat(0, grid.xm, 1)
    nodes = (x,)
    Ns = size(var)[2]
    varr = zeros(length(grid.xmr),Ns)
    tmpvar = vcat(transpose(zeros(Ns)), var, transpose(zeros(Ns)))
    for i=1:Ns
        itp = interpolate(nodes, tmpvar[:,i], Gridded(Linear()))
        varr[:,i] = itp.(grid.xmr)
    end
    return varr
end
function refine(var::Vector{Float64}, grid)
    x = vcat(0, grid.xm, 1)
    nodes = (x,)
    tmpvar = vcat(0, var[:], 0)
    itp = interpolate(nodes, tmpvar, Gridded(Linear()))
    varr = itp.(grid.xmr)
    return varr
end

end # module