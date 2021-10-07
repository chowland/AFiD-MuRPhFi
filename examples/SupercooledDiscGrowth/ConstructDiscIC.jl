# Load required packages
using HDF5, DelimitedFiles, SpecialFunctions

# Read input file
bou = readdlm("bou.in");

# Store variables from input needed for writing continua_master
istr3, str3, istr3r = bou[24,1:3];
nx, ny, nz = bou[3,1:3] .+ 1;
nxr,nyr,nzr= bou[6,2:4] .+ 1;
time = 0.0;
ylen, zlen = bou[21,2:3];

# Initial radius of solid disc
R0 = 0.1

# Check if output directory already exists, create if not
"outputdir" ∉ readdir(".") ? mkdir("outputdir") : nothing

# Write continua_master file
fname = "outputdir/continua_master.h5"
f = h5open(fname, "w")
f["istr3"] = istr3
f["istr3r"]= istr3r
f["nx"] = nx
f["nxr"]= nxr
f["ny"] = ny
f["nyr"]= nyr
f["nz"] = nz
f["nzr"]= nzr
f["str3"] = str3
f["time"] = time
f["ylen"] = ylen
f["zlen"] = zlen
close(f)

# Create velocity continua files as zero arrays
v = zeros(nx,ny-1,nz-1)
fname = "outputdir/continua_vx.h5"
f = h5open(fname, "w")
f["var"] = v
close(f)
fname = "outputdir/continua_vy.h5"
f = h5open(fname, "w")
f["var"] = v
close(f)
fname = "outputdir/continua_vz.h5"
f = h5open(fname, "w")
f["var"] = v
close(f)

# Create uniform grids in periodic dimensions
yc, ycr = LinRange(0, ylen, ny), LinRange(0, ylen, nyr);
zc, zcr = LinRange(0, zlen, nz), LinRange(0, zlen, nzr);
ym, ymr = 0.5*(yc[2:end] + yc[1:end-1]), 0.5*(ycr[2:end] + ycr[1:end-1]);
zm, zmr = 0.5*(zc[2:end] + zc[1:end-1]), 0.5*(zcr[2:end] + zcr[1:end-1]);

# Create wall-normal grid based on input parameters
function make_x_grid(nx::Int64, istr3::Int64, str3::Float64)
    if istr3==0
        xc = LinRange(0.0, 1.0, nx)
    elseif istr3==4
        xc = zeros(nx)
        for k=2:nx
            z2dp = (2*k - nx - 1)/(nx - 1)
            xc[k] = 0.5*(1 + tanh(str3*z2dp)/tanh(str3))
        end
    elseif istr3==6
        nclip=Int(str3)
        nxmo = nx + 2*nclip
        etazm = cos.(π*(1:nxmo - 0.5)/nxmo)
        etaz = etazm[nclip+1:nclip+nx]
        etaz = etaz ./ (0.5*(etaz[1] - etaz[nx]))
        xc = zeros(nx)
        xc[2:nx-1] = 0.5*(1 .- etaz[2:nx-1])
        xc[end] = 1.0
    elseif istr3==7
        nclip=Int(str3)
        nxmo = nx + nclip
        etazm = cos.(π*(1:nxmo)/nxmo/2.0)
        etaz = etazm[nclip+1:nclip+nx]
        etaz = etaz ./ etaz[1]
        xc = zeros(nx)
        xc[2:nx-1] = 0.5*(1 .- etaz[2:nx-1])
        xc[end] = 1.0
    end
    xm = 0.5*(xc[2:end] + xc[1:end-1])
    return xc, xm
end

# Construct wall-normal grid using above functions
xc, xm = make_x_grid(nx, istr3, str3);
xcr, xmr = make_x_grid(nxr, istr3r, str3);

# Create phase field variable condition: 1 in disc, 0 outside
ϕ = zeros(nxr,nyr-1,nzr-1);
ϵ = 1/(nxr-1)
for i=1:nzr-1
    for j=1:nyr-1
        for k=1:nxr-1
            r = sqrt((xmr[k] - 0.5)^2 + (ymr[j] - ylen/2)^2)
            ϕ[k,j,i] = 0.5*(1 - tanh((r - R0)/2/ϵ))
        end
    end
end

# Write phase field variable to continua file
fname = "outputdir/continua_phi.h5"
f = h5open(fname, "w")
f["var"] = ϕ
close(f)

# Read Stefan number from input data
Stefan = bou[42,3];

# Define function for radial temperature profile
function F(x)
    f = -expinti(-x^2/4)
    return f
end

# Use Newton iteration to find appropriate value for Lambda
function iterate_for_lambda(x, S, rtol)
    while true
        g = F(x) - 2/S/x^2*exp(-x^2 / 4)
        g′= exp(-x^2 / 4)*(8/S/x^3 + 2*(1-S)/S/x)
        println("g($x) = $g")
        Δx = -g/g′
        x +=Δx
        abs(Δx) ≤ abs(x)*rtol && break
    end
    return x
end

rtol = 1e-7   # Set tolerance for integral calculation
Λ = iterate_for_lambda(0.5,Stefan,1e-7)

# Read physical parameters from input data
RaT, PrT = bou[27,1:2]   # Rayleigh and Prandtl numbers
Pe = sqrt(RaT*PrT)       # Peclet number

t0 = Pe*(R0/Λ)^2        # Time value needed for similarity solution
println("t0 = ",t0)

# Create temperature field
T = zeros(nx,ny-1,nz-1);
for i=1:nz-1
    for j=1:ny-1
        for k=1:nx-1
            r = sqrt((xm[k] - 0.5)^2 + (ym[j] - ylen/2)^2)
            if r < R0
                T[k,j,i] = 1.0
            else
                T[k,j,i] = F(r*sqrt(Pe/t0))/F(Λ)
            end
        end
    end
end

# Write temperature field to continua file
fname = "outputdir/continua_temp.h5"
f = h5open(fname, "w")
f["var"] = T
close(f)