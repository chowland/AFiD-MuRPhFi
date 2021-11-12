from pandas import DataFrame
import numpy as np
import h5py
from .afidtools import Grid, InputParams, read_mean

def xmean(A, xc):
    """
    Returns the spatial average of the space-time data `A`.
    It is assumed the data of `A` lies on the cell-centred grid `xm`
    and that the grid `xc` to the function is offset from the data.
    """
    dx = xc[1:] - xc[:-1]
    Ame = np.sum(A*dx.reshape(dx.size,1), axis=0)
    return Ame

def mean_time(folder):
    """
    Returns a vector `t` that contains the time values for
    each sample in the `means.h5` file associated with the
    simulation in the directory `folder`.
    """
    inputs = InputParams(folder)
    dt = inputs.tout
    with h5py.File(folder+"/outputdir/means.h5","r") as f:
        samplist = list(f["Tbar"].keys())
    t = [float(samp)*dt for samp in samplist]
    return t

def Nusselt_numbers(folder):
    """
    Returns a DataFrame containing the time series of
    the Nusselt numbers from the upper and lower plates
    and that calculated from the scalar dissipation rates.
    Solutal Nusselt numbers are also calculated if the input
    file suggests salinity is being simulated.
    """
    inputs = InputParams(folder)
    grid = Grid(folder)
    Tbar = read_mean(folder, "Tbar")
    if inputs.RayT > 0:
        Tup, Tlo = -0.5, 0.5
    else:
        Tup, Tlo = 0.5, -0.5
    if inputs.TfixS:
        NuTlo = (Tbar[0,:] - Tlo)/grid.xm[0]
    else:
        NuTlo = np.zeros(Tbar.shape[1])
    if inputs.TfixN:
        NuTup = (Tup - Tbar[-1,:])/(1.0 - grid.xm[-1])
    else:
        NuTup = np.zeros(Tbar.shape[1])
    if inputs.RayT >0:
        NuTup, NuTlo = -NuTup, -NuTlo
    chiT = read_mean(folder, "chiT")
    if inputs.FFscaleS and inputs.flagsal:
        Re = np.sqrt(inputs.RayS/inputs.PraS)
    else:
        Re = np.sqrt(inputs.RayT/inputs.PraT)
    PeT = Re*inputs.PraT
    NuTdiss = PeT*xmean(chiT, grid.xc)

    if inputs.flagsal:
        Sbar = read_mean(folder, "Sbar")
        if inputs.RayS > 0:
            Sup, Slo = 0.5, -0.5
        else:
            Sup, Slo = -0.5, 0.5
        if inputs.SfixS:
            NuSlo = (Sbar[0,:] - Slo)/grid.xmr[0]
        else:
            NuSlo = np.zeros(Sbar.shape[1])
        if inputs.SfixN:
            NuSup = (Sup - Sbar[-1,:])/(1.0 - grid.xmr[-1])
        else:
            NuSup = np.zeros(Sbar.shape[1])
        if inputs.RayS < 0:
            NuSlo, NuSup = -NuSlo, -NuSup
        chiS = read_mean(folder, "chiS")
        PeS = Re*inputs.PraS
        NuSdiss = PeS*xmean(chiS, grid.xcr)

    t = mean_time(folder)

    df = DataFrame({
        "t": t, "NuTlo": NuTlo, "NuTup": NuTup,
        "NuTdiss": NuTdiss
    })
    if inputs.flagsal:
        df["NuSlo"] = NuSlo
        df["NuSup"] = NuSup
        df["NuSdiss"] = NuSdiss

    return df

def tmean_profile(A, t, t0=0.0):
    """
    Return the time-averaged mean profile of a space-time array `A`.
    Pass a vector `t` that describes the times recorded.
    Time averaging begins at `t0`.
    If `t0 < 0` then averaging is over the last `t0` time units.
    """
    if t0 < 0:
        t0 = t[-1] + t0
    i0 = np.nonzero(t > t0)[0][0]
    return np.mean(A[:,i0:], axis=1)

def wall_shear(folder):
    """
    Returns a DataFrame containing the time series of the
    shear Reynolds number calculated at both the upper and lower plates
    """
    inputs = InputParams(folder)
    grid = Grid(folder)
    vbar = read_mean(folder, "vybar")
    if inputs.inslwN:
        shear_up = (inputs.xplusU - vbar[-1,:])/(inputs.alx3 - grid.xm[-1])
    else:
        shear_up = np.zeros(vbar.shape[1])
    if inputs.inslwS:
        shear_lo = (vbar[0,:] - inputs.xminusU)/grid.xm[0]
    else:
        shear_lo = np.zeros(vbar.shape[1])
    if inputs.FFscaleS and inputs.flagsal:
        Re = np.sqrt(inputs.RayS/inputs.PraS)
    else:
        Re = np.sqrt(inputs.RayT/inputs.PraT)
    
    Vtau_up = np.sqrt(shear_up/Re)
    Vtau_lo = np.sqrt(shear_lo/Re)
    Retau_up, Retau_lo = Re*Vtau_up, Re*Vtau_lo

    t = mean_time(folder)

    df = DataFrame({
        "t": t, "Retau_up": Retau_up, "Retau_lo": Retau_lo
    })
    return df

def hermite_max(x4, f4):
    """
    Find the location `x` and value `f` of the turning point of the
    cubic Hermite spline generated from the vectors `x4` and `f4`
    describing a function $f(x)$ at four points.
    """
    xc, xp = x4[1:3]
    fc, fp = f4[1:3]
    dm = x4[1] - x4[0]
    dp = x4[2] - x4[1]
    dfc = (dm**2*(f4[2] - f4[1]) + dp**2*(f4[1] - f4[0]))/dp/dm/(dp+dm)
    dm = x4[2] - x4[1]
    dp = x4[3] - x4[2]
    dfp = (dm**2*(f4[3] - f4[2]) + dp**2*(f4[2] - f4[1]))/dp/dm/(dp+dm)
    a = 6*fc + dm*3*dfc - 6*fp + dm*3*dfp
    b =-6*fc - dm*4*dfc + 6*fp - dm*2*dfp
    c = dm*dfc
    tm = (-b - np.sqrt(b**2 - 4*a*c))/2/a
    tp = (-b + np.sqrt(b**2 - 4*a*c))/2/a
    if tm>=0 and tm<=1:
        t = tm
    elif tp>=0 and tp<=1:
        t = tp
    else:
        t = 0
    x = xc + (xp - xc)*t
    h00 = (1 + 2*t)*(1 - t)**2
    h10 = t*(1 - t)**2
    h01 = t**2*(3 - 2*t)
    h11 = t**2*(t - 1)
    v = h00*fc + h10*dfc*(xp-xc) + h01*fp + h11*dfp*(xp-xc)
    return x, v

def ddx(A, x, BClo=0.0, BCup=0.0):
    """
    Computes the wall-normal derivative of the space(-time) data `A`
    which sits on the grid `x`. Upper and lower boundary values, at
    x=1 and x=0 respectively, can be set by the optional arguments
    `BCup` and `BClo`.
    """
    dA = np.zeros(A.shape)
    A = np.array(A)
    x = np.array(x)
    if dA.ndim==2:
        dA[0,:] = (A[1,:] - BClo)/x[1]
        dA[1:-1,:] = (A[2:,:] - A[:-2,:])/(x[2:] - x[:-2]).reshape(x.size-2,1)
        dA[-1,:] = (BCup - A[-2,:])/(1.0 - x[-2])
    elif dA.ndim==1:
        dA[0] = (A[1] - BClo)/x[1]
        dA[1:-1] = (A[2:] - A[:-2])/(x[2:] - x[:-2])
        dA[-1] = (BCup - A[-2])/(1.0 - x[-2])
    return dA

def refine(A, grid, BClo=0.0, BCup=0.0):
    """
    Interpolates the space(-time) data `A` from the coarse grid `xm`
    to the refined grid `xmr`. Upper and lower boundary values, at
    x=1 and x=0 respectively, can be set by the optional arguments
    `BCup` and `BClo`.
    """
    xgrid = np.concatenate(([0], grid.xm, [1]))
    if A.ndim==2:
        Ar = np.zeros((grid.xmr.size, A.shape[1]))
        for i in range(A.shape[1]):
            xA = np.concatenate(([BClo], A[:,i], [BCup]))
            Ar[:,i] = np.interp(grid.xmr, xgrid, xA)
    elif A.ndim==1:
        xA = np.concatenate(([BClo], A, [BCup]))
        Ar = np.interp(grid.xmr, xgrid, xA)
    return Ar

def budget_time_series(folder):
    """
    Returns a DataFrame containing the time series of all the
    budget terms from the kinetic energy equation and the equations
    for the scalar variance.
    """
    inputs = InputParams(folder)
    ν = np.sqrt(inputs.PraS/inputs.RayS)
    κT = ν/inputs.PraT
    κS = ν/inputs.PraS
    grid = Grid(folder)
    uv = read_mean(folder, "vxvy")
    uw = read_mean(folder, "vxvz")
    v̄ = read_mean(folder, "vybar")
    w̄ = read_mean(folder, "vzbar")
    dv̄ = ddx(v̄, grid.xm)
    dw̄ = ddx(w̄, grid.xm)
    shear_prod = xmean(
        -uv*dv̄ - uw*dw̄, grid.xc
    )
    eps_bar = xmean(ν*(dv̄**2 + dw̄**2), grid.xc)
    epsilon = read_mean(folder, "epsilon")
    eps_p = xmean(epsilon, grid.xc) - eps_bar
    
    Tbar = read_mean(folder, "Tbar")
    qvT = read_mean(folder, "vyT")
    qT_bar = xmean(v̄*Tbar, grid.xc)
    qT_p = xmean(qvT,grid.xc) - qT_bar
    
    vr = refine(v̄, grid)
    Sbar = read_mean(folder, "Sbar")
    qvS = read_mean(folder, "vyS")
    qS_bar = xmean(vr*Sbar, grid.xcr)
    qS_p = xmean(qvS,grid.xcr) - qS_bar
    
    t = mean_time(folder)
    
    dTbar = ddx(Tbar, grid.xm, BClo=0.5, BCup=-0.5)
    qT = read_mean(folder, "vxT")
    buoy_prodT = xmean(-qT*dTbar, grid.xc)
    chiT_bar = xmean(κT*dTbar**2, grid.xc)
    chiT = read_mean(folder, "chiT")
    chiT_p = xmean(chiT, grid.xc) - chiT_bar
    qT = -0.5*κT*(dTbar[0,:] + dTbar[-1,:])
    
    dSbar = ddx(Sbar, grid.xmr, BClo=-0.5, BCup=0.5)
    qS = read_mean(folder, "vxS")
    buoy_prodS = xmean(-qS*dSbar, grid.xcr)
    chiS_bar = xmean(κS*dSbar**2, grid.xcr)
    chiS = read_mean(folder, "chiS")
    chiS_p = xmean(chiS, grid.xcr) - chiS_bar
    qS = 0.5*κS*(dSbar[0,:] + dSbar[-1,:])
    
    RaS = inputs.RayS*np.ones(len(t))
    PrS = inputs.PraS*np.ones(len(t))
    RaT = inputs.RayT*np.ones(len(t))
    PrT = inputs.PraT*np.ones(len(t))
    
    df = DataFrame({
        "t": t, "shear_prod": shear_prod,
        "eps_bar": eps_bar, "eps_p": eps_p,
        "qT_bar": qT_bar, "qT_p": qT_p,
        "qS_bar": qS_bar, "qS_p": qS_p,
        "chiT_bar": chiT_bar, "chiT_p": chiT_p,
        "chiS_bar": chiS_bar, "chiS_p": chiS_p,
        "qT": qT, "qS": qS,
        "buoy_prodT": buoy_prodT, "buoy_prodS": buoy_prodS,
        "RaT": RaT, "PrT": PrT, "RaS": RaS, "PrS": PrS
    })
    return df

def crossover_locations(A1, A2, x):
    """
    Returns the time series of the two crossover points of
    quantities A1 and A2 that are closest to the two walls
    """
    A = A1 - A2
    nt = A.shape[1]
    xlo = np.zeros(nt)
    xup = np.zeros(nt)
    for i in range(nt):
        zero_indices = np.where(np.diff(np.sign(A[:,i])))[0]
        if zero_indices.size>1:
            ix = zero_indices[0]
            xlo[i] = x[ix] + (x[ix+1] - x[ix])/(A[ix+1,i] - A[ix,i])*(0 - A[ix,i])
            ix = zero_indices[-1]
            xup[i] = x[ix] + (x[ix+1] - x[ix])/(A[ix+1,i] - A[ix,i])*(0 - A[ix,i])
        else:
            xlo[i] = np.nan
            xup[i] = np.nan
    return xlo, xup

def boundary_layers(folder):
    """
    Returns a DataFrame containing the time series of various
    definitions of boundary layer widths from each wall.
    """
    inputs = InputParams(folder)
    grid = Grid(folder)
    if inputs.FFscaleS and inputs.flagsal:
        Re = np.sqrt(inputs.RayS/inputs.PraS)
    else:
        Re = np.sqrt(inputs.RayT/inputs.PraT)
    PeT = Re*inputs.PraT
    if inputs.flagsal:
        PeS = Re*inputs.PraS
    
    # Compute velocity peak location and displacement thickness
    vybar = read_mean(folder, "vybar")
    dvdx = ddx(vybar, grid.xm)
    nt = vybar.shape[1]
    dmaxL, dmaxU = np.zeros(nt), np.zeros(nt)
    VmaxL, VmaxU = np.zeros(nt), np.zeros(nt)
    dstarL, dstarU = np.zeros(nt), np.zeros(nt)

    for i in range(nt):
        zero_indices = np.where(np.diff(np.sign(dvdx[:,i])))[0]
        if zero_indices.size > 0:
            ix = zero_indices[0]
            # Location and magnitude of velocity peak by lower wall
            dmaxL[i], VmaxL[i] = hermite_max(grid.xm[ix-1:ix+3], vybar[ix-1:ix+3,i])
            xx = np.append(grid.xm[:ix+1], dmaxL[i])
            vv = np.append(vybar[:ix+1,i], VmaxL[i])
            # Displacement thickness at lower wall
            dstarL[i] = np.trapz(vv,x=xx)/VmaxL[i]
            ix = zero_indices[-1]
            # Location and magnitude of velocity peak by upper wall
            dmaxU[i], VmaxU[i] = hermite_max(grid.xm[ix-1:ix+3], vybar[ix-1:ix+3,i])
            xx = np.insert(grid.xm[ix+1:], 0, dmaxU[i])
            vv = np.insert(vybar[ix+1:,i], 0, VmaxU[i])
            # Displacement thickness at lower wall
            dstarU[i] = np.trapz(vv,x=xx)/VmaxU[i]
        else:
            dmaxL[i], VmaxL[i], dstarL[i] = np.nan, np.nan, np.nan
            dmaxU[i], VmaxU[i], dstarU[i] = np.nan, np.nan, np.nan

    # Compute shear-based boundary layer
    # Wall shear:
    dvwL = vybar[0,:]/grid.xm[0]
    dvwU = vybar[-1,:]/(1.0 - grid.xm[-1])
    dVL, dVU = VmaxL/dvwL, VmaxU/dvwU

    # Compute dissipation crossover location
    epsilon = read_mean(folder, "epsilon")
    dwdx = ddx(read_mean(folder, "vzbar"), grid.xm)
    epsbar = (dvdx**2 + dwdx**2)/Re
    epsper = epsilon - epsbar
    dVdL, dVdU = crossover_locations(epsbar, epsper, grid.xm)

    # Temperature flux-based boundary layer
    # Set temperature boundary conditions:
    if inputs.RayT > 0:
        if inputs.flagpf or (not inputs.inslwN):
            Tup, Tlo = 0.0, 1.0
        else:
            Tup, Tlo = -0.5, 0.5
    else:
        if inputs.flagpf:
            Tup, Tlo = 1.0, 0.0
        else:
            Tup, Tlo = 0.5, -0.5
    Tbar = read_mean(folder, "Tbar")
    dTwL = (Tbar[0,:] - Tlo)/grid.xm[0]
    dTwU = (Tup - Tbar[-1,:])/(1 - grid.xm[-1])
    dTL, dTU = 1/2/np.abs(dTwL), 1/2/np.abs(dTwU)

    # Thermal flux crossover location
    dTdx = ddx(Tbar, grid.xm, BClo=Tlo, BCup=Tup)
    vxT = read_mean(folder, "vxT")
    dfTL, dfTU = crossover_locations(dTdx/PeT, -vxT, grid.xm)

    # Thermal dissipation crossover location
    chi = read_mean(folder, "chiT")
    chibar = dTdx**2/PeT
    chiper = chi - chibar
    ddTL, ddTU = crossover_locations(chibar, chiper, grid.xm)

    # Salinity boundary layers
    if inputs.flagsal:
        if inputs.flagpf:
            Sup, Slo = 0.0, 1.0
        elif inputs.RayS > 0:
            Sup, Slo = 0.5, -0.5
        else:
            Sup, Slo = -0.5, 0.5
        # Nusselt-based boundary layer
        Sbar = read_mean(folder, "Sbar")
        dSwL = (Sbar[0,:] - Slo)/grid.xmr[0]
        dSwU = (Sup - Sbar[-1,:])/(1 - grid.xmr[-1])
        dSL, dSU = 1/2/np.abs(dSwL), 1/2/np.abs(dSwU)

        # Solutal flux crossover location
        dSdx = ddx(Sbar, grid.xmr, BClo=Slo, BCup=Sup)
        vxS = read_mean(folder, "vxS")
        dfSL, dfSU = crossover_locations(dSdx/PeS, -vxS, grid.xmr)

        # Solutal dissipation crossover location
        chi = read_mean(folder, "chiS")
        chibar = dSdx**2/PeS
        chiper = chi - chibar
        ddSL, ddSU = crossover_locations(chibar, chiper, grid.xmr)
    
    t = mean_time(folder)

    df = DataFrame({
        "t": t,
        "dvL": dVL, "dvU": dVU, # Shear-based boundary layers
        "dstarL": dstarL, "dstarU": dstarU, # Displacement thickness
        "depsL": dVdL, "depsU": 1 - dVdU, # KE dissipation crossover
        "dmaxL": dmaxL, "dmaxU": 1 - dmaxU, # Velocity maximum location
        "dTL": dTL, "dTU": dTU, # Nusselt-based BL
        "dfluxTL": dfTL, "dfluxTU": 1 - dfTU, # Thermal flux crossover
        "dchiTL": ddTL, "dchiTU": 1 - ddTU # Thermal dissipation crossover
    })

    if inputs.flagsal:
        df["dSL"], df["dSU"] = dSL, dSU # Solutal Nusselt BL
        df["dfluxSL"], df["dfluxSU"] = dfSL, 1 - dfSU # Solutal flux crossover
        df["dchiSL"], df["dchiSU"] = ddSL, 1 - ddSU # Solutal dissipation crossover

    return df