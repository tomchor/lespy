def radial_dist(data, cond0=lambda x, y: x<=np.percentile(y,5), 
        condr=lambda x, y: x>=np.percentile(y,95), simulation=None, bins=None):
    """
    Calculates the radial distribution for data

    data: np.ndarray
        indices are:
            0: time
            1: x
            2: y
    """
    import numpy as np
    sim=simulation
    nt, Lx1, Ly1 = data.shape
    Lx, Ly = (np.array([Lx1, Ly1])/2).astype(int)
    x = np.arange(-Lx,-Lx+Lx1,1)*sim.domain.dx
    y = np.arange(-Ly,-Ly+Ly1,1)*sim.domain.dy
    xx, yy = np.meshgrid(x, y)
    r = np.sqrt(xx**2. + yy**2.)

    if type(bins)==type(None):
        bins = np.arange(0, 700, 10)

    x = np.arange(-Lx,-Lx+Lx1,1)
    y = np.arange(-Ly,-Ly+Ly1,1)

    full_hist = np.zeros((nt, Lx1, Ly1, len(bins)-1))
    for it in range(nt):
        origins = np.where(cond0(data[it], data[it]))
        for ix,iy in zip(*origins):
            rolled = np.roll(data[it], -x[ix], axis=0)
            rolled = np.roll(rolled,-y[iy],axis=1)
            high_r = r[ condr(rolled, rolled) ]
            full_hist[it,ix,iy,:] = np.histogram(high_r, bins=bins)[0]
    hist = full_hist.mean(axis=(1,2))
    summ = hist.sum(axis=(1), keepdims=True)
    summ[ summ==0 ] = 1.
    hist = hist/summ
    norm = np.histogram(r, bins=bins)[0]
    hist = hist/norm
    centers = (bins[:-1]+bins[1:])/2
    return hist, centers


def radial_homogFunction(Vars, simulation=None, nc=None, func=None):
    """ Calculates the normalized conditional density as a function of radius """
    from . import utils
    import numpy as np
    sim=simulation
    timelength=Vars.shape[1]

    if type(func)==type(None):
        func=condnorm2d_fft
    if type(nc)==type(None):
        nc=sim.nx//2

    x, y, condC = func(Vars, simulation=sim)
    nv, nt, nnx, nny = condC.shape

    print('Calculating phi(r) from phi(x,y) ... ')
    Rs = utils.radial_prof(condC[0,0], simulation=sim, axes=(0,1))[0]
    rCond = np.zeros((nv, nt, Rs.shape[0]), dtype=np.float64)
    for iv in range(nv):
        print('For variable {} of {}'.format(iv+1,nv))
        rCond[iv,:,:] = utils.radial_prof3D(condC[iv], simulation=sim)[1]
    print('done')
    return np.array(Rs), np.array(rCond)



def condnorm2d_fft(Vars, simulation=None, dx=None, dy=None):
    """
    Calculates 2D (modified) correlations in the arrays contained in Vars

    Parameters
    ----------
    Vars: np.array
        4-dimensional array with the data. Dimensions are
            0: index for variables (obviously separate calculation for each variable)
            1: time (one correlation for each time as well)
            2: x (used in the correlation)
            3: y (used in the correlation)
    simulation: lespy.Simulation
    nyc, nyx: int, int
        number of delta_ys and delta_xs to sample in each direction.
        if nyc=10, for example. The x-dim of the output matrix will be 21.
    dx, dy: float, float
        resolution. Overwriten by simulation keyword.

    Returns
    -------
    The ouput will be 4-dimensional. The dimensions will be
    0: variables
    1: time
    2: delta_y
    3: delta_x
    """
    import numpy as _np
    from numpy.fft import fft2, ifft2
    sim = simulation

    #------
    # Resolution and size
    if dx==None and dy==None:
        if sim!=None:
            dx = sim.domain.dx
            dy = sim.domain.dy
        else:
            dx, dy = 1., 1.
    nv, nt, nx, ny = Vars.shape
    #------


    #------
    # useful to get y and x sizes of the correlation array
    nxc = int(_np.floor(nx/2))
    nyc = int(_np.floor(ny/2))
    #------

    #------
    # Calculate mean, var and correlation matrix for calculations
    vAvg = Vars.mean(axis=(2,3), keepdims=True)
    Fluct = Vars - vAvg
    vVar = Vars.var(axis=(2,3))
    vCorr = _np.zeros((nv, nt, 2*nxc+1, 2*nyc+1))
    #------

    #---------
    # Calculate the correlation
    for iv in range(nv):
        print('Calculating phi(x,y) for variable {} of {}...'.format(iv+1, nv), end='')
        for it in range(nt):
            a=Vars[iv,it]
            fftv=fft2(a)
            aux = _np.roll(_np.roll(ifft2(fftv.conj()*fftv).real, nxc, axis=0), nyc, axis=1)
            vCorr[iv,it,:-1,:-1] = (aux/(nx*ny))/a.mean()**2.
            vCorr[iv,it,-1] = vCorr[iv,it,0]
            vCorr[iv,it,:,-1] = vCorr[iv,it,:,0]
        print('done')
    #---------
  
    print('calculating grid...', end='')
    y, x = _np.mgrid[-dy*nyc:dy*nyc+1:dy, -dx*nxc:dx*nxc+1:dx]
    print('done')
    return x, y, vCorr



def correlate2d_fft(Vars, simulation=None, dx=None, dy=None):
    """
    Calculates 2D correlations in the arrays contained in Vars

    Parameters
    ----------
    Vars: np.array
        4-dimensional array with the data. Dimensions are
            0: index for variables (obviously separate calculation for each variable)
            1: time (one correlation for each time as well)
            2: x (used in the correlation)
            3: y (used in the correlation)
    simulation: lespy.Simulation
    nyc, nyx: int, int
        number of delta_ys and delta_xs to sample in each direction.
        if nyc=10, for example. The x-dim of the output matrix will be 21.
    dx, dy: float, float
        resolution. Overwriten by simulation keyword.

    Returns
    -------
    The ouput will be 4-dimensional. The dimensions will be
    0: variables
    1: time
    2: delta_y
    3: delta_x
    """
    import numpy as _np
    from numpy.fft import fft2, ifft2
    sim = simulation

    #------
    # Resolution and size
    if dx==None and dy==None:
        if sim!=None:
            dx = sim.domain.dx
            dy = sim.domain.dy
        else:
            dx, dy = 1., 1.
    nv, nt, nx, ny = Vars.shape
    #------


    #------
    # useful to get y and x sizes of the correlation array
    nxc = int(_np.floor(nx/2))
    nyc = int(_np.floor(ny/2))
    #------

    #------
    # Calculate mean, var and correlation matrix for calculations
    vAvg = Vars.mean(axis=(2,3), keepdims=True)
    Fluct = Vars - vAvg
    vVar = Vars.var(axis=(2,3))
    vCorr = _np.zeros((nv, nt, 2*nxc+1, 2*nyc+1))
    #------

    #---------
    # Calculate the correlation
    for iv in range(nv):
        print('Calculating separately for variable {} of {}...'.format(iv, nv-1), end='')
        for it in range(nt):
            a=Vars[iv,it]
            fftv=fft2(a)
            aux = _np.roll(_np.roll(ifft2(fftv.conj()*fftv).real, nxc, axis=0), nyc, axis=1)
            vCorr[iv,it,:-1,:-1] = ((aux/(nx*ny))-a.mean()**2.)/a.var()
            vCorr[iv,it,-1] = vCorr[iv,it,0]
            vCorr[iv,it,:,-1] = vCorr[iv,it,:,0]
        print('done')
    #---------
  
    y, x = _np.mgrid[-dy*nyc:dy*nyc+1:dy, -dx*nxc:dx*nxc+1:dx]
    return x, y, vCorr


def moving_average(x, window, **kwargs):
    """
    Redirects to pandas.rolling_mean

    Parameters
    ----------
    axis: float 
        keyword works just as in numpy
    center: bool
        center or forward?
    min_periods: int, default None
        Minimum number of observations in window required to have a value (otherwise result is NaN).
    """
    import pandas as pd
    return pd.rolling_mean(x, window=window, **kwargs)

