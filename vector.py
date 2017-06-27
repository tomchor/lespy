
def vorticity(u,v,w, simulation=None, domain=None, axes=[1,2,3], as_dataarray=True):
    """Calculates the 3D relative vorticity"""
    import numpy as np
    diff = np.gradient
    sim = simulation

    if domain==None:
        domain=simulation.domain
    dx=domain.dx
    dy=domain.dy
    dz=domain.dz
    x,y,z = axes

    om_x = diff(w, axis=y)
    om_x-= diff(v, axis=z)
    om_y = diff(u, axis=z)
    om_y-= diff(w, axis=x)
    om_z = diff(v, axis=x)
    om_z-= diff(u, axis=y)
    
    if as_dataarray:
        om_x = sim.DataArray(om_x)
        om_y = sim.DataArray(om_y)
        om_z = sim.DataArray(om_z)

    return om_x, om_y, om_z


def div_2d(u, v, axes=(0,1)):
    """Calculate the 2d divergence of an ndarray"""
    from . import numerical as nm
    div = nm.diff_fft(u, axis=axes[0]) + nm.diff_fft(v, axis=axes[1])
    return div


def cluster_coeff(cons, simulation=None, axes=(0,1), total=None):
    """Clustering coefficient"""
    from . import numerical as nm
    import numpy as np
    from matplotlib import pyplot as plt

    sim=simulation
    dx=sim.domain.dx
    dy=sim.domain.dy

    if type(total)==type(None):
        total = np.trapz(np.trapz(cons, axis=axes[1], dx=sim.domain.dy), axis=axes[0], dx=sim.domain.dx)

    dCdx = nm.diff_fft(cons, axis=axes[0], dx=dx)
    dCdy = nm.diff_fft(cons, axis=axes[1], dx=dy)
   
    divC2 = dCdx**2. + dCdy**2.
    Ccoeff = np.trapz(np.trapz(divC2, axis=axes[1], dx=sim.domain.dy), axis=axes[0], dx=sim.domain.dx)
    return Ccoeff/total#, Ccoeff2/total


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



def detect_local_minima(arr):
    """
    Takes an array and detects the troughs using the local maximum filter.
    Returns a boolean mask of the troughs (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """
    import numpy as np
    import scipy.ndimage.filters as filters
    import scipy.ndimage.morphology as morphology
    # http://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#generate_binary_structure
    neighborhood = morphology.generate_binary_structure(len(arr.shape),2)
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.filters.html#minimum_filter
    local_min = (filters.minimum_filter(arr, footprint=neighborhood)==arr)
    # In order to isolate the peaks we must remove the background from the mask.
    background = (arr==0)
    # we must erode the background in order to subtract it from local_min, or a line will 
    # appear along the background border (artifact of the local minimum filter)
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#binary_erosion
    eroded_background = morphology.binary_erosion(background, structure=neighborhood, border_value=1)
    detected_minima = local_min - eroded_background
    return np.where(detected_minima)  


def moving_average(a, window=3, axis=0, method='forward'):
    """simple forward Moving average"""
    import numpy as np
    n=window
    ret = np.cumsum(a, dtype=float, axis=axis)
    if method=='forward':
        if axis==0:
            ret[n:] = ret[n:] - ret[:-n]
            return ret[n - 1:] / n
        elif axis==1:
            ret[:,n:] = ret[:,n:] - ret[:,:-n]
            return ret[:,n - 1:] / n
        elif axis==2:
            ret[:,:,n:] = ret[:,:,n:] - ret[:,:,:-n]
            return ret[:,:,n - 1:] / n
    else:
        raise Exception('Not implemented')



