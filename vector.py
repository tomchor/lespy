
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



def correlate_2d(Vars, simulation=None):
    """
    Calculates 2D correlations in the arrays contained in Vars
    
    Vars should be five-dimensional. Dimensions are
    1: variable (different variables are allowed)
    2: time (the output will be a mean in this dimension)
    3: z (one correlation matrix for each index of this dimension)
    4: y (used in the correlation)
    5: x (used in the correlation)

    So the ouput will be 4-dimensional. The dimensions will be
    1: variable
    2: z
    3: delta_y
    4: delta_x
    """
    import numpy as np
    sim = simulation
    key0 = list(Vars.keys())[0]

    dx = sim.domain.dx
    dy = sim.domain.dy
    dz = sim.domain.dz

    nx = sim.nx
    ny = sim.ny
    nz = Vars[key0].shape[-1]

    #------
    # useful to get y and x sizes of the correlation array
    nyc = int(np.floor(ny/2)+1)
    nxc = int(np.floor(nx/2)+1)
    #------

    #------
    # Initiate dicts (should be improved later)
    vel = dict()
    vAvg = dict()
    vVar = dict()
    vCorr = dict()
    #------

    #------
    # Transpose (because first version worked with z,y,x) and calculate mean and std
    for key in Vars:
        Vars[key] = Vars[key].T
        vAvg[key] = Vars[key].mean(axis=(1,2,3))#.values
        vVar[key] = Vars[key].std(axis=(1,2,3))#.values
    #------

    for key in Vars:
        vCorr[key] = np.zeros((nz, 2*nyc-1, 2*nxc-1))

    for its in range(Vars[key0].shape[-1]):
    #for its in range(Vars[list(Vars.keys())[0]].shape[-1]):
    #for its,ts in enumerate(w.indexes['time']):
        print('Avg for:', its)
        # variance
        for key in Vars:
            #vel[key] = Vars[key].values[:,:,:,its]
            vel[key] = Vars[key][:,:,:,its]
            vVar[key] = vel[key].std(axis=(1,2))

        # correlation
        for key in Vars:
            print(key)
            for iz in range(Vars[key].shape[0]):
                print('Calculating separately for:',iz)
                for iy in range(2*nyc-1):
                    for ix in range(2*nxc-1):
                        vCorr[key][iz, iy, ix] = np.average([vCorr[key][iz, iy, ix],
                            np.average((vel[key][iz, :, :] - vAvg[key][iz]) * np.roll(
                            np.roll(vel[key][iz, :, :] - vAvg[key][iz],
                                nxc-1-ix, axis=1), nyc-1-iy, axis=0))\
                                        / np.maximum(vVar[key][iz], 10e-40) ** 2],
                                        axis=0, weights=(its, 1))
    #-----
    # From z,y,x to x,y,z
    for key in Vars:
        vCorr[key] = vCorr[key].T
    #-----

    y, x = np.mgrid[dy*(-nyc+1):dy*nyc:dy, dx*(-nxc+1):dx*nxc:dx]
    return x, y, vCorr






def _corr_field(cons, simulation=None, dr=[], axes=(0,1)):
    """Calculates the clustering strength"""
    import numpy as np
    from matplotlib import pyplot as plt
    from scipy.signal import correlate

    sim=simulation
    dx=sim.domain.dx
    dy=sim.domain.dy

    nx=cons.shape[axes[0]]
    ny=cons.shape[axes[1]]
    nx2=nx//10
    ny2=ny//10
    xrang = np.arange(-nx2, nx2+1)
    yrang = np.arange(-ny2, ny2+1)

    Dx=np.arange(-nx2, nx2+1)*dx
    Dy=np.arange(-ny2, ny2+1)*dy
    Dxv, Dyv = np.meshgrid(Dx, Dy)
    Dr = Dxv**2. + Dyv**2.
    expd = np.exp(-Dr/L**2.)

    try:
        cons = cons.values
    except:
        pass

    Cst = np.full(cons.shape, np.nan)
    for i0 in range(Cst.shape[0]):
        for j0 in range(Cst.shape[1]):
            rcon = cons.take(xrang+i0, axis=axes[0], mode='wrap').take(yrang+j0, axis=axes[1], mode='wrap')
            cexp = rcon*expd
            Cst[i0, j0] = np.trapz(np.trapz(cexp, axis=axes[1], dx=dy), axis=axes[0], dx=dx)

    return Cst


def _cluster_field(cons, simulation=None, dr=[], axes=(0,1), L=1.):
    """Calculates the clustering strength"""
    import numpy as np
    from matplotlib import pyplot as plt
    #import scipy.integrate as sp

    sim=simulation
    dx=sim.domain.dx
    dy=sim.domain.dy

    nx=cons.shape[axes[0]]
    ny=cons.shape[axes[1]]
    nx2=nx//10
    ny2=ny//10
    xrang = np.arange(-nx2, nx2+1)
    yrang = np.arange(-ny2, ny2+1)

    Dx=np.arange(-nx2, nx2+1)*dx
    Dy=np.arange(-ny2, ny2+1)*dy
    Dxv, Dyv = np.meshgrid(Dx, Dy)
    Dr = Dxv**2. + Dyv**2.
    expd = np.exp(-Dr/L**2.)

    try:
        cons = cons.values
    except:
        pass

    Cst = np.full(cons.shape, np.nan)
    for i0 in range(Cst.shape[0]):
        for j0 in range(Cst.shape[1]):
            rcon = cons.take(xrang+i0, axis=axes[0], mode='wrap').take(yrang+j0, axis=axes[1], mode='wrap')
            cexp = rcon*expd
            Cst[i0, j0] = np.trapz(np.trapz(cexp, axis=axes[1], dx=dy), axis=axes[0], dx=dx)

    return Cst

