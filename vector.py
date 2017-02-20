
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


def div_2d(u, v, axes=(1,2)):
    """Calculate the 2d divergence of an ndarray"""
    from . import numerical as nm
    div = nm.diff_fft(u, axis=1) + nm.diff_fft(v, axis=2)
    return div


def cluster_coeff(cons, simulation=None, axes=(0,1), total=None):
    """Clustering coefficient"""
    from . import numerical as nm
    import numpy as np

    sim=simulation
    dx=sim.domain.dx
    dy=sim.domain.dy

    if type(total)==type(None):
        total = np.trapz(np.trapz(cons**2., axis=axes[1], dx=sim.domain.dy), axis=axes[0], dx=sim.domain.dx)

    dCdx = nm.diff_fft(cons, axis=axes[0])
    dCdy = nm.diff_fft(cons, axis=axes[1])
   
    divC2 = dCdx**2. + dCdy**2.
    Ccoeff = np.trapz(np.trapz(divC2, axis=axes[1], dx=sim.domain.dy), axis=axes[0], dx=sim.domain.dx)
    return Ccoeff/total



def cluster_field(cons, simulation=None, dr=[], axes=(0,1), L=1.):
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

