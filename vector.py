
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


def _cluster_coeff(cons, simulation=None, axes=(0,1), total=None):
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




