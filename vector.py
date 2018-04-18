import numpy as _np

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

def velgrad_tensor(u_vec, simulation=None, trim=True, w_int=None):
    """
    u_vec is a list with [u, v, w]
    u, v, and w should be xarrays
    """
    from . import numerical
    sim=simulation
    u_st = ['u', 'v', 'w']
    dims = ['x', 'y', 'z']

    R=_np.full((3,3), _np.nan, dtype=_np.object)
    for i, (vec, vec_s) in enumerate(zip(u_vec, u_st)):
        for j, dim in enumerate(dims[:2]):
            print('R[{},{}] = {}_{}'.format(i, j, vec_s, dim))
            if type(w_int)!=type(None) and vec_s=='w':
                R[i,j] = numerical.diff_fft_xr(w_int, dim=dim)
            else:
                R[i,j] = numerical.diff_fft_xr(vec, dim=dim)

        print('R[{},2] = {}_z'.format(i, vec_s))
        R[i,2] = -vec.diff(dim='z')/sim.domain.dz
        R[i,2].coords['z'] = R[i,2].coords['z'] + sim.domain.dz/2

    if trim:
        zlen = len(R[2,2].coords['z']) # The shortest
        for i in range(3):
            for j in range(3):
                R[i,j] = R[i,j].isel(z=slice(None,zlen-1))
    return R




def get_QD(uv, verbose=True):
    """
    uv is a list of [u, v]
    u, v are xarray.DataArrays
    """
    from . import numerical
    u, v, = uv
    if verbose: print("Start differentiating ...")
    dxu = numerical.diff_fft_xr(u, dim="x")
    dyu = numerical.diff_fft_xr(u, dim="y")
    dxv = numerical.diff_fft_xr(v, dim="x")
    dyv = numerical.diff_fft_xr(v, dim="y")
    if verbose: print('end differentiating')

    hdiv = dxu + dyv
    asym = dxu**2 - 2*dxu*dyv + 4*dxv*dyu + dyv**2
    return asym, hdiv




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




