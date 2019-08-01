import numpy as _np
import xarray as _xr

def vorticity(u,v,w, simulation=None, domain=None, axes=[1,2,3], as_dataarray=True):
    """Calculates the 3D relative vorticity"""
    diff = _np.gradient
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

def velgrad_tensor2d(u_vec, simulation=None, trim=True, as_DA=True, real=True):
    """
    u_vec is a list with [u, v]
    u, v should be xarrays
    """
    from . import numerical
    sim=simulation
    vec_st = ['u', 'v']

    R=_np.full((2,2), _np.nan, dtype=_np.object)
    for i, (comp, comp_st) in enumerate(zip(u_vec, vec_st)):
        for j, dim_st in enumerate("xy"):
            print('R[{},{}] = {}_{}'.format(i, j, comp_st, dim_st))
            R[i,j] = numerical.diff_fft_xr(comp, dim=dim_st, real=real)
            R[i,j].name = "d{}/d{}".format(comp_st, dim_st)

    if as_DA:
        import xarray as xr
        Dir = xr.DataArray(["x", "y"], dims=["dir"])
        Dir.name="dir"
        Comp = xr.DataArray(["u", "v"], dims=["comp"])
        Comp.name="comp"
        R = xr.concat([ xr.concat(col, dim=Dir) for col in R ], dim=Comp )

    return R

def velgrad_tensor(u_vec, simulation=None, trim=True, w_int=None, as_DA=True):
    """
    u_vec is a list with [u, v, w]
    u, v, and w should be xarrays

    trim:
        w is in different nodes, so the w components are longer. This gets rids of that
    w_int:
        w interpolated to u-nodes
    """
    from . import numerical
    sim=simulation
    vec_st = ['u', 'v', 'w']

    R=_np.full((3,3), _np.nan, dtype=_np.object)
    for i, (comp, comp_st) in enumerate(zip(u_vec, vec_st)):
        for j, dim_st in enumerate("xy"):
            print('R[{},{}] = {}_{}'.format(i, j, comp_st, dim_st))
            if comp_st=='w' and type(w_int)!=type(None):
                R[i,j] = numerical.diff_fft_xr(w_int, dim=dim_st)
            else:
                R[i,j] = numerical.diff_fft_xr(comp, dim=dim_st)
            R[i,j].name = "d{}/d{}".format(comp_st, dim_st)

        print('R[{},2] = {}_z'.format(i, comp_st))
        if comp_st=="w":
            R[i,2] = -comp.diff(dim='z')/sim.domain.dz
            R[i,2].coords['z'] = R[i,2].coords['z'] + sim.domain.dz/2
        else:
            R[i,2] = gradient_xr(comp, "z")
        R[i,2].name = "d{}/dz".format(comp_st)

    if trim:
        zlen = len(R[2,2].coords['z']) # The shortest
        for i in range(3):
            for j in range(3):
                R[i,j] = R[i,j].isel(z=slice(None,zlen-1))

    if as_DA:
        import xarray as xr
        Dir = xr.DataArray(["x", "y", "z"], dims=["dir"])
        Dir.name="dir"
        Comp = xr.DataArray(["u", "v", "w"], dims=["comp"])
        Comp.name="comp"
        R = xr.concat([ xr.concat(col, dim=Dir) for col in R ], dim=Comp )

    return R


def gradient_xr(da, dim):
    ds = da.coords[dim].diff(dim)[0].item()
    grad_da = da.copy(deep=True)
    grad_da = grad_da.where(False, _np.gradient(da, ds, axis=da.dims.index(dim)))
    return grad_da



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




def gaussian_conv(da, delta=1, dims=["x", "y"], truncate=4, how="manual", 
                  full_output=False, **kwargs):
    """ 
    Applies convolution with gaussian kernel 

    delta: is the std of the guassian (in meters)
    how: if manual we calculate kernel manually using Pope's Table 13.2. If auto, we do it using scipy's gaussian_filter1d
    """
    from scipy.ndimage import convolve1d
    from scipy.ndimage.filters import gaussian_filter1d

    da = da.copy(deep=True)
    for dim in dims:
        axis = da.dims.index(dim)
        if dim=="z":
            bc="constant"
        else:
            bc="wrap"

        #----
        # Calculate kernel resolution
        s = da.coords[dim]
        ds = abs(s.diff(dim)[0].item())
        #----
        if how=="manual":
            r = _np.arange(0, delta*truncate+ds, ds)
            r = _np.sort(_np.concatenate([-r[1:], r]))

            # calculate 1d kernel
            G = _np.sqrt(6/(_np.pi*delta**2)) * _np.exp(-6*(r/delta)**2)
            G /= G.sum()
            if delta==0: G=_np.array([1,])

            da = _xr.apply_ufunc(convolve1d, da, G, kwargs=dict(axis=axis, mode=bc, cval=0), **kwargs)

        elif how=="auto":
            da = _xr.apply_ufunc(gaussian_filter1d, da, kwargs=dict(axis=axis, mode=bc, cval=0, sigma=delta/ds, truncate=truncate), **kwargs)

    if full_output and how=="manual":
        return r, G, da
    else:
        return da




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




