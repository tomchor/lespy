
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



def integrate(arr, axis=0, dx=1., chunks=1):
    """
    Re-implementation of numpy.trapz but saving RAM memory

    axis can't be -1, as opposed to numpy.trapz
    """
    arr = arr.values
    if chunks==1:
        return dx*(0.5*(arr.take(0, axis=axis) + arr.take(-1, axis=axis)) + arr.take(range(1,arr.shape[axis]-1), axis=axis).sum(axis=axis))
    else:
        import numpy as np
        newshape = [ el for index, el in enumerate(arr.shape) if index != axis ]
        integ = np.zeros(newshape)
        limits = np.linspace(0, arr.shape[axis], chunks+1, dtype=int)
        for ini, end in zip(limits,limits[1:]):
            print('Integrating axis {} from {} to {}'.format(axis,ini,end))
            arr2 = arr.take(np.arange(ini, end), axis=axis, mode='clip')
            integ += dx*(0.5*(arr2.take(0, axis=axis) + arr2.take(-1, axis=axis)) + arr2.take(range(1,arr2.shape[axis]-1), axis=axis).sum(axis=axis))
        del arr2
        return integ



def vorticity_readable(u,v,w, simulation=None, domain=None, axes=[1,2,3]):
    """Calculates the 3D relative vorticity vector (first readable version, but uses too much ram probably)"""
    import numpy as np
    diff = np.gradient

    if domain==None:
        domain=simulation.domain
    dx=domain.dx
    dy=domain.dy
    dz=domain.dz
    x,y,z = axes
    
    dwdy = diff(w, axis=y)
    dvdz = diff(v, axis=z)
    dudz = diff(u, axis=z)
    dwdx = diff(w, axis=x)
    dvdx = diff(v, axis=x)
    dudy = diff(u, axis=y)

    om_x = dwdy - dvdz
    om_y = dudz - dwdx
    om_z = dvdx - dudy

    return om_x, om_y, om_z
