
def vorticity(u,v,w, simulation=None, domain=None, axes=[1,2,3]):
    """Calculates the 3D relative vorticity"""
    import numpy as np
    diff = np.gradient

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

    return om_x, om_y, om_z


def vorticity_readable(u,v,w, simulation=None, domain=None, axes=[1,2,3]):
    """Calculates the 3D relative vorticity vector"""
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
