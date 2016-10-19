#from decorators import autoassign as _auto
from .decorators import autoargs as _auto

class Simulation(object):
    """class for simulation parameters"""
    #@_auto
    def __init__(self, domain=None, dt=None, timelength=None, u_scale=None, inversion_depth=None, **kwargs):
        self.domain = domain
        self.dt = dt
        self.timelength = timelength
        self.u_scale = u_scale
        self.inversion_depth = inversion_depth
        self.__dict__.update(kwargs)

    def __str__(self):
        buff=' Simulation Parameters\n'+ '-'*22
        return buff

    __repr__ = __str__


def simulation(namelist, tlength_from_ke=True):
    from .utils import paramParser
    from .dmClass import domain as Dom
    params = paramParser(namelist)
    dmn = Dom(nx=params['nx'], ny=params['ny'], nz_tot=params['nz_tot'],
            nz=params['nz_tot']-1,
            lx=params['lx'], ly=params['ly'], lz=params['lz_tot'])
    
    if tlength_from_ke:
        from os import path
        import numpy as _np
        kefile = path.join(path.dirname(namelist), '../check_ke.out')
        print('opening',kefile)
        kearray = _np.loadtxt(kefile)
        print(kearray[:,0])
        kelast = int(kearray[-1,0]+1)
        kecount = len(kearray[:,0])
        if kelast == kecount:
            tlength = kelast
        else:
            print('Warning: linescount in ke_check.out is different from index label in the file.')
            print('Setting timelength to the line count.')
            tlength = kecount
    else:
        print('Warining: getting length solely from param.nml, which can be flawed.')
        tlength = params['nsteps']

    out = Simulation(domain=dmn,
            timelength=tlength,
            avglength=params['p_count'],
            dt=params['dt'],
            u_scale=params['u_star'],
            inversion_depth=params['z_i'],
            T_scale=params['t_scale'])

    return out
