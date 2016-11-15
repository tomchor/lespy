from .decorators import autoargs as _auto

class Simulation(object):
    """class for simulation parameters"""
    def __init__(self, domain=None, dt=None, timelength=None, u_scale=None, inversion_depth=None, **kwargs):
        self.domain = domain
        self.dt = dt
        self.timelength = timelength
        self.u_scale = u_scale
        self.inversion_depth = inversion_depth
        self.__dict__.update(kwargs)

    def check(self, nproc=None):
        """Check characteristics of the simulation
        """
        CFL_x=self.u_scale*self.dt/self.domain.dx
        print('CFL (u_scale*dt/dx)          : {:.2e}'.format(CFL_x))
        print('dx/dz                        : {:2.1f}\t\t{}'.format(self.domain.dx/self.domain.dz,'-- Should be < 5 in practice'))
        print('Lx/z_inv                     : {:2.1f}\t\t{}'.format(self.domain.Lx/self.inversion_depth,'-- Should be > 6. At *least* 4.'))
        divs = []
        for i in range(2,140):
            if self.domain.Nz%i == 0:
                divs.append(i)
        print('Nz = {:03d} and is divisible by : {}'.format(self.domain.Nz, divs))


    def __str__(self):
        buff='Simulation Parameters\n'+ '-'*21
        buff+='\ndt: {} s'.format(self.dt)
        buff+= self.domain.__str__()
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
        #import numpy as _np
        import pandas as _pd
        kefile = path.join(path.dirname(namelist), '../check_ke.out')
        kefile2 = path.join(path.dirname(namelist), 'output/check_ke.out')
        print('opening',kefile)
        try:
            #kearray = _np.loadtxt(kefile)
            kearray = _pd.read_csv(kefile, delim_whitespace=True, squeeze=True, index_col=0, header=None)
            #kelast = int(kearray[-1,0]+1)
            #kecount = len(kearray[:,0])
            kelast = kearray.index[-1]
            kecount = len(kearray.index)
            if kelast == kecount:
                tlength = kelast
            else:
                print('Warning: linescount in ke_check.out is different from index label in the file.')
                print('Setting timelength to the line count.')
                tlength = kecount
        except FileNotFoundError:
            print("Coulndn't open check_ke.out.")
            tlength = params['nsteps']
    else:
        print('Warining: getting length solely from param.nml, which can be flawed.')
        tlength = params['nsteps']

    out = Simulation(domain=dmn,
            timelength=tlength,
            avglength=params['p_count'],
            u_scale=params['u_star'],
            inversion_depth=params['z_i'],
            T_scale=params['t_scale'],
            **params)
#            pcon_flag=params['pcon_flag'],
#            s_flag=params['s_flag'],
#            ix_src=params['ix_src'],
#            iy_src=params['iy_src'],
#            iz_src=params['iz_src'],
#            vel_settling=params['vel_settling'],
#            stokes_flag=['stokes_flag'])

    return out
