
class Simulation(object):
    """class for simulation parameters"""
    def __init__(self, from_file=None, domain=None, timelength=None, u_scale=None, inversion_depth=None, **kwargs):
        #------------
        # Start simulation from param.nml file
        if from_file!=None:
            if isinstance(from_file, str):
                aux = sim_from_file(from_file)
                self.__dict__.update(vars(aux))
                return
            else:
                raise
        #------------
        #------------
        # Start simulation from call arguments
        else:
            self.domain = domain
            self.timelength = timelength
            self.u_scale = u_scale
            self.inversion_depth = inversion_depth
            self.__dict__.update(kwargs)
        #------------


    def check(self, full=True):
        """
        Check important characteristics of the simulation
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
        if full:
            print('Coriolis timescale           : {:1.1e} timesteps'.format(int(1./self.freq_coriolis/self.dt)))

    def __str__(self):
        buff='Simulation Parameters\n'+ '-'*21
        buff += 'Endless    : {}\n'.format(self.flag_endless)
        buff += 'dt:        : {} s'.format(self.dt)
        buff+= self.domain.__str__()
        return buff

    __repr__ = __str__



def sim_from_file(namelist, tlength_from_ke=True, check_ke_file=None):
    """
    Reads and parses namelist and then calls Simulation class
    """
    from ..utils import paramParser
    from .dmClass import domain as Dom
    params = paramParser(namelist)
    dmn = Dom(nx=params['nx'], ny=params['ny'], nz_tot=params['nz_tot'],
            nz=params['nz_tot']-1,
            lx=params['lx'], ly=params['ly'], lz=params['lz_tot'])
    
    if tlength_from_ke:
        from os import path
        import pandas as _pd

        nml_dir = path.dirname(namelist)
        if check_ke_file:
            kefile = check_ke_file
        else:
            atmpt = path.join(nml_dir, 'check_ke.out')
            if path.isfile(atmpt):
                kefile = atmpt
            else:
                atmpt = path.join(nml_dir, '../check_ke.out')
                if path.isfile(atmpt):
                    kefile = atmpt
                else:
                    atmpt = path.join(nml_dir, 'output/check_ke.out')
                    if path.isfile(atmpt):
                        kefile = atmpt
                    else:
                        raise
        
        print('opening', kefile)
        try:
            kearray = _pd.read_csv(kefile, delim_whitespace=True, squeeze=True, index_col=0, header=None)
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
        print('Warning: getting length solely from param.nml, which can be flawed.')
        tlength = params['nsteps']

    out = Simulation(domain=dmn,
            timelength=tlength,
            avglength=params['p_count'],
            u_scale=params['u_star'],
            inversion_depth=params['z_i'],
            T_scale=params['t_scale'],
            **params)

    return out
