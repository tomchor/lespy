
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
        buff += '\nEndless    : {}\n'.format(self.flag_endless)
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


class Output(object):
    """
    """
    def __init__(self, oppath, basename=False):
        """
        """
        from os import path
        from glob import glob
        import pandas as pd
        from .. import utils

        opfiles = pd.DataFrame(columns=['vel_sc', 'con_tt'])

        if path.isfile(oppath):
            oppath = path.dirname(oppath)
        elif path.isdir(oppath):
            pass
        else:
            print('No wildcards')
            raise Exception

        vel_sc = sorted(glob(path.join(oppath, 'vel_sc*out')))
        for fname in vel_sc:
            ndtime = utils.nameParser(fname)
            opfiles.loc[ndtime, 'vel_sc'] = fname


        con_tt = sorted(glob(path.join(oppath, 'con_tt*out')))
        if basename:
            con_tt = map(path.basename, con_tt)
        for fname in con_tt:
            ndtime, ncon, row, col = utils.nameParser(fname)

            if ndtime not in opfiles.index:
                opfiles.loc[ndtime, 'con_tt'] = [fname]
            else:
                if isinstance(opfiles.loc[ndtime, 'con_tt'], str):
                    opfiles.loc[ndtime, 'con_tt'] = [ opfiles.loc[ndtime, 'con_tt'], fname ]
                elif isinstance(opfiles.loc[ndtime, 'con_tt'], list):
                    opfiles.loc[ndtime, 'con_tt'].append(fname)
                else:
                    opfiles.loc[ndtime, 'con_tt'] = [fname]

        self.df = opfiles



    def compose_pcon(self, t_ini=0, t_end=None, simulation=None):
        """Puts together ENDLESS patches
        """
        from .. import utils, routines
        import numpy as np

        if simulation:
            sim=simulation
        else:
            raise ValueError

        labels = ['con_tt', 'vel_sc']

        #--------------
        # Adjust intervals
        if t_end is None:
            cons = self.df.loc[t_ini:, labels].dropna(axis=1, how='all')
        else:
            cons = self.df.loc[t_ini:t_end, labels].dropna(axis=1, how='all')
        #--------------

        if len(cons.columns) not in [1,2]:
            print('Columns are:', cons.columns)
            raise ValueError('Should be only con_tt or vel_sc left')

        elif (cons.columns.tolist() == ['con_tt']) or (cons.columns.tolist() == ['con_tt', 'vel_sc']):
            cons = cons.con_tt.dropna()
            rows = []
            cols = []
            for patches in cons:
                for patch in patches:
                    ndtime, ncon, row, col = utils.nameParser(patch)
                    rows.append(row)
                    cols.append(col)
            delta_rows = max(rows) - min(rows) + 1
            delta_cols = max(cols) - min(cols) + 1

            pcons = np.full((len(cons), delta_cols*sim.nx, delta_rows*sim.ny, sim.nz_tot, sim.n_con), np.nan)
            for i, ndtime in enumerate(cons.index):
                print(i, ndtime)
                for patch in cons.loc[ndtime]:
                    ndtime, ncon, row, col = utils.nameParser(patch)
                    con = routines.readBinary(patch, simulation=sim)
    
                    min_yj = (row - min(rows))*sim.ny
                    max_yj = min_yj + sim.ny

                    min_xi = (col - min(cols))*sim.nx
                    max_xi = min_xi + sim.nx

                    pcons[i, min_xi:max_xi, min_yj:max_yj, :, :] = con[:sim.nx,:,:,:]

        #---------
        # In this case there's no endless!
        elif cons.columns.tolist() == ['vel_sc']:
            cons = cons.vel_sc
            pcons = np.full((len(cons), sim.nx, sim.ny, sim.nz_tot, sim.n_con), np.nan)
            for i, fname in enumerate(cons):
                print(i, cons.index[i])
                u,v,w,T,con = routines.readBinary(fname, simulation=sim)

                pcons[i,:,:,:,:] = con[:sim.nx,:,:,:]
        #---------

        else:
            print(cons.columns.tolist())

        return pcons



    def put_together(self, t_ini=0, t_end=None, simulation=None):
        """Puts together ENDLESS patches
        """
        from .. import utils, routines
        import numpy as np

        if simulation:
            sim=simulation
        else:
            raise ValueError

        labels = ['con_tt', 'vel_sc']

        #--------------
        # Adjust intervals
        if t_end is None:
            cons = self.df.loc[t_ini:, labels].dropna(axis=1, how='all')
        else:
            cons = self.df.loc[t_ini:t_end, labels].dropna(axis=1, how='all')
        #--------------

        #--------------
        # Create output matrix
        if 'con_tt' in cons.columns:
            ttcons = cons['con_tt'].dropna(axis=0, how='all')
            rows = []
            cols = []
            for patches in ttcons:
                print(patches)
                for patch in patches:
                    ndtime, ncon, row, col = utils.nameParser(patch)
                    rows.append(row)
                    cols.append(col)
            delta_rows = max(rows) - min(rows) + 1
            delta_cols = max(cols) - min(cols) + 1
            pcons = np.full((len(cons), delta_cols*sim.nx, delta_rows*sim.ny, sim.nz_tot, sim.n_con), np.nan)
        else:
            pcons = np.full((len(cons), sim.nx, sim.ny, sim.nz_tot, sim.n_con), np.nan)
        #--------------

        if 'vel_sc' in cons.columns:
            sccons = cons['con_tt'].dropna(axis=0, how='all')

            for i, ndtime in enumerate(cons.index):
                if ndtime not in sccons.index:
                    continue
                print(i, ndtime)
                for patch in sccons.loc[ndtime]:
                    ndtime, ncon, row, col = utils.nameParser(patch)
                    con = routines.readBinary(patch, simulation=sim)

                    min_yj = (row - min(rows))*sim.ny
                    max_yj = min_yj + sim.ny

                    min_xi = (col - min(cols))*sim.nx
                    max_xi = min_xi + sim.nx

                    pcons[i, min_xi:max_xi, min_yj:max_yj, :, :] = con[:sim.nx,:,:,:]

        if 'con_tt' in cons.columns:
            ttcons = cons['con_tt'].dropna(axis=0, how='all')
            rows = []
            cols = []
            for patches in ttcons:
                print(patches)
                for patch in patches:
                    ndtime, ncon, row, col = utils.nameParser(patch)
                    rows.append(row)
                    cols.append(col)
            delta_rows = max(rows) - min(rows) + 1
            delta_cols = max(cols) - min(cols) + 1

            for i, ndtime in enumerate(cons.index):
                if ndtime not in tcons.index:
                    continue
                print(i, ndtime)
                for patch in ttcons.loc[ndtime]:
                    ndtime, ncon, row, col = utils.nameParser(patch)
                    con = routines.readBinary(patch, simulation=sim)

                    min_yj = (row - min(rows))*sim.ny
                    max_yj = min_yj + sim.ny

                    min_xi = (col - min(cols))*sim.nx
                    max_xi = min_xi + sim.nx

                    pcons[i, min_xi:max_xi, min_yj:max_yj, :, :] = con[:sim.nx,:,:,:]
        return pcons



