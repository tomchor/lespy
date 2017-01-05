
class Output(object):
    """
    """
    def __init__(self, oppath, apply_basename=False):
        """
        Lists every output file from a simulation (so far only vel_sc, vel_t, temp_t and con_tt)

        oppath can be either the directory of the output, or one of the outuput binary files
        """
        from os import path
        from glob import glob
        import pandas as pd
        from .. import utils

        #----------
        # We use pandas to organize all the files (very handy)
        opfiles = pd.DataFrame(columns=['vel_sc', 'con_tt', 'temp_t', 'vel_t'])
        #----------

        #----------
        # If we get a file, we search in the directory of that file
        if path.isfile(oppath):
            oppath = path.dirname(oppath)
        elif path.isdir(oppath):
            pass
        else:
            print('No wildcards')
            raise Exception
        #----------

        #----------
        # List vel_sc files
        vel_sc = sorted(glob(path.join(oppath, 'vel_sc*out')))
        for fname in vel_sc:
            ndtime = utils.nameParser(fname)
            opfiles.loc[ndtime, 'vel_sc'] = fname
        #----------

        #----------
        # list temp_t files
        temp_t = sorted(glob(path.join(oppath, 'temp_t*out')))
        for fname in temp_t:
            ndtime = utils.nameParser(fname)
            opfiles.loc[ndtime, 'temp_t'] = fname
        #----------

        #----------
        # list vel_t files
        vel_t = sorted(glob(path.join(oppath, 'vel_t*out')))
        for fname in vel_t:
            ndtime = utils.nameParser(fname)
            opfiles.loc[ndtime, 'vel_t'] = fname
        #----------

        #----------
        # list con_tt files (each antry is a list, since there can be lots for each timestep
        con_tt = sorted(glob(path.join(oppath, 'con_tt*out')))
        if apply_basename:
            con_tt = map(path.basename, con_tt)

        for fname in con_tt:
            ndtime, ncon, row, col = utils.nameParser(fname)
            try:
                opfiles.loc[ndtime, 'con_tt'].append(fname)
            except (KeyError, AttributeError):
                opfiles.loc[ndtime, 'con_tt'] = [fname]

#            if ndtime not in opfiles.index:
#                opfiles.loc[ndtime, 'con_tt'] = [fname]
                #if isinstance(opfiles.loc[ndtime, 'con_tt'], str):
                #    opfiles.loc[ndtime, 'con_tt'] = [ opfiles.loc[ndtime, 'con_tt'], fname ]
#            else:
#                if isinstance(opfiles.loc[ndtime, 'con_tt'], str):
#                    opfiles.loc[ndtime, 'con_tt'] = [ opfiles.loc[ndtime, 'con_tt'], fname ]
#                elif isinstance(opfiles.loc[ndtime, 'con_tt'], list):
#                    opfiles.loc[ndtime, 'con_tt'].append(fname)
#                else:
#                    opfiles.loc[ndtime, 'con_tt'] = [fname]
        #----------

        self.binaries = opfiles
        return


    def compose_pcon(self, t_ini=0, t_end=None, simulation=None):
        """
        Puts together particle outputs in space (for ENDLESS patches) and in time
        creating one big 5-dimensional numpy array in return, with the axes being
        time, x, y, z, n_con
        where n_con is the number of the particle
        """
        from .. import utils, routines
        import numpy as np

        if simulation:
            sim=simulation
        else:
            raise ValueError("You need to provide a simulation object here")

        #--------------
        # Only these types of file have pcon output (maybe con_t also)
        labels = ['con_tt', 'vel_sc']
        #--------------

        #--------------
        # Adjust intervals
        if t_end is None:
            cons = self.binaries.loc[t_ini:, labels].dropna(axis=1, how='all')
        else:
            cons = self.binaries.loc[t_ini:t_end, labels].dropna(axis=1, how='all')
        #--------------

        #--------------
        # Making sure we have the right result
        if len(cons.columns) not in [1,2]:
            raise ValueError('Columns are: {} and there should be only con_tt and/or vel_sc left'.format(cons.columns.tolist()))
        #--------------

        #--------------
        # Case with endless (with endless, all pcon should be in con_tt)
        elif 'con_tt' in cons.columns.tolist():
            cons = cons.con_tt.dropna()

            #--------------
            # We find the max and min row and col to find the dimension of matrix
            rows = []
            cols = []
            for patches in cons:
                for patch in patches:
                    ndtime, ncon, row, col = utils.nameParser(patch)
                    rows.append(row)
                    cols.append(col)
            delta_rows = max(rows) - min(rows) + 1
            delta_cols = max(cols) - min(cols) + 1
            #--------------

            pcons = np.full((len(cons), delta_cols*sim.nx, delta_rows*sim.ny, sim.nz_tot, sim.n_con), np.nan)
            for i, ndtime in enumerate(cons.index):
                print(i, ndtime)
                for patch in cons.loc[ndtime]:
                    ndtime, ncon, row, col = utils.nameParser(patch)
                    con = routines.readBinary2(patch, simulation=sim)
    
                    min_yj = (row - min(rows))*sim.ny
                    max_yj = min_yj + sim.ny

                    min_xi = (col - min(cols))*sim.nx
                    max_xi = min_xi + sim.nx

                    pcons[i, min_xi:max_xi, min_yj:max_yj, :, :] = con[:sim.nx,:,:,:]
        #---------

        #---------
        # In this case there's no endless!
        elif cons.columns.tolist() == ['vel_sc']:
            cons = cons.vel_sc
            pcons = np.full((len(cons), sim.nx, sim.ny, sim.nz_tot, sim.n_con), np.nan)
            for i, fname in enumerate(cons):
                print(i, cons.index[i])
                u,v,w,T,con = routines.readBinary2(fname, simulation=sim)

                pcons[i,:,:,:,:] = con[:sim.nx,:,:,:]
        #---------

        else:
            print(cons.columns.tolist())

        return pcons

    def compose_time(self, t_ini=0, t_end=None, simulation=None, trim=True):
        """
        Puts together everything in time, but can't compose endless patches. Output
        matrices indexes are
        time, x, y, z
        """
        from .. import utils, routines
        import numpy as np
        import pandas as pd

        if simulation:
            sim=simulation
        else:
            raise ValueError("You need to provide a simulation object here")

        #--------------
        # Only these types of file we deal with here
        labels = ['vel_sc', 'temp_t', 'vel_t']
        #--------------

        #--------------
        # Adjust intervals
        if t_end is None:
            bins = self.binaries.loc[t_ini:, labels].dropna(axis=1, how='all')
        else:
            bins = self.binaries.loc[t_ini:t_end, labels].dropna(axis=1, how='all')
        #--------------

        #---------
        # Definition of output with time, x, y and z
        u = np.full((len(bins), sim.domain.ld, sim.ny, sim.nz_tot), np.nan)
        v = np.full((len(bins), sim.domain.ld, sim.ny, sim.nz_tot), np.nan)
        w = np.full((len(bins), sim.domain.ld, sim.ny, sim.nz_tot), np.nan)
        T = np.full((len(bins), sim.domain.ld, sim.ny, sim.nz_tot), np.nan)
        #---------

        #---------
        for i, tstep in enumerate(bins.index):
            iSeries = bins.loc[tstep]

            #---------
            # Iterate between vel_sc, vel_t and temp_t
            # Currently only works with vel_sc (I think)
            for col in iSeries:
                if not isinstance(col, str): continue
                print(col)
                aux = routines.readBinary2(col, simulation=sim)
                try:
                    ui,vi,wi,Ti = aux
                    u[i,:,:,:] = ui
                    v[i,:,:,:] = vi
                    w[i,:,:,:] = wi
                    T[i,:,:,:] = Ti

                except ValueError:
                    try: 
                        ui,vi,wi = aux
                        u[i,:,:,:] = ui
                        v[i,:,:,:] = vi
                        w[i,:,:,:] = wi
                    except ValueError:
                        Ti = aux
                        T[i,:,:,:] = Ti
            #---------
        #---------

        if trim:
            u = u[:, :sim.nx, :, :]
            v = v[:, :sim.nx, :, :]
            w = w[:, :sim.nx, :, :]
            T = T[:, :sim.nx, :, :]
        return u,v,w,T


    def put_together_old(self, t_ini=0, t_end=None, simulation=None):
        """
        Puts together ENDLESS patches
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
            cons = self.binaries.loc[t_ini:, labels].dropna(axis=1, how='all')
        else:
            cons = self.binaries.loc[t_ini:t_end, labels].dropna(axis=1, how='all')
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



