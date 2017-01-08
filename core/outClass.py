
class Output(object):
    """
    Class that holds the output of the model and processes it.
    It's diferent from the Simulation class, but it can't do stuff without it.
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
        #----------

        self.binaries = opfiles
        return


    def compose_pcon(self, t_ini=0, t_end=None, simulation=None, as_dataarray=True):
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

        if as_dataarray:
            print(pcons.shape)
            dims = ['time', 'x', 'y', 'z', 'size']
            coords = {'time':cons.index.tolist(),
                    'x':sim.domain.x, 'y':sim.domain.y, 'z':sim.domain.z,
                    'size':np.arange(pcons.shape[-1])}
            return sim.DataArray(pcons, dims=dims, coords=coords)
        else:
            return pcons

    def compose_time(self, t_ini=0, t_end=None, simulation=None, trim=True, as_dataarray=True):
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
        # Iterates on the timestamps to put the time-indexed array together
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

        #---------
        # Trims the extra node(s) at the end of the x coordinate
        if trim:
            u = u[:, :sim.nx, :, :]
            v = v[:, :sim.nx, :, :]
            w = w[:, :sim.nx, :, :]
            T = T[:, :sim.nx, :, :]
        #---------

        #---------
        # Passes from numpy.array to xarray.DataArray, so that the coordinates go with the data
        if as_dataarray:
            timestamps=bins.index.tolist()
            u = sim.DataArray(u, timestamps=timestamps)
            v = sim.DataArray(v, timestamps=timestamps)
            w = sim.DataArray(w, timestamps=timestamps)
            T = sim.DataArray(T, timestamps=timestamps)
        #---------

        return u,v,w,T



