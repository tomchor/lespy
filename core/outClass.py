
class Output(object):
    """
    Class that holds the output of the model and processes it.
    It's diferent from the Simulation class, but it can't do stuff without it.
    """
    def __init__(self, oppath, apply_basename=False, n_cons=[1], separate_ncon=True):
        """
        Lists every output file from a simulation (so far only vel_sc, vel_t, temp_t and con_tt)

        oppath can be either the directory of the output, or one of the outuput binary files
        """
        from os import path
        from glob import glob
        import pandas as pd
        from .. import utils
        print('Starting to read output as ', end='')

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
        if vel_sc: print('vel_sc', end=' ')
        for fname in vel_sc:
            ndtime = utils.nameParser(fname)
            opfiles.loc[ndtime, 'vel_sc'] = fname
        #----------

        #----------
        # list temp_t files
        temp_t = sorted(glob(path.join(oppath, 'temp_t*out')))
        if temp_t: print('temp_t', end=' ')
        for fname in temp_t:
            ndtime = utils.nameParser(fname)
            opfiles.loc[ndtime, 'temp_t'] = fname
        #----------

        #----------
        # list vel_t files
        vel_t = sorted(glob(path.join(oppath, 'vel_t*out')))
        if vel_t: print('vel_t', end=' ')
        for fname in vel_t:
            ndtime = utils.nameParser(fname)
            opfiles.loc[ndtime, 'vel_t'] = fname
        #----------

        #----------
        # list con_tt files (each entry is a list, since there can be lots for each timestep
        con_tt = sorted(glob(path.join(oppath, 'con_tt*out')))
        if con_tt: print('con_tt', end=' ')
        if apply_basename:
            con_tt = map(path.basename, con_tt)

        for fname in con_tt:
            ndtime, ncon, row, col = utils.nameParser(fname)
            try:
                opfiles.loc[ndtime, 'con_tt'].append(fname)
            except (KeyError, AttributeError):
                opfiles.loc[ndtime, 'con_tt'] = [fname]
        print()
        #----------

        #---------
        # Separate each list entry of con_tt into list of lists (one for each n_con)
        if separate_ncon:
            nconlist = [ 's{:03d}'.format(el) for el in n_cons ]
            self.n_cons = nconlist
            for i,entries in opfiles.loc[:, 'con_tt'].dropna().iteritems():
                flist = [ [ entry for entry in entries if ncon in entry ] for ncon in nconlist ]
                opfiles.loc[i, 'con_tt'] = flist
        #---------

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

            #--------
            # We iterate over the different sizes of pcon.
            # Here n_i is a local integer and n_cons is a lit with the files for one size
            pcons = []
            for n_i, n_con in enumerate(self.n_cons):
                print(n_i, n_con)
                #---------
                # First now for each pcon size, we find the max and min row 
                # and col to find the dimension of the array that will hold the results.
                rows = []
                cols = []
                #for patches in n_cons:
                for patches in cons:
                    for patch in patches[n_i]:
                        ndtime, ncon, row, col = utils.nameParser(patch)
                        rows.append(row)
                        cols.append(col)
                delta_rows = max(rows) - min(rows) + 1
                delta_cols = max(cols) - min(cols) + 1
                #---------
    
                #---------
                # Now we iterate again over the time steps to read the files of one pcon size
                pcon = np.full((len(cons), delta_cols*sim.nx, delta_rows*sim.ny, sim.nz_tot), np.nan)
                for i, ndtime in enumerate(cons.index):
                    print(i, ndtime, n_con)
                    for patch in cons.loc[ndtime][n_i]:
                        ndtime, ncon, row, col = utils.nameParser(patch)
                        con = routines.readBinary2(patch, simulation=sim, read_pcon=True, read_just_pcon=True)
        
                        min_yj = (row - min(rows))*sim.ny
                        max_yj = min_yj + sim.ny
    
                        min_xi = (col - min(cols))*sim.nx
                        max_xi = min_xi + sim.nx
    
                        pcon[i, min_xi:max_xi, min_yj:max_yj, :] = con[:sim.nx,:,:]
                pcons.append(pcon)
                #---------
            #--------
        #---------

        #---------
        # In this case there's no endless!
        elif cons.columns.tolist() == ['vel_sc']:
            cons = cons.vel_sc
            pcons = np.full((len(cons), sim.nx, sim.ny, sim.nz_tot, sim.n_con), np.nan)
            for i, fname in enumerate(cons):
                print(i, cons.index[i], fname)
                con = routines.readBinary2(fname, simulation=sim, n_con=sim.n_con, read_just_pcon=True)

                pcons[i,:,:,:,:] = con[:sim.nx,:,:,:]
        #---------

        else:
            print(cons.columns.tolist())

        if as_dataarray:
            #----------
            # If there is more than one droplet, each can have a different domain (endless)
            if isinstance(pcons, list):
                dims = ['time', 'x', 'y', 'z']
                pcons_da = []
                for pcon in pcons:
                    x,y,z = sim.domain.makeAxes(pcon[0])
                    coords = {'time':cons.index.tolist(),
                        'x':x, 'y':y, 'z':z}
                    pcons_da.append(sim.DataArray(pcon, dims=dims, coords=coords))
                del pcons
                return pcons_da
            #----------

            #----------
            # Else, we do it just for one
            else:
                dims = ['time', 'x', 'y', 'z', 'size']
                coords = {'time':cons.index.tolist(),
                    'x':sim.domain.x, 'y':sim.domain.y, 'z':sim.domain.z,
                    'size':np.arange(pcons.shape[-1])}
                return sim.DataArray(pcons, dims=dims, coords=coords)
            #----------
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
                aux = routines.readBinary2(col, simulation=sim, read_pcon=False)
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



