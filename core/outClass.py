
class Output(object):
    """
    Class that holds the output of the model and processes it.
    It's diferent from the Simulation class, but it can't do stuff without it.
    """
    def __init__(self, oppath, apply_basename=False, n_cons=[1], separate_ncon=True, verbose=False):
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
        opfiles = pd.DataFrame(columns=['vel_sc', 'con_tt', 'temp_t', 'vel_t', 'div_z0_t'])
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

        print('Starting to read output as ', end='')
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
        # list div_z0_t files
        div_z0_t = sorted(glob(path.join(oppath, 'div_z0_t*out')))
        if div_z0_t: print('div_z0_t', end=' ')
        for fname in div_z0_t:
            ndtime = utils.nameParser(fname)
            opfiles.loc[ndtime, 'div_z0_t'] = fname
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
        print('... done reading and organizing.')
        return


    def compose_pcon(self, times=None, t_ini=0, t_end=None, simulation=None, as_dataarray=True,
            apply_to_z=False, z_function=lambda x: x[:,:,0], dtype=None):
        """
        Puts together particle outputs in space (for ENDLESS patches) and in time
        creating one big 5-dimensional numpy array in return, with the axes being
        time, x, y, z, n_con
        where n_con is the number of the particle

        t_ini, t_end: int
            initial and final time steps
        simulation: lespy.Simulation object
            simulation to consider when processing the results
        as_dataarray: bool
            if true, return xarray.DataArray. If false, return np.ndarray object
        apply_to_z: bool
            If true, z_function is applied to each instantaneous 3D patch before
            merging into final array and final array has shape (time, x, y[, n_con]).
            If false, full 3D patch is merged and final array has shape (time,x,y,z[,n_con]).
            Currently only tested for endless (con_tt outputs).
        z_function: function
            Function (should be manually set to apply on z axis) to be applied on z axis
            if apply_to_z in True.
        dtype: numpy.type, python.type
            type to initialize the array with.
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
        if type(times)!=type(None):
            cons = self.binaries.loc[times, labels].dropna(axis=1, how='all')
        else:
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
                for patches in cons:
                    for patch in patches[n_i]:
                        ndtime, ncon, row, col = utils.nameParser(patch)
                        rows.append(row)
                        cols.append(col)
                delta_rows = max(rows) - min(rows) + 1
                delta_cols = max(cols) - min(cols) + 1
                #---------
    
                #---------
                # For too-large sizes, it's better to integrate over z patch-by-patch
                if apply_to_z:
                    print('Creating array of ',(len(cons), delta_cols*sim.nx, delta_rows*sim.ny))
                    pcon = np.full((len(cons), delta_cols*sim.nx, delta_rows*sim.ny), np.nan, dtype=dtype)
                else:
                    print('Creating array of ',(len(cons), delta_cols*sim.nx, delta_rows*sim.ny, sim.nz_tot))
                    pcon = np.full((len(cons), delta_cols*sim.nx, delta_rows*sim.ny, sim.nz_tot), np.nan, dtype=dtype)
                #---------

                #---------
                # Now we iterate again over the time steps to read the files of one pcon size
                for i, ndtime in enumerate(cons.index):
                    print(i, ndtime, n_con)
                    for patch in cons.loc[ndtime][n_i]:
                        ndtime, ncon, row, col = utils.nameParser(patch)
                        con = routines.readBinary2(patch, simulation=sim, read_pcon=True, read_just_pcon=True)
        
                        min_yj = (row - min(rows))*sim.ny
                        max_yj = min_yj + sim.ny
    
                        min_xi = (col - min(cols))*sim.nx
                        max_xi = min_xi + sim.nx
    
                        #--------
                        # Here we apply the z_function to 3D patch, making it 2D (x, y), and merge
                        if apply_to_z:
                            pcon[i, min_xi:max_xi, min_yj:max_yj] = z_function(con[:sim.nx,:,:])
                        else:
                            pcon[i, min_xi:max_xi, min_yj:max_yj, :] = con[:sim.nx,:,:]
                        #--------
                pcons.append(pcon)
                #---------
            #--------
        #---------

        #---------
        # In this case there's no endless!
        elif cons.columns.tolist() == ['vel_sc']:
            cons = cons.vel_sc

            #---------
            # For too-large sizes, it's better to integrate over z patch-by-patch
            if apply_to_z:
                print('Creating array of ',(len(cons), sim.nx, sim.ny, sim.n_con))
                pcons = np.full((len(cons), sim.nx, sim.ny, sim.n_con), np.nan, dtype=dtype)
            else:
                print('Creating array of ',(len(cons), sim.nx, sim.ny, sim.nz_tot, sim.n_con))
                pcons = np.full((len(cons), sim.nx, sim.ny, sim.nz_tot, sim.n_con), np.nan, dtype=dtype)
            #---------

            for i, fname in enumerate(cons):
                print(i, cons.index[i], fname)
                con = routines.readBinary2(fname, simulation=sim, n_con=sim.n_con, read_just_pcon=True)

                #--------
                # Here we apply the z_function to 4D patch, making it 3D (x, y, n_con), and merge
                if apply_to_z:
                    pcons[i,:,:,:] = z_function(con[:sim.nx,:,:,:])
                else:
                    pcons[i,:,:,:,:] = con[:sim.nx,:,:,:]
                #--------
        #---------

        else:
            print(cons.columns.tolist())

        if as_dataarray:
            return utils.get_dataarray(pcons, simulation=sim, with_time=cons.index.tolist())
        else:
            return pcons


    def compose_time(self, times=None, t_ini=0, t_end=None, simulation=None, trim=True, as_dataarray=True,
            apply_to_z=False, z_function=lambda x: x[:,:,0], dtype=None):
        """
        Puts together everything in time, but can't compose endless patches. Output
        matrices indexes are
        time, x, y, z

        Parameters
        ----------
        self: lp.Output
        times: slice
            slice with which times to read.
        simulation: lp.Simulation
            to be used as base
        trim: bool
            whether to trim two extra points in x axis
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
        if type(times)!=type(None):
            bins = self.binaries.loc[times, labels].dropna(axis=1, how='all')
        else:
            if t_end is None:
                bins = self.binaries.loc[t_ini:, labels].dropna(axis=1, how='all')
            else:
                bins = self.binaries.loc[t_ini:t_end, labels].dropna(axis=1, how='all')
        #--------------

        #---------
        # Definition of output with time, x, y[ and z]
        if apply_to_z:
            print('Creating 4 arrays of {}, {}, {}...'.format(len(bins), sim.domain.ld, sim.ny), end='')
            u = np.full((len(bins), sim.domain.ld, sim.ny), np.nan)
            v = np.full((len(bins), sim.domain.ld, sim.ny), np.nan)
            w = np.full((len(bins), sim.domain.ld, sim.ny), np.nan)
            T = np.full((len(bins), sim.domain.ld, sim.ny), np.nan)
        else: 
            print('Creating 4 arrays of {}, {}, {}, {}...'.format(len(bins), sim.domain.ld, sim.ny, sim.nz_tot), end='')
            u = np.full((len(bins), sim.domain.ld, sim.ny, sim.nz_tot), np.nan)
            v = np.full((len(bins), sim.domain.ld, sim.ny, sim.nz_tot), np.nan)
            w = np.full((len(bins), sim.domain.ld, sim.ny, sim.nz_tot), np.nan)
            T = np.full((len(bins), sim.domain.ld, sim.ny, sim.nz_tot), np.nan)
        print(' done.')
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

                #------
                # Reduce z coordinate if theres a z_function
                if apply_to_z:
                    aux = [ z_function(el) for el in aux ]
                #------

                try:
                    ui,vi,wi,Ti = aux
                    u[i] = ui
                    v[i] = vi
                    w[i] = wi
                    T[i] = Ti
                except ValueError:
                    try: 
                        ui,vi,wi = aux
                        u[i] = ui
                        v[i] = vi
                        w[i] = wi
                    except ValueError:
                        Ti = aux
                        T[i] = Ti
            #---------
        #---------

        #---------
        # Trims the extra node(s) at the end of the x coordinate
        if trim:
            u = u[:, :sim.nx]
            v = v[:, :sim.nx]
            w = w[:, :sim.nx]
            T = T[:, :sim.nx]
        #---------

        #---------
        # Passes from numpy.array to xarray.DataArray, so that the coordinates go with the data
        if as_dataarray:
            u,v,w,T = utils.get_dataarray([u, v, w, T], simulation=sim, with_time=bins.index.tolist())
        #---------

        return u,v,w,T



    def compose_divs(self, times=None, t_ini=0, t_end=None, simulation=None, trim=True, as_dataarray=True):
        """
        Puts together everything in time, but can't compose endless patches. Output
        matrices indexes are
        time, x, y, z

        Parameters
        ----------
        self: lp.Output
        times: slice
            slice with which times to read.
        simulation: lp.Simulation
            to be used as base
        trim: bool
            whether to trim two extra points in x axis
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
        labels = ['div_z0_t']
        #--------------

        #--------------
        # Adjust intervals
        if type(times)!=type(None):
            bins = self.binaries.loc[times, labels].dropna(axis=1, how='all')
        else:
            if t_end is None:
                bins = self.binaries.loc[t_ini:, labels].dropna(axis=1, how='all')
            else:
                bins = self.binaries.loc[t_ini:t_end, labels].dropna(axis=1, how='all')
        #--------------

        #---------
        # Definition of output with time, x, y[ and z]
        print('Creating 4 arrays of {}, {}, {}...'.format(len(bins), sim.domain.ld, sim.ny), end='')
        du = np.full((len(bins), sim.domain.ld, sim.ny), np.nan)
        dv = np.full((len(bins), sim.domain.ld, sim.ny), np.nan)
        dw = np.full((len(bins), sim.domain.ld, sim.ny), np.nan)
        print(' done.')
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
                ui,vi,wi = aux
                du[i] = ui
                dv[i] = vi
                dw[i] = wi
            #---------
        #---------

        #---------
        # Trims the extra node(s) at the end of the x coordinate
        if trim:
            du = du[:, :sim.nx]
            dv = dv[:, :sim.nx]
            dw = dw[:, :sim.nx]
        #---------

        #---------
        # Passes from numpy.array to xarray.DataArray, so that the coordinates go with the data
        if as_dataarray:
            du, dv, dw = utils.get_dataarray([du, dv, dw], simulation=sim, with_time=bins.index.tolist())
        #---------

        return du, dv, dw


    def _compose_par(self, nprocs=10, t_ini=None, t_end=None, **kwargs):
        import multiprocessing
        from multiprocessing import Pool
        bounds = [t_ini, t_end]
        def par_compose_time(x):
            return out.compose_time(simulation=sim, t_ini=x[0], t_end=x[1], **kwargs)
        pool = Pool(processes=nprocs)
        outs = pool.map(par_compose_time, bounds)
        return zip(*outs)

