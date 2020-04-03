
class Output(object):
    """
    Class that holds the output of the model and processes it.
    It's diferent from the Simulation class, but it can't do stuff without it.
    """
    def __init__(self, oppath, apply_basename=False, n_cons=[1], separate_ncon=True, verbose=False):
        """
        Lists every output file from a simulation (so far only uvw_jt, theta_jt and pcon_jt)

        oppath can be either the directory of the output, or one of the outuput binary files
        """
        from os import path
        from glob import glob
        import pandas as pd
        from .. import utils

        #----------
        # We use pandas to organize all the files (very handy)
        opfiles = pd.DataFrame(columns=['uvw_jt', 'theta_jt', 'pcon_jt', 'uv0_jt'])
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
        self.oppath = oppath
        #----------

        print('Starting to read output as ', end='')
        #----------
        # List uvw_jt files
        uvw_jt = sorted(glob(path.join(oppath, 'uvw_jt*bin')))
        if uvw_jt: print('uvw_jt', end=' ')
        for fname in uvw_jt:
            ndtime = utils.nameParser(fname)
            opfiles.loc[ndtime, 'uvw_jt'] = fname
        #----------

        #----------
        # list theta_jt files
        theta_jt = sorted(glob(path.join(oppath, 'theta_jt*bin')))
        if theta_jt: print('theta_jt', end=' ')
        for fname in theta_jt:
            ndtime = utils.nameParser(fname)
            opfiles.loc[ndtime, 'theta_jt'] = fname
        #----------


        #----------
        # list uv0_jt files
        uv0_jt = sorted(glob(path.join(oppath, 'uv0_jt*bin')))
        if uv0_jt: print('uv0_jt', end=' ')
        for fname in uv0_jt:
            ndtime = utils.nameParser(fname)
            opfiles.loc[ndtime, 'uv0_jt'] = fname
        #----------


        #----------
        # list pcon_jt files (each entry is a list, since there can be lots for each timestep
        pcon_jt = sorted(glob(path.join(oppath, 'pcon_jt*bin')))
        if pcon_jt: print('pcon_jt', end=' ')
        if apply_basename:
            pcon_jt = map(path.basename, pcon_jt)

        for fname in pcon_jt:
            ndtime = utils.nameParser(fname)
            opfiles.loc[ndtime, 'pcon_jt'] = fname
        #----------

        self.binaries = opfiles.sort_index()
        print('... done reading and organizing.')
        return


    def compose_pcon(self, times=None, t_ini=0, t_end=None, simulation=None,
                     apply_to_z=False, z_function=lambda x: x[:,:,0], 
                     chunksize=None, pcon_index="index", n_con=None,
                     dtype=None, nz=None, nz_full=None):
        """
        Puts together particle outputs in space (for ENDLESS patches) and in time
        creating one big 5-dimensional numpy array in return, with the axes being
        time, x, y, z, n_con
        where n_con is the number of the particle

        t_ini, t_end: int
            initial and final time steps
        simulation: lespy.Simulation object
            simulation to consider when processing the results
        apply_to_z: bool
            If true, z_function is applied to each instantaneous 3D patch before
            merging into final array and final array has shape (time, x, y[, n_con]).
            If false, full 3D patch is merged and final array has shape (time,x,y,z[,n_con]).
            Currently only tested for endless (pcon_jt outputs).
        z_function: function
            Function (should be manually set to apply on z axis) to be applied on z axis
            if apply_to_z in True.
        dtype: numpy.type, python.type
            type to initialize the array with.
        """
        from .. import utils, routines, io
        import numpy as np
        import xarray as xr

        if simulation:
            sim=simulation
        else:
            raise ValueError("You need to provide a simulation object here")

        #--------------
        # Only these types of file have pcon output (maybe con_t also)
        label = 'pcon_jt'
        #--------------

        #--------------
        # Adjust intervals
        if type(times)!=type(None):
            bins = self.binaries.loc[:,label].dropna(how='all').loc[times]
        else:
            if t_end is None:
                bins = self.binaries.loc[t_ini:, label].dropna(how='all')
            else:
                bins = self.binaries.loc[t_ini:t_end, label].dropna(how='all')
        #--------------

        #---------
        # For too-large sizes, it's better to integrate over z patch-by-patch
        if type(nz)==type(None):
            nz=sim.domain.nz_tot-1
        #---------

        pcons = []
        for i, fname in enumerate(bins):
            print(i, fname)
            con = io.readBinary(fname, simulation=sim,
                                as_DA=True, nz=nz, nz_full=nz_full, 
                                pcon_index=pcon_index, n_con=n_con)

            #--------
            # Here we apply the z_function to 4D patch, making it 3D (x, y, n_con), and merge
            if apply_to_z:
                pcons.append(z_function(con))
            else:
                pcons.append(con)
            #--------


        #--------
        out = xr.concat(pcons, dim="itime").assign_coords(itime=bins.index.tolist())
        from ..utils import add_units
        out = add_units(out)
        if chunksize is not None:
            print("Chunking with Dask")
            out = out.chunk(dict(itime=chunksize))
        return out


    def compose_uvw(self, simulation=None, times=None, t_ini=0, t_end=None,
                    apply_to_z=False, z_function=lambda x: x[:,:,0], 
                    chunksize=None,
                    dtype=None, nz=None, nz_full=None):
        """
        Puts together everything in time
        Output array indexes are
        time, x, y, z

        Parameters
        ----------
        self: lp.Output
        times: slice
            slice with which times to read.
        simulation: lp.Simulation
            to be used as base
        """
        from .. import utils, routines, io
        import numpy as np
        import xarray as xr

        if simulation:
            sim=simulation
        else:
            raise ValueError("You need to provide a simulation object here")

        #--------------
        # Only these types of file we deal with here
        label = 'uvw_jt'
        Nx=sim.domain.nx
        #--------------

        #--------------
        # Adjust intervals
        if type(times)!=type(None):
            bins = self.binaries.loc[:,label].dropna(how='all').loc[times]
        else:
            if t_end is None:
                bins = self.binaries.loc[t_ini:, label].dropna(how='all')
            else:
                bins = self.binaries.loc[t_ini:t_end, label].dropna(how='all')
        #--------------
        
        #---------
        # Definition of output with time, x, y[ and z]
        if type(nz)==type(None):
            nz=sim.domain.nz_tot-1
        #---------

        #---------
        # Iterate between uvw_jt files
        u,v,w, = [],[],[],
        for i,col in enumerate(bins):
            if not isinstance(col, str): continue
            print(i, col)
            aux = io.readBinary(col, simulation=sim, as_DA=True, nz=nz, nz_full=nz_full)

            #------
            # Reduce z coordinate if theres a z_function
            if apply_to_z:
                aux = [ z_function(el) for el in aux ]
            #------

            u.append(aux[0])
            v.append(aux[1])
            w.append(aux[2])
            #---------
        #---------

        #---------
        # Passes from numpy.array to xarray.DataArray, so that the coordinates go with the data
        u = xr.concat(u, dim="itime").assign_coords(itime=bins.index.tolist())
        v = xr.concat(v, dim="itime").assign_coords(itime=bins.index.tolist())
        w = xr.concat(w, dim="itime").assign_coords(itime=bins.index.tolist())
        out = [u, v, w]
        from ..utils import add_units
        for i in range(len(out)):
            out[i] = add_units(out[i])
            if chunksize is not None:
                print("Chunking with Dask")
                out[i] = out[i].chunk(dict(itime=chunksize))
        #---------

        return out


    def compose_theta(self, simulation=None, times=None, t_ini=0, t_end=None,
                      apply_to_z=False, z_function=lambda x: x[:,:,0], 
                      chunksize=None,
                      dtype=None, nz=None, nz_full=None):
        """
        Puts together everything in time
        Output array indexes are
        time, x, y, z

        Parameters
        ----------
        self: lp.Output
        times: slice
            slice with which times to read.
        simulation: lp.Simulation
            to be used as base
        """
        from .. import utils, routines, io
        import numpy as np
        import xarray as xr

        if simulation:
            sim=simulation
        else:
            raise ValueError("You need to provide a simulation object here")

        #--------------
        # Only these types of file we deal with here
        label = 'theta_jt'
        Nx=sim.domain.nx
        #--------------

        #--------------
        # Adjust intervals
        if type(times)!=type(None):
            bins = self.binaries.loc[:,label].dropna(how='all').loc[times]
        else:
            if t_end is None:
                bins = self.binaries.loc[t_ini:, label].dropna(how='all')
            else:
                bins = self.binaries.loc[t_ini:t_end, label].dropna(how='all')
        #--------------
        
        #---------
        # Definition of output with time, x, y[ and z]
        if type(nz)==type(None):
            nz=sim.domain.nz_tot-1
        #---------

        #---------
        # Iterate between uvw_jt files
        theta = []
        for i,col in enumerate(bins):
            if not isinstance(col, str): continue
            print(i, col)
            aux = io.readBinary(col, simulation=sim, as_DA=True, nz=nz, nz_full=nz_full)

            #------
            # Reduce z coordinate if theres a z_function
            if apply_to_z:
                aux = [ z_function(el) for el in aux ]
            #------

            theta.append(aux)
            #---------
        #---------

        #---------
        # Passes from numpy.array to xarray.DataArray, so that the coordinates go with the data
        out = xr.concat(theta, dim="itime").assign_coords(itime=bins.index.tolist())
        from ..utils import add_units
        out = add_units(out)
        if chunksize is not None:
            print("Chunking with dask")
            out = out.chunk(dict(itime=chunksize))

        #---------

        return out



    def compose_uv0(self, times=None, t_ini=0, t_end=None, simulation=None, trim=True, as_dataarray=True):
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
        from .. import utils, io
        import numpy as np
        import pandas as pd

        if simulation:
            sim=simulation
        else:
            raise ValueError("You need to provide a simulation object here")

        #--------------
        # Only these types of file we deal with here
        label = ['uv0_jt']
        #--------------

        #--------------
        # Adjust intervals
        if type(times)!=type(None):
            bins = self.binaries.loc[:, label].dropna(axis=1, how='all').loc[times]
        else:
            if t_end is None:
                bins = self.binaries.loc[t_ini:, label].dropna(axis=1, how='all')
            else:
                bins = self.binaries.loc[t_ini:t_end, label].dropna(axis=1, how='all')
        #--------------

        #---------
        # Definition of output with time, x, y[ and z]
        print('Creating 2 arrays of {}, {}, {}...'.format(len(bins), sim.domain.nx, sim.ny), end='')
        u0 = np.full((len(bins), sim.nx, sim.ny), np.nan)
        v0 = np.full((len(bins), sim.nx, sim.ny), np.nan)
        print(' done.')
        #---------

        #---------
        # Iterates on the timestamps to put the time-indexed array together
        for i, tstep in enumerate(bins.index):
            iSeries = bins.loc[tstep]

            #---------
            # Iterate between uvw_jt, vel_t and theta_jt
            # Currently only works with uvw_jt (I think)
            for col in iSeries:
                if not isinstance(col, str): continue
                print(col)
                aux = io.readBinary(col, simulation=sim)
                ui,vi = aux
                u0[i] = ui
                v0[i] = vi
            #---------
        #---------

        #---------
        # Trims the extra node(s) at the end of the x coordinate
        if trim:
            u0 = u0[:, :sim.nx]
            v0 = v0[:, :sim.nx]
        #---------

        #---------
        # Passes from numpy.array to xarray.DataArray, so that the coordinates go with the data
        if as_dataarray:
            u0 = utils.get_DA(u0, simulation=sim, dims=['itime', 'x', 'y'], itime=bins.index.tolist())
            v0 = utils.get_DA(v0, simulation=sim, dims=['itime', 'x', 'y'], itime=bins.index.tolist())
        #---------

        return u0, v0


    def compose_var(self, simulation=None, name="", times=None, t_ini=0, t_end=None,
                    apply_to_z=False, z_function=lambda x: x[:,:,0], 
                    chunksize=None,
                    dtype=None, nz=None, nz_full=None):
        """
        Puts together everything in time
        Output array indexes are
        time, x, y, z

        Parameters
        ----------
        self: lp.Output
        times: slice
            slice with which times to read.
        simulation: lp.Simulation
            to be used as base
        """
        from .. import utils, routines, io
        import numpy as np
        import xarray as xr
        import pandas as pd
        from glob import glob
        from os import path

        if simulation:
            sim=simulation
        else:
            raise ValueError("You need to provide a simulation object here")

        #----------
        # list name_jt files
        name_jt = sorted(glob(path.join(self.oppath, f'{name}_jt*bin')))
        if name_jt: print(f'Reading {name}_jt files')
        ndtime = np.vectorize(utils.nameParser)(name_jt)
        bins = pd.Series(data=name_jt, index=ndtime)
        #--------------

        #--------------
        # Adjust intervals
        if type(times)!=type(None):
            bins = bins.loc[times]
        else:
            if t_end is None:
                bins = bins.loc[t_ini:]
            else:
                bins = bins.loc[t_ini:t_end]
        #--------------
 
        #---------
        # Definition of output with time, x, y[ and z]
        if type(nz)==type(None):
            nz=sim.domain.nz_tot-1
        #---------

        #---------
        # Iterate between uvw_jt files
        var = []
        for i,col in enumerate(bins):
            if not isinstance(col, str): continue
            print(i, col)
            aux = io.readBinary(col, simulation=sim, as_DA=True, nz=nz, nz_full=nz_full)

            #------
            # Reduce z coordinate if theres a z_function
            if apply_to_z:
                aux = [ z_function(el) for el in aux ]
            #------

            var.append(aux)
            #---------
        #---------

        #---------
        # Passes from numpy.array to xarray.DataArray, so that the coordinates go with the data
        out = xr.concat(var, dim="itime").assign_coords(itime=bins.index.tolist())
        from ..utils import add_units
        out = add_units(out)
        if chunksize is not None:
            print("Chunking with dask")
            out = out.chunk(dict(itime=chunksize))

        #---------

        return out






class Output_sp(object):
    """
    Class that holds the output of the model and processes it.
    It's diferent from the Simulation class, but it can't do stuff without it.
    """
    def __init__(self, oppath, apply_basename=False, n_cons=[1], separate_ncon=True, verbose=False):
        """
        Lists every output file from a simulation (so far only uvw_jt, theta_jt and pcon_jt)

        oppath can be either the directory of the output, or one of the outuput binary files
        """
        from os import path
        from glob import glob
        import pandas as pd
        from .. import utils

        #----------
        # We use pandas to organize all the files (very handy)
        opfiles = pd.DataFrame(columns=['field_', 'pcon_jt'])
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
        # List field_ files
        field_ = sorted(glob(path.join(oppath, 'field_*out')))
        if field_: print('field_', end=' ')
        for fname in field_:
            ndtime = utils.nameParser(fname)
            opfiles.loc[ndtime, 'field_'] = fname
        #----------

        #----------
        # list pcon_jt files (each entry is a list, since there can be lots for each timestep
        pcon_jt = sorted(glob(path.join(oppath, 'pcon_jt*bin')))
        if pcon_jt: print('pcon_jt', end=' ')
        if apply_basename:
            pcon_jt = map(path.basename, pcon_jt)

        for fname in pcon_jt:
            ndtime = utils.nameParser(fname)
            opfiles.loc[ndtime, 'pcon_jt'] = fname
        #----------

        self.binaries = opfiles.sort_index()
        print('... done reading and organizing.')
        return


    def compose_pcon(self, times=None, t_ini=0, t_end=None, simulation=None,
                     apply_to_z=False, z_function=lambda x: x[:,:,0], 
                     chunksize=None, pcon_index="w_r", 
                     dtype=None, nz=None, nz_full=None):
        """
        Puts together particle outputs in space (for ENDLESS patches) and in time
        creating one big 5-dimensional numpy array in return, with the axes being
        time, x, y, z, n_con
        where n_con is the number of the particle

        t_ini, t_end: int
            initial and final time steps
        simulation: lespy.Simulation object
            simulation to consider when processing the results
        apply_to_z: bool
            If true, z_function is applied to each instantaneous 3D patch before
            merging into final array and final array has shape (time, x, y[, n_con]).
            If false, full 3D patch is merged and final array has shape (time,x,y,z[,n_con]).
            Currently only tested for endless (pcon_jt outputs).
        z_function: function
            Function (should be manually set to apply on z axis) to be applied on z axis
            if apply_to_z in True.
        dtype: numpy.type, python.type
            type to initialize the array with.
        """
        from .. import utils, routines, io
        import numpy as np

        if simulation:
            sim=simulation
        else:
            raise ValueError("You need to provide a simulation object here")

        #--------------
        # Only these types of file have pcon output (maybe con_t also)
        label = 'pcon_jt'
        #--------------

        #--------------
        # Adjust intervals
        if type(times)!=type(None):
            cons = self.binaries.loc[:,label].dropna(how='all').loc[times]
        else:
            if t_end is None:
                cons = self.binaries.loc[t_ini:, label].dropna(how='all')
            else:
                cons = self.binaries.loc[t_ini:t_end, label].dropna(how='all')
        #--------------

        #---------
        # For too-large sizes, it's better to integrate over z patch-by-patch
        if type(nz)==type(None):
            nz=sim.domain.nz_tot-1
        if apply_to_z:
            print('Creating array of ',(len(cons), sim.nx, sim.ny, sim.n_con))
            pcons = np.full((len(cons), sim.nx, sim.ny, sim.n_con), np.nan, dtype=dtype)
            dims = ['itime', 'x', 'y', pcon_index]
        else:
            print('Creating array of ',(len(cons), sim.nx, sim.ny, nz, sim.n_con))
            pcons = np.full((len(cons), sim.nx, sim.ny, nz, sim.n_con), np.nan, dtype=dtype)
            dims = ['itime', 'x', 'y', 'z_u', pcon_index]
        #---------

        for i, fname in enumerate(cons):
            print(i, fname)
            con = io.readBinary(fname, simulation=sim, n_con=sim.n_con, as_DA=False, nz=nz, nz_full=nz_full)

            #--------
            # Here we apply the z_function to 4D patch, making it 3D (x, y, n_con), and merge
            if apply_to_z:
                pcons[i,:,:,:] = z_function(con[:sim.nx,:,:,:])
            else:
                pcons[i,:,:,:,:] = con[:sim.nx,:,:,:]
            #--------


        from ..utils import add_units
        out = utils.get_DA(pcons, simulation=sim, dims=dims, time=cons.index.tolist())
        out = add_units(out)
        if chunksize is not None:
            out = out.chunk(dict(itime=chunksize))
        return out


    def compose_uvw(self, simulation=None, times=None, t_ini=0, t_end=None,
                    apply_to_z=False, z_function=lambda x: x[:,:,0], 
                    chunksize=None,
                    dtype=None, nz=None, nz_full=None):
        """
        Puts together everything in time
        Output array indexes are
        time, x, y, z

        Parameters
        ----------
        self: lp.Output
        times: slice
            slice with which times to read.
        simulation: lp.Simulation
            to be used as base
        """
        from .. import utils, routines, io
        import numpy as np
        import pandas as pd
        import xarray as xr

        if simulation:
            sim=simulation
        else:
            raise ValueError("You need to provide a simulation object here")

        #--------------
        # Only these types of file we deal with here
        label = 'field_'
        Nx=sim.domain.nx
        #--------------

        #--------------
        # Adjust intervals
        if type(times)!=type(None):
            bins = self.binaries.loc[:,label].dropna(how='all').loc[times]
        else:
            if t_end is None:
                bins = self.binaries.loc[t_ini:, label].dropna(how='all')
            else:
                bins = self.binaries.loc[t_ini:t_end, label].dropna(how='all')
        #--------------
        
        #---------
        # Definition of output with time, x, y[ and z]
        if type(nz)==type(None):
            nz=sim.domain.nz_tot-1
        #---------

        #---------
        # Prepare for reading fortran file
        vdict = utils.get_dicts(sim, nz=nz)
        vdict["v"]["jumpto"] = sim.nx*sim.ny*(sim.nz_tot-1)*4*1
        vdict["w"]["jumpto"] = sim.nx*sim.ny*(sim.nz_tot-1)*4*2
        vlist = [vdict["u"], vdict["v"], vdict["w"]]
        #---------

        #---------
        # Iterate between field_ files
        u,v,w, = [],[],[],
        for i,col in enumerate(bins):
            if not isinstance(col, str): continue
            print(i, col)
            aux = io.fortran2xr(col, vlist, dtype=np.float32)
            aux = [ a*sim.u_scale for a in aux ]

            #------
            # Reduce z coordinate if theres a z_function
            if apply_to_z:
                aux = [ z_function(el) for el in aux ]
            #------

            #---------
            u.append(aux[0])
            v.append(aux[1])
            w.append(aux[2])
            #---------
        #---------

        #---------
        # Concatenates the DataArrays
        from ..utils import add_units
        u = xr.concat(u, dim="itime").assign_coords(itime=bins.index.tolist())
        v = xr.concat(v, dim="itime").assign_coords(itime=bins.index.tolist())
        w = xr.concat(w, dim="itime").assign_coords(itime=bins.index.tolist())
        out = [u, v, w]
        for i in range(len(out)):
            out[i] = add_units(out[i])
            if chunksize is not None:
                out[i] = out[i].chunk(dict(itime=chunksize))
        #---------

        return out


    def compose_theta(self, simulation=None, times=None, t_ini=0, t_end=None,
                    apply_to_z=False, z_function=lambda x: x[:,:,0], 
                    chunksize=None,
                    dtype=None, nz=None, nz_full=None):
        """
        Puts together everything in time
        Output array indexes are
        time, x, y, z

        Parameters
        ----------
        self: lp.Output
        times: slice
            slice with which times to read.
        simulation: lp.Simulation
            to be used as base
        """
        from .. import utils, routines, io
        import numpy as np
        import pandas as pd
        import xarray as xr

        if simulation:
            sim=simulation
        else:
            raise ValueError("You need to provide a simulation object here")

        #--------------
        # Only these types of file we deal with here
        label = 'field_'
        Nx=sim.domain.nx
        #--------------

        #--------------
        # Adjust intervals
        if type(times)!=type(None):
            bins = self.binaries.loc[:,label].dropna(how='all').loc[times]
        else:
            if t_end is None:
                bins = self.binaries.loc[t_ini:, label].dropna(how='all')
            else:
                bins = self.binaries.loc[t_ini:t_end, label].dropna(how='all')
        #--------------
        
        #---------
        # Definition of output with time, x, y[ and z]
        if type(nz)==type(None):
            nz=sim.domain.nz_tot-1
        #---------

        #---------
        # Prepare for reading fortran file
        vdict = utils.get_dicts(sim, nz=nz)
        vdict["θ"]["jumpto"] = sim.nx*sim.ny*(sim.nz_tot-1)*4*4
        vlist = [vdict["θ"]]
        #---------

        #---------
        # Iterate between field_ files
        theta = []
        for i,col in enumerate(bins):
            if not isinstance(col, str): continue
            print(i, col)
            aux, = io.fortran2xr(col, vlist, dtype=np.float32)
            aux = 2.*sim.t_init - aux*sim.t_scale

            #------
            # Reduce z coordinate if theres a z_function
            if apply_to_z:
                aux = [ z_function(el) for el in aux ]
            #------

            theta.append(aux)
            #---------
        #---------

        #---------
        # Concatenates in itime
        from ..utils import add_units
        theta = xr.concat(theta, dim="itime").assign_coords(itime=bins.index.tolist())
        out = theta
        for i in range(len(out)):
            out = add_units(out)
            if chunksize is not None:
                out = out.chunk(dict(itime=chunksize))
        #---------

        return out






class Output_old(object):
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

        self.binaries = opfiles.sort_index()
        print('... done reading and organizing.')
        return

    def compose_edls(self, times=None, t_ini=0, t_end=None, simulation=None, as_dataarray=True,
            apply_to_z=False, z_function=lambda x: x[:,:,0], dtype=None, trim=True):
        """
        Puts together particle outputs in space (for ENDLESS patches) and in time

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
        # Adjust intervals
        labels = ['con_tt']
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
        else:# 'con_tt' in cons.columns.tolist():
            cons = cons.con_tt.dropna()
       
            #---------
            # For too-large sizes, it's better to integrate over z patch-by-patch
            print(self)
            print(self.n_cons)
            shape=(len(cons), len(self.n_cons))
            print('Creating object array of shape', shape)
            pcons = np.full(shape, np.nan, dtype='object')
            #---------

            #---------
            # First we find the max and min row and col
            rows = []
            cols = []
            for patches in cons:
                for n_i in range(len(self.n_cons)):
                    for patch in patches[n_i]:
                        ndtime, ncon, row, col = utils.nameParser(patch)
                        rows.append(row)
                        cols.append(col)
            min_row, min_col = min(rows), min(cols)
            #---------

            #--------
            # We iterate over the different sizes of pcon.
            # Here n_i is a local integer and n_cons is a lit with the files for one size
            for n_i, n_con in enumerate(self.n_cons):
                print(n_i, n_con)

                #---------
                # Now we iterate again over the time steps to read the files of one pcon size
                for i, ndtime in enumerate(cons.index):
                    print(i, ndtime, n_con)
                    clist = []
                    for patch in cons.loc[ndtime][n_i]:
                        ndtime, ncon, row, col = utils.nameParser(patch)
                        con = routines.readBinary2(patch, simulation=sim, read_pcon=True, only_pcon=True, trim=trim)
                        
                        #X = sim.domain.x + sim.lx*(col-min_col)
                        #Y = sim.domain.y + sim.ly*(row-min_row)
                        X = sim.domain.x + sim.lx*(col)
                        Y = sim.domain.y + sim.ly*(row)

                        #--------
                        # Here we apply the z_function to 3D patch, making it 2D (x, y), and merge
                        if apply_to_z:
                            pcon = z_function(con[:sim.nx,:,:])
                            if as_dataarray:
                                pcon = sim.DataArray(pcon, dims=['x', 'y'], coords=dict(x=X, y=Y))
                        else:
                            pcon = con[:sim.nx,:,:]
                            if as_dataarray:
                                pcon = sim.DataArray(pcon, dims=['x', 'y', 'z'], coords=dict(x=X, y=Y, z=sim.domain.z_u))

                        clist.append(pcon)
                        #--------
                    pcons[i, n_i] = clist
                #---------
            #--------
            if as_dataarray==True:
                import xarray as xr
                print(pcons.shape)
                pcons = xr.DataArray(pcons, dims=['time', 'size'], 
                        coords=dict(time=cons.index.tolist(), size=range(sim.n_con)))
            return pcons



    def compose_pcon(self, times=None, t_ini=0, t_end=None, simulation=None, as_dataarray=True,
            apply_to_z=False, z_function=lambda x: x[:,:,0], dtype=None, trim=True):
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
            raise NotImplementedError('Use compose_edls() function')
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
                dims=['time', 'x', 'y', 'size']
            else:
                print('Creating array of ',(len(cons), sim.nx, sim.ny, sim.nz_tot, sim.n_con))
                dims=['time', 'x', 'y', 'z', 'size']
                pcons = np.full((len(cons), sim.nx, sim.ny, sim.nz_tot, sim.n_con), np.nan, dtype=dtype)
            #---------

            for i, fname in enumerate(cons):
                print(i, cons.index[i], fname)
                con = routines.readBinary2(fname, simulation=sim, n_con=sim.n_con, only_pcon=True)

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
            #output = utils.get_dataarray(pcons, simulation=sim, with_time=cons.index.tolist())
            output = utils.get_DA(pcons, simulation=sim, dims=dims, time=cons.index.tolist())
            if type(output)==list:
                return output[0]
            else: return output
        else:
            if type(pcons)==list:
                return pcons[0]
            else:
                return pcons


    def compose_uvwT(self, times=None, t_ini=0, t_end=None, simulation=None, trim=True, as_dataarray=True,
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
        if trim:
            Nx=sim.domain.nx
        else:
            Nx=sim.domain.ld
        #---------
        # Definition of output with time, x, y[ and z]
        if apply_to_z:
            print('Creating 4 arrays of {}, {}, {}...'.format(len(bins), Nx, sim.ny), end='')
            u = np.full((len(bins), Nx, sim.ny), np.nan)
            v = np.full((len(bins), Nx, sim.ny), np.nan)
            w = np.full((len(bins), Nx, sim.ny), np.nan)
            T = np.full((len(bins), Nx, sim.ny), np.nan)
            dims_u = ['time', 'x', 'y']
            dims_w = ['time', 'x', 'y']
        else: 
            print('Creating 4 arrays of {}, {}, {}, {}...'.format(len(bins), Nx, sim.ny, sim.nz_tot), end='')
            u = np.full((len(bins), Nx, sim.ny, sim.nz_tot), np.nan)
            v = np.full((len(bins), Nx, sim.ny, sim.nz_tot), np.nan)
            w = np.full((len(bins), Nx, sim.ny, sim.nz_tot), np.nan)
            T = np.full((len(bins), Nx, sim.ny, sim.nz_tot), np.nan)
            dims_u = ['time', 'x', 'y', 'z_u']
            dims_w = ['time', 'x', 'y', 'z_w']
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
                aux = routines.readBinary2(col, simulation=sim, read_pcon=False, trim=trim)

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
            #u,v,w,T = utils.get_dataarray([u, v, w, T], simulation=sim, with_time=bins.index.tolist())
            u,v,w,T = [ utils.get_DA(u, simulation=sim, dims=dims_u, time=bins.index.tolist()),
                utils.get_DA(v, simulation=sim, dims=dims_u, time=bins.index.tolist()),
                utils.get_DA(w, simulation=sim, dims=dims_w, time=bins.index.tolist()),
                utils.get_DA(T, simulation=sim, dims=dims_u, time=bins.index.tolist()), ]
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
        if trim:
            du = np.full((len(bins), sim.domain.nx, sim.ny), np.nan)
            dv = np.full((len(bins), sim.domain.nx, sim.ny), np.nan)
            dw = np.full((len(bins), sim.domain.nx, sim.ny), np.nan)
        else:
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
                aux = routines.readBinary2(col, simulation=sim, read_pcon=False, trim=trim)
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
        def par_compose_uvwT(x):
            return out.compose_uvwT(simulation=sim, t_ini=x[0], t_end=x[1], **kwargs)
        pool = Pool(processes=nprocs)
        outs = pool.map(par_compose_uvwT, bounds)
        return zip(*outs)

