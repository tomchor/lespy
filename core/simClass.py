
class Simulation(object):
    """class for simulation parameters"""
    def __init__(self, from_file=None, domain=None, timelength=None, inversion_depth=None, **kwargs):
        import numpy as np
        from .. import physics
        #------------
        # Start simulation from param.nml file
        if from_file!=None:
            if isinstance(from_file, str):
                aux = sim_from_file(from_file)
                self.__dict__.update(vars(aux))
                return
            else:
                raise ValueError('fromfile keyword should be path to param.nml')
        #------------

        #------------
        # Start simulation from call arguments
        else:
            self.domain = domain
            self.timelength = timelength
            self.inversion_depth = inversion_depth
            self.inv_depth = inversion_depth
            self.__dict__.update(kwargs)
        #------------

        #-----
        try:
            self.u_scale*=1
        except:
            self.u_scale=self.u_star
        #-----

        #-----
        self.w_star = self.get_w_star()
        self.vel_settling=np.array(self.vel_settling)
        try:
            self.relax_freq=np.array(self.relax_freq)
        except:
            pass
        self.droplet_sizes=physics.get_dropletSize(self.vel_settling, nominal=True, nowarning=True).astype(int)
        if type(self.droplet_sizes)!=np.ndarray:
            self.droplet_sizes=np.array([self.droplet_sizes])
        #-----

        #-----
        # Make backwards compatible with older version of code
        try:
            self.s_flag=self.theta_flag
        except AttributeError:
            self.theta_flag=self.s_flag
        #-----

        #-----
        # Make it easy to access resolution
        self.Δx=self.domain.dx
        self.Δy=self.domain.dy
        self.Δz=self.domain.dz
        #-----


    def check(self, full=True):
        """
        Check important characteristics of the simulation
        """
        CFL_x=self.u_scale*self.dt/self.domain.dx
        print('CFL (u_scale*dt/dx)          : {:.2e}'.format(CFL_x))
        print('dx/dz                        : {:2.1f}\t\t{}'.format(self.domain.dx/self.domain.dz,'-- Should be < 5 in practice'))
        print('lx/z_inv                     : {:2.1f}\t\t{}'.format(self.domain.lx/self.inversion_depth,'-- Should be > 6. At *least* 4.'))
        divs = []
        for i in range(2,140):
            if self.domain.nz%i == 0:
                divs.append(i)
        print('Nz = {:03d} and is divisible by : {}'.format(self.domain.nz, divs))
        if full:
            print('Coriolis timescale           : {:1.1e} timesteps'.format(int(1./self.freq_coriolis/self.dt)))

    def get_w_star(self):
        """Calculates the convective scale"""
        from .. import physics as phys
        return phys.w_star(self)


    def DataArray(self, array, **kwargs):
        """
        Creates a DataArray specifically for this simulation.
        
        Currently does not work with xarray yet because of different domain.
        """
        #import xarray as xr
        from ..utils import get_DA
        da = get_DA(array, simulation=self, **kwargs)
        return da

    def to_hours(self, timesteps, to_label=False):
        """Transforms from simulation timesteps to hours"""
        out = timesteps*self.dt/(60*60)
        if to_label:
            out = [ '{:.2f} hours'.format(el) for el in out ]
        return out

    def __str__(self):
        buff='Simulation Parameters\n'+ '-'*21
        buff += '\nEndless    : {}\n'.format(self.flag_endless)
        buff += 'dt:        : {} s\n'.format(self.dt)
        buff+= self.domain.__str__()
        return buff
    def __repr__(self):
        aux = """<lespy.Simulation object>. Domain: {}""".format(self.domain.__repr__())
        return aux




class hSimulation(object):
    """class for simulation parameters"""
    def __init__(self, from_file=None, domain=None, timelength=None, inversion_depth=None, **kwargs):
        import numpy as np
        from .. import physics
        from ..utils import paramParser, find_in_tree, nameParser
        from .dmClass import hDomain as hDom

        #------------
        # Start simulation from param.nml file
        if isinstance(from_file, str):
            params = paramParser(from_file)
            self.__dict__.update(params)
        else:
            raise ValueError('fromfile keyword should be path to param.nml')
        #-----

        #-----
        dmn = hDom(nx=params['nx'], ny=params['ny'], nz=params['nz_tot'],
                   lx=params['lx_tot'], ly=params['ly_tot'], lz=params['lz_tot'],
                   environment=params['environment'])
        self.domain = dmn
        #-----

        #-----
        if self.environment=="ocean":
            self.ocean_flag=True
        else:
            self.ocean_flag=False
        self.inversion_depth = self.z_i * self.prop_mixed
        #-----

        #-----
        self.w_star = self.get_w_star()
        self.vel_settling = np.array(self.vel_settling)
        self.droplet_sizes = physics.get_dropletSize(self.vel_settling, nominal=True, nowarning=True).astype(int)
        if type(self.droplet_sizes)!=np.ndarray:
            self.droplet_sizes=np.array([self.droplet_sizes])
        #-----

        #-----
        # Make backwards compatible with older version of code
        try:
            self.s_flag=self.theta_flag
        except AttributeError:
            self.theta_flag=self.s_flag
        #-----

        #-----
        # Make it easy to access resolution
        self.Δx=self.domain.dx
        self.Δy=self.domain.dy
        self.Δz=self.domain.dz
        #-----


    def check(self, full=True):
        """
        Check important characteristics of the simulation
        """
        CFL_x=self.u_scale*self.dt/self.domain.dx
        print('CFL (u_scale*dt/dx)          : {:.2e}'.format(CFL_x))
        print('dx/dz                        : {:2.1f}\t\t{}'.format(self.domain.dx/self.domain.dz,'-- Should be < 5 in practice'))
        print('lx/z_inv                     : {:2.1f}\t\t{}'.format(self.domain.lx/self.inversion_depth,'-- Should be > 6. At *least* 4.'))
        divs = []
        for i in range(2,140):
            if self.domain.nz%i == 0:
                divs.append(i)
        print('Nz = {:03d} and is divisible by : {}'.format(self.domain.nz, divs))
        if full:
            print('Coriolis timescale           : {:1.1e} timesteps'.format(int(1./self.freq_coriolis/self.dt)))

    def get_w_star(self):
        """Calculates the convective scale"""
        from .. import physics as phys
        return phys.w_star(self)


    def DataArray(self, array, **kwargs):
        """
        Creates a DataArray specifically for this simulation.
        
        Currently does not work with xarray yet because of different domain.
        """
        #import xarray as xr
        from ..utils import get_DA
        da = get_DA(array, simulation=self, **kwargs)
        return da

    def to_hours(self, timesteps, to_label=False):
        """Transforms from simulation timesteps to hours"""
        out = timesteps*self.dt/(60*60)
        if to_label:
            out = [ '{:.2f} hours'.format(el) for el in out ]
        return out

    def __str__(self):
        buff='Simulation Parameters\n'+ '-'*21
        buff += '\nEndless    : {}\n'.format(self.flag_endless)
        buff += 'dt:        : {} s\n'.format(self.dt)
        buff+= self.domain.__str__()
        return buff
    def __repr__(self):
        aux = """<lespy.Simulation object>. Domain: {}""".format(self.domain.__repr__())
        return aux





class Simulation_sp(object):
    """class for simulation parameters"""
    def __init__(self, from_file=None, **kwargs):
        import numpy as np
        from .. import physics
        #------------
        # Start simulation from param.nml file
        from ..utils import paramParser
        params = paramParser(from_file)
        #------------

        #------------
        params["nx"] = params["nxt"]
        params["ny"] = params["nyt"]
        params["nz_tot"] = params["nzt"]
        params["lx"] = params["lx_tot"]
        params["ly"] = params["ly_tot"]
        #------------


        #------------
        aux = sim_from_file(from_file, params=params)
        self.__dict__.update(vars(aux))
        #------------

        try:
            self.u_scale*=1
        except:
            self.u_scale=self.u_star

        self.w_star = self.get_w_star()
        self.vel_settling=np.array(self.vel_settling)
        try:
            self.relax_freq=np.array(self.relax_freq)
        except:
            pass
        self.droplet_sizes=physics.get_dropletSize(self.vel_settling, nominal=True, nowarning=True).astype(int)
        if type(self.droplet_sizes)!=np.ndarray:
            self.droplet_sizes=np.array([self.droplet_sizes])

        try:
            self.s_flag=self.theta_flag
        except AttributeError:
            self.theta_flag=self.s_flag


    def check(self, full=True):
        """
        Check important characteristics of the simulation
        """
        CFL_x=self.u_scale*self.dt/self.domain.dx
        print('CFL (u_scale*dt/dx)          : {:.2e}'.format(CFL_x))
        print('dx/dz                        : {:2.1f}\t\t{}'.format(self.domain.dx/self.domain.dz,'-- Should be < 5 in practice'))
        print('lx/z_inv                     : {:2.1f}\t\t{}'.format(self.domain.lx/self.inversion_depth,'-- Should be > 6. At *least* 4.'))
        divs = []
        for i in range(2,140):
            if self.domain.nz%i == 0:
                divs.append(i)
        print('Nz = {:03d} and is divisible by : {}'.format(self.domain.nz, divs))
        if full:
            print('Coriolis timescale           : {:1.1e} timesteps'.format(int(1./self.freq_coriolis/self.dt)))

    def get_w_star(self):
        """Calculates the convective scale"""
        from .. import physics as phys
        return phys.w_star(self)


    def DataArray(self, array, **kwargs):
        """
        Creates a DataArray specifically for this simulation.
        
        Currently does not work with xarray yet because of different domain.
        """
        #import xarray as xr
        from ..utils import get_DA
        da = get_DA(array, simulation=self, **kwargs)
        return da

    def to_hours(self, timesteps, to_label=False):
        """Transforms from simulation timesteps to hours"""
        out = timesteps*self.dt/(60*60)
        if to_label:
            out = [ '{:.2f} hours'.format(el) for el in out ]
        return out

    def __str__(self):
        buff='Simulation Parameters\n'+ '-'*21
        buff += '\nEndless    : {}\n'.format(self.flag_endless)
        buff += 'dt:        : {} s\n'.format(self.dt)
        buff+= self.domain.__str__()
        return buff
    def __repr__(self):
        aux = """<lespy.Simulation object>. Domain: {}""".format(self.domain.__repr__())
        return aux







class Simulation_old(object):
    """class for simulation parameters"""
    def __init__(self, from_file=None, domain=None, timelength=None, u_scale=None, inversion_depth=None, **kwargs):
        import numpy as np
        from .. import physics
        #------------
        # Start simulation from param.nml file
        if from_file!=None:
            if isinstance(from_file, str):
                aux = sim_from_file(from_file)
                self.__dict__.update(vars(aux))
                return
            else:
                raise ValueError('fromfile keyword should be path to param.nml')
        #------------

        #------------
        # Start simulation from call arguments
        else:
            self.domain = domain
            self.timelength = timelength
            self.u_scale = u_scale
            self.inversion_depth = inversion_depth
            self.inv_depth = inversion_depth
            self.__dict__.update(kwargs)
        #------------

        self.w_star = self.get_w_star()
        self.vel_settling=np.array(self.vel_settling)
        self.droplet_sizes=physics.get_dropletSize(self.vel_settling, nominal=True, nowarning=True).astype(int)
        if type(self.droplet_sizes)!=np.ndarray:
            self.droplet_sizes=np.array([self.droplet_sizes])

        try:
            self.s_flag=self.theta_flag
        except AttributeError:
            self.theta_flag=self.s_flag


    def check(self, full=True):
        """
        Check important characteristics of the simulation
        """
        CFL_x=self.u_scale*self.dt/self.domain.dx
        print('CFL (u_scale*dt/dx)          : {:.2e}'.format(CFL_x))
        print('dx/dz                        : {:2.1f}\t\t{}'.format(self.domain.dx/self.domain.dz,'-- Should be < 5 in practice'))
        print('lx/z_inv                     : {:2.1f}\t\t{}'.format(self.domain.lx/self.inversion_depth,'-- Should be > 6. At *least* 4.'))
        divs = []
        for i in range(2,140):
            if self.domain.nz%i == 0:
                divs.append(i)
        print('Nz = {:03d} and is divisible by : {}'.format(self.domain.nz, divs))
        if full:
            print('Coriolis timescale           : {:1.1e} timesteps'.format(int(1./self.freq_coriolis/self.dt)))

    def get_w_star(self):
        """Calculates the convective scale"""
        from .. import physics as phys
        return phys.w_star(self)


    def DataArray(self, array, timestamps=None, attrs=None, dims=None, coords=None):
        """
        Creates a DataArray specifically for this simulation.
        
        Currently does not work with xarray yet because of different domain.
        """
        import xarray as xr

        #--------
        # Define coordinates and dimensions if not provided
        if coords==None:
            if len(array.shape)==4:
                if timestamps==None:
                    import numpy as np
                    timestamps = np.arange(array.shape[0])
                coords = [ timestamps, self.domain.x, self.domain.y, self.domain.z ]
                dims = [ 'time', 'x', 'y', 'z' ]
            elif len(array.shape)==3:
                coords = [ self.domain.x, self.domain.y, self.domain.z ]
                dims = [ 'x', 'y', 'z' ]
            else:
                raise ValueError('Works only with 3D or 4D array.')
        #--------

        da = xr.DataArray(array, coords=coords, dims=dims, name=None, attrs=attrs, encoding=None, fastpath=False)
        return da

    def to_hours(self, timesteps, to_label=False):
        """Transforms from simulation timesteps to hours"""
        out = timesteps*self.dt/(60*60)
        if to_label:
            out = [ '{:.2f} hours'.format(el) for el in out ]
        return out

    def __str__(self):
        buff='Simulation Parameters\n'+ '-'*21
        buff += '\nEndless    : {}\n'.format(self.flag_endless)
        buff += 'dt:        : {} s\n'.format(self.dt)
        buff+= self.domain.__str__()
        return buff
    def __repr__(self):
        aux = """<lespy.Simulation object>. Domain: {}""".format(self.domain.__repr__())
        return aux







def sim_from_file(namelist, tlength_from_ke=False, check_ke_file=None, params=None):
    """Reads and parses namelist and then calls Simulation class"""
    from ..utils import paramParser, find_in_tree, nameParser
    from .dmClass import Domain as Dom
    from .dmClass import hDomain as hDom
    from os import path

    #---------
    # Reads parameters from param.nml and creates domain class
    if params is None:
        params = paramParser(namelist)
    try:
        dmn = Dom(nx=params['nx'], ny=params['ny'], nz=params['nz_tot'], nz_tot=params['nz_tot'],
                  lx=params['lx'], ly=params['ly'], lz=params['lz_tot'],
                  ocean_flag=params['ocean_flag'])
    except KeyError:
        dmn = hDom(nx=params['nx'], ny=params['ny'], nz=params['nz_tot'],
                   lx=params['lx_tot'], ly=params['ly_tot'], lz=params['lz_tot'],
                   environment=params['environment'])
    #---------
    
    #---------
    # Tries to find check_ke.out file
    namelist = path.abspath(namelist)
    if path.isfile(namelist):
        nml_dir = path.dirname(namelist)
    else:
        nml_dir = namelist

    if check_ke_file:
        kefile = path.abspath(check_ke_file)
    else:
        kefile = find_in_tree('check_ke.out', nml_dir)
        if len(kefile)==0:
            print('No check_ke file was found')
            kefile=None
        else:
            kefile = path.abspath(kefile[0])
    #---------

    #---------
    # This is most precise, but slow. Finds out complete length of simulation with check_ke.out file
    if tlength_from_ke:
        #---------
        # We try to open kefile then
        if kefile==None:
            print('Warning: getting length solely from param.nml, which can be flawed.')
            tlength = None
        else:
            import pandas as _pd
            print('Opening', kefile)
            try:
                kearray = _pd.read_csv(kefile, delim_whitespace=True, squeeze=True, index_col=0, header=None)
                kelast = kearray.index[-1]
                kecount = len(kearray.index)
                if kelast == kecount:
                    tlength = kelast
                else:
                    print('Warning: linescount in ke_check.out is different from index label in the file.')
                    print('Setting timelength to the line count of check_ke.')
                    tlength = kecount
            except (FileNotFoundError, OSError) as e:
                print("Coulndn't open check_ke.out.")
                tlength = None
        #---------
    #---------

    #---------
    # This is relatively less accurate, but much faster.
    else:
        from glob import glob
        print('Getting number of timesteps from last timestep of vel*.out output...', end='')
        
        out_dir = find_in_tree('output', nml_dir, type='d')
        if len(out_dir) == 0:
            tlength = None
            print('failed.')
        else:
            out_dir = out_dir[0]
            try:
                lastout = sorted(glob(path.join(out_dir, 'vel_*.out')))[-1]
                tstep = nameParser(lastout)
                if isinstance(tstep, list):
                    tstep = tstep[0]
            except IndexError:
                tstep = None
            tlength = tstep
            print('done.')
    #---------

    #---------
    # re-call to Simulation class but with all the parameters known
    out = Simulation(domain=dmn,
            check_ke=kefile,
            timelength=tlength,
            avglength=params['p_count'],
            inversion_depth=params['z_i']*params['prop_mixed'],
            T_scale=params['t_scale'],
            **params)
    #---------

    return out


