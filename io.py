import numpy as np
import xarray as xr

def empty_array(simulation, n_con=False, trimmed=True, nz_full=None):
    """ Creates an empty array with the shape dictated by simulation """
    sim=simulation

    #----
    # Write extra two rows in x or not
    if trimmed:
        Nx = sim.nx
    else:
        Nx = sim.domain.ld
    #----

    #----
    # Is it the new version of io (Jan 2019) or the version before?
    # Version before that read until nz_tot, after that only nz_tot-1 should be written
    if type(nz_full)==type(None):
        nz_full=sim.domain.nz_tot-1
    #----

    #----
    # Create blank array
    if n_con:
        blank=np.full((Nx, sim.ny, nz_full, sim.n_con), np.nan, dtype=np.float64)
    else:
        blank=np.full((Nx, sim.ny, nz_full), np.nan, dtype=np.float64)
    #----

    return blank


def write_to_les(array, fname, simulation=None, **kwargs):
    """
    Writes array into a file fname in a format that LES can easily understand
    """
    sim=simulation
    try:
        array[sim.nx:] = 0.
    except:
        pass
    array.T.tofile(fname, **kwargs)
    return


def read_aver(fname, simulation, squeeze=True, return_times=False, dims=[], unique=None, **kwargs):
    """
    Reads aver_* files from LES

    for reading aver_PCon.out: dims=["ndtime", "z", "index"]
    """
    import tempfile
    from .utils import fix_E
    sim=simulation

    #-----
    # Read file with the fix for wrong fortran formatting
    with tempfile.TemporaryFile("r+t") as ftemp:
        fix_E(fname, ftemp)
        aver=np.loadtxt(ftemp, **kwargs)
    #-----

    #-----
    # If pcon file, the shape is different
    if 'pcon' in fname.lower():
        aver = aver.reshape(-1, sim.n_con, aver.shape[-1], order='C').transpose(0,2,1)
        ndtimes = aver[:,0,0]
    else:
        ndtimes = aver[:,0]
    #-----

    #-----
    # Get the actual averages (without ndtime)
    aver = aver[:,1:]
    if squeeze:
        aver = np.squeeze(aver)
    #-----

    #-----
    # If dims is passed, return dataArray
    if dims:
        from . import utils
        aver = utils.get_DA(aver, simulation=sim, dims=dims, time=ndtimes)
    #-----

    #-----
    if (unique is not None):
        from . import utils
        aver = utils.unique_xr(aver, dim="ndtime")
    #-----

    if return_times:
        return ndtimes, aver
    else:
        return aver



def readBinary(fname, simulation=None, domain=None, n_con=None, as_DA=True, 
               pcon_index='index', nz=None, nz_full=None, dtype=np.float64):
    """
    Reads a binary file according to the simulation or domain object passed

    Parameters
    ----------
    fname: string
        path of the binary file you want to open
    simulation: lespy.Simulation
        simulation that contains important information to assemble the file. Mostly
        it's used to get number of points, if temperature is in the file or not and n_con.
    domain: lespy.Domain
        depending on what you're doing, just a domain file suffices
    read_pcon: bool
        if the file is vel_sc, try to read pcon after it finishes reading u,v,w,T.
        It just gives a warning if it fails
    n_con: int
        number of pcon sizes to read. Overwrites n_con from simulation
    only_pcon: bool
        if file is vel_sc, skip u,v,w,T and read just pcon. Makes reading pcon a lot faster.
    """
    from os import path
    from . import Simulation
    import numpy as np

    if isinstance(simulation, str):
        from ..simClass import Simulation as Sim
        simulation = Sim(simulation)

    #---------
    # Trying to be general
    if simulation!=None:
        sim = simulation
        domain = sim.domain
    else:
        if domain==None:
            sim = Simulation(path.dirname(path.abspath(fname)))
            domain = sim.domain
        else:
            sim = Simulation(domain=domain, n_con=n_con)
    #---------

    #---------
    # Useful for later. Might as well do it just once here
    if type(nz)==type(None):
        nz=domain.nz_tot-1
    if type(nz_full)==type(None):
        nz_full=domain.nz_tot-1
    
    kfill = nz_full - nz
    u_nd = domain.nx*domain.ny*nz
    u_fill = domain.nx*domain.ny*kfill
    u_nd2 = domain.ld*domain.ny
    #---------

    #---------
    if dtype==np.float64:
        byteprec = 8
    elif dtype==np.float32:
        byteprec = 4
    #---------

    #---------
    bfile = open(fname, 'rb')
    #---------

    if path.basename(fname).startswith('uv0_jt'):
        u0, v0 = np.fromfile(bfile, dtype=np.float64).reshape(2,-1, order='C')
        u0 = u0.reshape((domain.nx, domain.ny), order='F')
        v0 = v0.reshape((domain.nx, domain.ny), order='F')
        if as_DA:
            u0=sim.DataArray(u0*sim.u_scale, dims=['x', 'y'])
            v0=sim.DataArray(v0*sim.u_scale, dims=['x', 'y'])
        outlist = [u0, v0]

    #--------------
    # For fortran unformatted direct files you have to skip first 4 bytes
    if fname.endswith('.bin'): bfile.read(4)
    #--------------
    
   
    #---------
    # Straightforward
    elif path.basename(fname).startswith('theta_') or path.basename(fname).startswith('temp_tt'):
        T = np.fromfile(bfile, dtype=np.float64, count=u_nd).reshape((domain.nx, domain.ny, nz), order='F')

        #----
        # Temperature depends on environment
        if sim.ocean_flag:
            T = 2.*sim.t_init - T*sim.t_scale
        else:
            T = T*sim.t_scale
        #----

        #----
        if as_DA:
            T=sim.DataArray(T, dims=['x', 'y', 'z_u'])
            T.attrs=dict(long_name="Temperature", units="K")
        outlist = [T]
        #----
    #---------

    #---------
    # Straightforward
    elif path.basename(fname).startswith('uvw_') or path.basename(fname).startswith('vel_tt'):
        u = np.fromfile(bfile, dtype=dtype, count=u_nd).reshape((domain.nx, domain.ny, nz), order='F')
        bfile.seek(byteprec*(u_nd+u_fill))
        v = np.fromfile(bfile, dtype=dtype, count=u_nd).reshape((domain.nx, domain.ny, nz), order='F')
        bfile.seek(2*byteprec*(u_nd+u_fill))
        w = np.fromfile(bfile, dtype=dtype, count=u_nd).reshape((domain.nx, domain.ny, nz), order='F')
        u = u*sim.u_scale
        v = v*sim.u_scale
        w = w*sim.u_scale
        if sim.ocean_flag:
            w = -w
        if as_DA:
            u=sim.DataArray(u, dims=['x', 'y', 'z_u'])
            v=sim.DataArray(v, dims=['x', 'y', 'z_u'])
            w=sim.DataArray(w, dims=['x', 'y', 'z_w'])

            u.attrs=dict(long_name="$u$", units="m/s")
            v.attrs=dict(long_name="$v$", units="m/s")
            w.attrs=dict(long_name="$w$", units="m/s")

        outlist = [u,v,w]
    #---------

    #---------
    # Straightforward
    elif path.basename(fname).startswith('nu_'):
        nu = np.fromfile(bfile, dtype=np.float64, count=u_nd).reshape((domain.nx, domain.ny, nz), order='F')
        nu = nu * (sim.u_scale*sim.z_i)

        #----
        if as_DA:
            nu=sim.DataArray(nu, dims=['x', 'y', 'z_u'])
            nu.attrs=dict(long_name="SGS viscosity", units="m**2/s")
        outlist = [nu]
        #----
    #---------

    #---------
    # Straightforward
    if path.basename(fname).startswith('pcon_jt'):
        if n_con==None:
            n_con = sim.n_con
        pcon = []
        for i in range(n_con):
            bfile.seek(byteprec*i*(u_nd+u_fill))
            pcon.append(np.fromfile(bfile, dtype=dtype, count=u_nd).reshape((domain.nx, domain.ny, nz), order='F'))
        pcon = np.stack(pcon, axis=-1)
        pcon *= sim.pcon_scale
        if as_DA:
            pcon=sim.DataArray(pcon, dims=['x', 'y', 'z_u', pcon_index], n_con=n_con)
            pcon.attrs=dict(long_name="$C$", units="kg/m$^3$")
        outlist = [pcon]
    #---------



    #---------
    # Straightforward
    if path.basename(fname).startswith('wc_jt'):
        if n_con==None:
            n_con = sim.n_con
        wc = []
        for i in range(n_con):
            bfile.seek(byteprec*i*(u_nd+u_fill))
            wc.append(np.fromfile(bfile, dtype=dtype, count=u_nd).reshape((domain.nx, domain.ny, nz), order='F'))
        wc = np.stack(wc, axis=-1)
        wc *= sim.pcon_scale * sim.u_scale
        if as_DA:
            wc=sim.DataArray(wc, dims=['x', 'y', 'z_w', pcon_index])
            wc.attrs=dict(long_name="$wc$", units="kg/m$^3$ * m/s")
        outlist = [wc]
    #---------


    #---------
    # Spectra
    elif path.basename(fname).startswith('spec_uvwT'):
        ndz=4*(domain.nx//2)*(domain.nz_tot-1)
        u,v,w,T = np.fromfile(bfile, dtype=np.float32, count=ndz).reshape((4,domain.nx//2, domain.nz_tot-1), order='F')
        print('Not normalized')
        outlist = [u,v,w,T]
    #---------

    #---------
    # Add units as attributes
    if as_DA:
        from .utils import add_units
        for i in range(len(outlist)):
            outlist[i] = add_units(outlist[i])
    if len(outlist)==1:
        outlist=outlist[0]
    #---------

    return outlist


def fortran2xr(fname, vlist, padded=False, dtype=np.float64, verbose=False):
    """
    Reads a binary file according to the simulation or domain object passed

    Parameters
    ----------
    fname: string
        path of the binary file you want to open
    vlist: list 
        list of dictionaries with relevant parameters for each variable (see loop)
    """
    import xarray as xr

    #---------
    # Open file as binary
    bfile = open(fname, 'rb')
    #---------

    #--------------
    # For some fortran unformatted direct files you have to skip first 4 bytes
    if padded: bfile.read(4)
    #--------------
    
    #---------
    # Read each variable into a DataArray and put it in a list.
    # Each variable is assumed to be in sequence, unless `jumpto`
    # keyword is added to the dict. That is the byte count
    # at which each variable starts in the file. `jumpto` is
    # useful if you don't wanna read the whole thing and know
    # exactly where your variable is.
    dalist = []
    for i, vdict in enumerate(vlist):
        if verbose: print(i)

        #------
        # Shape can be infered from the coords. dtype can't. Both are needed to
        # calculate the amount of bytes to read
        shape = [ len(el) for el in vdict["coords"].values() ]
        if "dtype" not in vdict.keys(): 
            idtype=dtype
        else:
            idtype=vdict["dtype"]
        #------

        #------
        # Jump to given byte
        if "jumpto" in vdict:
            bfile.seek(vdict["jumpto"])
        #------

        #------
        # Read based on size of var, make it into DataArray and append
        arr = np.fromfile(bfile, dtype=idtype, count=np.prod(shape)).reshape(shape, order='F')
        dalist.append(xr.DataArray(arr, dims=vdict["dims"], coords=vdict["coords"]))
        #------

    bfile.close()
    return dalist


def readBinary2(fname, simulation=None, n_con=None, pcon_index='w_r',
                nz=None, nz_full=None):
    """
    Reads a binary file according to the simulation or domain object passed

    Parameters
    ----------
    fname: string
        path of the binary file you want to open
    simulation: lespy.Simulation
        simulation that contains important information to assemble the file. Mostly
        it's used to get number of points, if temperature is in the file or not and n_con.
    n_con: int
        number of pcon sizes to read. Overwrites n_con from simulation
    """
    from os import path
    from . import Simulation
    import numpy as np

    #---------
    # Trying to be general
    if simulation!=None:
        sim = simulation
    else:
        sim = Simulation(path.dirname(path.abspath(fname)))
    domain = sim.domain
    #---------

    #---------
    # Useful for later. Might as well do it just once here
    if type(nz)==type(None):
        nz=domain.nz_tot-1
    if type(nz_full)==type(None):
        nz_full=domain.nz_tot-1
    
    kfill = nz_full - nz
    u_nd = domain.nx*domain.ny*nz
    u_fill = domain.nx*domain.ny*kfill
    #---------

    #--------------
    # For fortran unformatted direct files you have to skip first 4 bytes
    if fname.endswith('.bin'): 
        padded=True
    else:
        padded=False
    #--------------

    #--------------
    from .utils import get_dicts
    vdict = get_dicts(sim, nz=nz)
    #--------------

    #---------
    # Straightforward
    if path.basename(fname).startswith('pcon'):
        vlist = [vdict["c"]]
        pcon = fortran2xr(fname, vlist, padded=padded, dtype=np.float64)
        pcon *= sim.pcon_scale
        outlist = [pcon]
    #---------
    
    #---------
    # Straightforward
    elif path.basename(fname).startswith('theta_jt') or path.basename(fname).startswith('temp_tt'):
        vlist = [vdict["θ"]]
        T, = fortran2xr(fname, vlist, padded=padded, dtype=np.float64)
        T = 2.*sim.t_init - T*sim.t_scale
        T.attrs=dict(long_name="Pot temperature", units="K")
        outlist = [T]
    #---------

    #---------
    # Straightforward
    elif path.basename(fname).startswith('uvw_jt') or path.basename(fname).startswith('vel_tt'):
        vdict["v"]["jumpto"] = sim.nx*sim.ny*nz_full*8
        vdict["w"]["jumpto"] = sim.nx*sim.ny*nz_full*8*2
        vlist = [vdict["u"], vdict["v"], vdict["w"],] 
        u,v,w = fortran2xr(fname, vlist, padded=padded, dtype=np.float64)
        u = u*sim.u_scale; v = v*sim.u_scale; w = w*sim.u_scale
        if sim.ocean_flag:
            w = -w

        u.attrs=dict(long_name="$u$", units="m/s")
        v.attrs=dict(long_name="$v$", units="m/s")
        w.attrs=dict(long_name="$w$", units="m/s")

        outlist = [u,v,w]
    #---------

    #---------
    # Spectra
    elif path.basename(fname).startswith('spec_uvwT'):
        ndz=4*(domain.nx//2)*(domain.nz_tot-1)
        u,v,w,T = np.fromfile(bfile, dtype=np.float32, count=ndz).reshape((4,domain.nx//2, domain.nz_tot-1), order='F')
        print('Not normalized')
        outlist = [u,v,w,T]
    #---------

    #---------
    # Add units as attributes
    from .utils import add_units
    for i in range(len(outlist)):
        outlist[i] = add_units(outlist[i])
    if len(outlist)==1:
        outlist=outlist[0]
    #---------

    return outlist


def readall_aver(outdir, sim=None, pcon_index="index", verbose=False, 
                 unique="ndtime", **kwargs):
    """
    Reads (almost) all output in the aver_*.out files.
    Simulation argument is mandatory
    """
    pci = pcon_index
    read_aver2 = lambda f, **kwargs: read_aver(f, simulation=sim, unique=unique, **kwargs)

    if verbose: print("Reading averaged fluxes")
    #------
    # θ fluxes
    wθ_res = read_aver2(outdir+"/aver_wt.out", 
                           dims=["ndtime", "z_w"])*sim.u_scale*sim.t_scale
    wθ_sgs = read_aver2(outdir+"/aver_sgs_t3.out", 
                           dims=["ndtime", "z_w"])*sim.u_scale*sim.t_scale
    #------

    #------
    # Covariances
    uw_res = read_aver2(outdir+"/aver_uw.out", dims=["ndtime", "z_w"])*sim.u_scale**2
    vw_res = read_aver2(outdir+"/aver_vw.out", dims=["ndtime", "z_w"])*sim.u_scale**2

    u2_res = read_aver2(outdir+"/aver_u2.out", dims=["ndtime", "z_u"])*sim.u_scale**2
    v2_res = read_aver2(outdir+"/aver_v2.out", dims=["ndtime", "z_u"])*sim.u_scale**2
    w2_res = read_aver2(outdir+"/aver_w2.out", dims=["ndtime", "z_w"])*sim.u_scale**2

    uw_sgs = read_aver2(outdir+"/aver_txz.out", dims=["ndtime", "z_w"])*sim.u_scale**2
    vw_sgs = read_aver2(outdir+"/aver_tyz.out", dims=["ndtime", "z_w"])*sim.u_scale**2

    u2_sgs = read_aver2(outdir+"/aver_txx.out", dims=["ndtime", "z_u"])*sim.u_scale**2
    v2_sgs = read_aver2(outdir+"/aver_tyy.out", dims=["ndtime", "z_u"])*sim.u_scale**2
    w2_sgs = read_aver2(outdir+"/aver_tzz.out", dims=["ndtime", "z_w"])*sim.u_scale**2
    #------

    #------
    # Change sign if it's ocean
    if sim.ocean_flag:
        uw_res = -uw_res; uw_sgs = -uw_sgs
        vw_res = -vw_res; vw_sgs = -vw_sgs
    #------

    #------
    # Scalar fluxes
    if sim.pcon_flag:
        wc_res = read_aver2(outdir+"/aver_wPCon.out",
                              dims=["ndtime", "z_w", pci])*sim.u_scale*sim.pcon_scale
        wc_sgs = read_aver2(outdir+"/aver_sgs_PCon3.out",
                              dims=["ndtime", "z_w", pci])*sim.u_scale*sim.pcon_scale
        if sim.ocean_flag:
            wc_res = -wc_res; wc_sgs = -wc_sgs
    #------


    if verbose: print("Reading averaged profiles")
    #------
    # Averaged profiles
    u = read_aver2(outdir+"/aver_u.out", dims=["ndtime", "z_u",])*sim.u_scale
    v = read_aver2(outdir+"/aver_v.out", dims=["ndtime", "z_u",])*sim.u_scale
    w = read_aver2(outdir+"/aver_w.out", dims=["ndtime", "z_w",])*sim.u_scale
    θ_nd = read_aver2(outdir+"/aver_theta.out", dims=["ndtime", "z_u"])
    if sim.ocean_flag:
        θ = 2.*sim.t_init - θ_nd*sim.t_scale
    else:
        θ = θ_nd * sim.t_scale

    if sim.pcon_flag:
        C = read_aver2(outdir+"/aver_PCon.out", dims=["ndtime", "z_u", pci])*sim.pcon_scale
    #------

    if verbose: print("Reading averaged derivatives")

    #------
    # Derivatives
#    dudz = read_aver2(outdir+"/aver_dudz.out", dims=["ndtime", "z_w"])*sim.u_scale/sim.z_i
#    dvdz = read_aver2(outdir+"/aver_dvdz.out", dims=["ndtime", "z_w"])*sim.u_scale/sim.z_i
#    dθdz = read_aver2(outdir+"/aver_dTdz.out", dims=["ndtime", "z_w"])*sim.t_scale/sim.z_i
#    if sim.pcon_flag:
#        dC0dz = read_aver2(outdir+"/aver_dPCondz.out", dims=["ndtime", "z_w", pci])*sim.pcon_scale/sim.z_i
    #------


    if verbose: print("Reading dissipation")
    #------
    # Dissipation
    ε = read_aver2(outdir+"/aver_dissip.out", dims=["ndtime", "z_u",])*sim.u_scale**3/sim.z_i
    #------

    #------
    uw = uw_res + uw_sgs
    vw = vw_res + vw_sgs
    wθ = wθ_res + wθ_sgs

    u2 = u2_res + u2_sgs
    v2 = v2_res + v2_sgs
    w2 = w2_res + w2_sgs
    if sim.pcon_flag:
        wc = wc_res + wc_sgs
    #------

    #------
    # Prepare output
    dict_unode = dict(ε=ε,
                      u=u, v=v, θ=θ, 
                      u2=u2, v2=v2,
                      u2_res=u2_res, v2_res=v2_res,
                      u2_sgs=u2_sgs, v2_sgs=v2_sgs,
                      )
    dict_wnode = dict(w=w, 
                      uw=uw, vw=vw, w2=w2, wθ=wθ, 
                      uw_res=uw_res, vw_res=vw_res, w2_res=w2_res, wθ_res=wθ_res, 
                      uw_sgs=uw_sgs, vw_sgs=vw_sgs, w2_sgs=w2_sgs, wθ_sgs=wθ_sgs, 
#                      dudz=dudz, dvdz=dvdz, dθdz=dθdz,
                      )
    if sim.pcon_flag:
        dict_wnode["wc_res"] = wc_res
        dict_wnode["wc_sgs"] = wc_sgs
        dict_wnode["wc"] = wc
#        dict_wnode["dC0dz"] = dC0dz
        dict_unode["C"] = C
    ds_unode = xr.Dataset(dict_unode)
    ds_wnode = xr.Dataset(dict_wnode)
    #------

    return ds_unode, ds_wnode
