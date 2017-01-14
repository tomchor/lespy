#from .decorators import autoassign

def readMeanPP(fpath):
    """Open post processing averaged file produced by postproc-avgs.f90
    """
    import pandas as pd
    avgs = pd.read_csv(fpath, index_col=None, header=None, delim_whitespace=True)

    columns = ['z_uv', 'z_w', 'z_cs', '<U>/u*', '<V>/u*', '<W>/u*', '<dUdz>/u**dz*',
                '<dVdz>/u**dz*', '<U>', '<V>', '<Theta>', '<u^2>/u*^2', '<v^2>/u*^2',
                '<w^2>/u*^2', '<uw>/u*^2', '<vw>/u*^2', '<wT>/u*T*', 'cs',
                'beta', 'beta_clip', 'cs_rns', '<txz>/u*^2', '<tyz>/u*^2']
    avgs.columns=columns
    return avgs


def postProcess2D(model_outputdir, t_ini=100000, t_end=None, simulation=None, return_df=False):
    """
    Postprocess average results from LES output
    Compiled for python from the fpostproc.f90 subroutine file
    """
    from .fpostproc import postproc
    from ..utils import paramParser
    from .. import simClass
    import numpy as np
    if not simulation:
        from os import path
        simulation = simClass.Simulation(path.join(model_outputdir, 'codebkp/param.nml'))

    nz = simulation.domain.nz
    Lz = simulation.domain.lz
    u_star = simulation.u_scale
    dt_dim = simulation.dt
    if not t_end:
        t_end = simulation.timelength
    T_scale = simulation.T_scale
    nt = int(simulation.timelength/simulation.avglength)

    outdat = postproc(model_outputdir, nz, Lz, u_star, dt_dim, nt, t_ini, t_end, T_scale)

    if return_df:
        import pandas as pd
        outdat = pd.DataFrame(outdat.T)
        columns = ['z_uv', 'z_w', 'z_cs', '<U>/u*', '<V>/u*', '<W>/u*', '<dUdz>/u**dz*',
                '<dVdz>/u**dz*', '<U>', '<V>', '<Theta>', '<u^2>/u*^2', '<v^2>/u*^2',
                '<w^2>/u*^2', '<uw>/u*^2', '<vw>/u*^2', '<wT>/u*T*', 'cs',
                'beta', 'beta_clip', 'cs_rns', '<txz>/u*^2', '<tyz>/u*^2']
        outdat.columns = columns

    return outdat


def readBinary3(fname, simulation=None, domain=None, engine='fortran', n_con=None):
    """
    Reads a binary file according to the simulation or domain object passed

    Passing a simulation might not be trustworthy because it might refer to different files
    """
    from os import path
    from .. import Simulation
    import numpy as np
    from struct import unpack, error

    if isinstance(simulation, str):
        from ..simClass import Simulation as Sim
        simulation = Sim(simulation)

    #---------
    # Dealing with pre-requesites. Using only domain is more general
    if (simulation==None) and (domain==None):
        domain = Simulation(path.dirname(path.abspath(fname))).domain
        sim=None
    elif (simulation==None) and (domain!=None):
        sim=None
    elif (simulation!=None) and (domain==None):
        sim = simulation
        domain = sim.domain
    else:
        sim = simulation
    #---------

    #--------------
    # For fortran unformatted output you have to skip first 4 bytes
    bfile = open(fname, 'rb')
    bfile.read(4)
    #--------------
    
    if path.basename(fname).startswith('con_tt'):
        p_nd = domain.ld*domain.ny*domain.nz_tot*sim.n_con
        pcon = unpack('d'*p_nd, bfile.read(8*p_nd))
        pcon = np.array(pcon).reshape((sim.domain.ld, sim.ny, sim.nz_tot, sim.n_con), order='F')
        return pcon
    
    elif path.basename(fname).startswith('vel_sc'):
        u_nd = domain.ld*domain.ny*domain.nz_tot
        u = np.fromfile(bfile, dtype=np.float64, count=u_nd).reshape((domain.ld, domain.ny, domain.nz_tot), order='F')
        v = np.fromfile(bfile, dtype=np.float64, count=u_nd).reshape((domain.ld, domain.ny, domain.nz_tot), order='F')
        w = np.fromfile(bfile, dtype=np.float64, count=u_nd).reshape((domain.ld, domain.ny, domain.nz_tot), order='F')
        T = np.fromfile(bfile, dtype=np.float64, count=u_nd).reshape((domain.ld, domain.ny, domain.nz_tot), order='F')
        return u,v,w, T
        #--------

    return
    #---------
    

def readBinary2(fname, simulation=None, domain=None, read_pcon=True, n_con=None, read_just_pcon=False):
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
    read_just_pcon: bool
        if file is vel_sc, skip u,v,w,T and read just pcon. Makes reading pcon a lot faster.
    """
    from os import path
    from .. import Simulation
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
    u,v,w,T,pcon = [None]*5
    u_nd = domain.ld*domain.ny*domain.nz_tot
    if read_pcon or read_just_pcon:
        if n_con==None:
            n_con = sim.n_con
        p_nd = u_nd*n_con
    #---------

    #--------------
    # For fortran unformatted output you have to skip first 4 bytes
    bfile = open(fname, 'rb')
    bfile.read(4)
    #--------------
    
    #---------
    # Straightforward
    if path.basename(fname).startswith('con_tt'):
        pcon = np.fromfile(bfile, dtype=np.float64, count=u_nd).reshape((domain.ld, domain.ny, domain.nz_tot), order='F')
    #---------
    
    #---------
    # Straightforward
    elif path.basename(fname).startswith('temp_t'):
        T = np.fromfile(bfile, dtype=np.float64, count=u_nd).reshape((domain.ld, domain.ny, domain.nz_tot), order='F')
    #---------

    #---------
    # Straightforward
    elif path.basename(fname).startswith('vel_t'):
        u = np.fromfile(bfile, dtype=np.float64, count=u_nd).reshape((domain.ld, domain.ny, domain.nz_tot), order='F')
        v = np.fromfile(bfile, dtype=np.float64, count=u_nd).reshape((domain.ld, domain.ny, domain.nz_tot), order='F')
        w = np.fromfile(bfile, dtype=np.float64, count=u_nd).reshape((domain.ld, domain.ny, domain.nz_tot), order='F')
    #---------

    #---------
    # Here lies the problem! The same file name for a bunch of formats! (should consider change)
    elif path.basename(fname).startswith('vel_sc'):
        #---------
        # If you want just the concentration, this will skip everything else and will read faster
        if read_just_pcon:
            if sim.s_flag:
                bfile.read(8*u_nd*4)
            else:
                bfile.read(8*u_nd*3)
            pcon = np.fromfile(bfile, dtype=np.float64, count=p_nd).reshape((domain.ld, domain.ny, domain.nz_tot, sim.n_con), order='F')
            return pcon*sim.pcon_scale
        #---------

        #---------
        # Every vel_sc file has at least u,v,w, so we start with that
        u = np.fromfile(bfile, dtype=np.float64, count=u_nd).reshape((domain.ld, domain.ny, domain.nz_tot), order='F')
        v = np.fromfile(bfile, dtype=np.float64, count=u_nd).reshape((domain.ld, domain.ny, domain.nz_tot), order='F')
        w = np.fromfile(bfile, dtype=np.float64, count=u_nd).reshape((domain.ld, domain.ny, domain.nz_tot), order='F')
        #---------

        #---------
        # If the scalar flag is on, means temperature is output
        # It's very unlikely that this flag is on in the beginning of a simnulation and
        # off in other parts, so we won't try to correct for that
        if sim.s_flag:
            T = np.fromfile(bfile, dtype=np.float64, count=u_nd).reshape((domain.ld, domain.ny, domain.nz_tot), order='F')
        #---------

        #---------
        # Finally, if the pcon_flag is on, we try to read the pcon field. This flag might be off
        # on the beginning of the vel_sc outputs and then on when pcon is released, so we do it as a try clause
        if read_pcon:
            try:
                pcon = np.fromfile(bfile, dtype=np.float64, count=p_nd).reshape((domain.ld, domain.ny, domain.nz_tot, sim.n_con), order='F')
            except ValueError as e:
                if n_con!=None:
                    raise e
                else:
                    print("Couldn't read pcon value at the end. Life goes on.")
        #---------
    #---------

    #---------
    # Now we set up the output and re-scale one by one of they exist
    outlist = []
    if isinstance(u, np.ndarray):
        u *= sim.u_scale
        v *= sim.u_scale
        w *= sim.u_scale
        outlist+=[u,v,w]
    if isinstance(T, np.ndarray):
        T = 2.*sim.t_init - T*sim.t_scale
        outlist.append(T)
    if isinstance(pcon, np.ndarray):
        pcon *= sim.pcon_scale
        outlist.append(pcon)
    #---------

    #---------
    # We simplify if there's only one output
    if len(outlist)==1:
        return outlist[0]
    else:
        return outlist
    #---------
    

def readBinary(fname, simulation=None, domain=None, engine='fortran'):
    """Reads a binary file according to the simulation object passed
    """
    from os import path
    from .. import Simulation

    if isinstance(simulation, str):
        from ..simClass import simulation as Sim
        simulation = Sim(simulation)

    #---------
    # Dealing with pre-requesites. Using only domain is more general
    if (simulation==None) and (domain==None):
        domain = Simulation(fname).domain
    elif (simulation==None) and (domain!=None):
        pass
    elif (simulation!=None) and (domain==None):
        sim = simulation
        domain = sim.domain
    else:
        sim = simulation
    #---------

    #---------
    # If reading with python, it's more flexible but around 20 times slower
    if engine=='python':
        import numpy as np
        from struct import unpack
        #--------------
        # For fortran unformatted output you have to skip first 4 bytes
        bfile = open(fname, 'rb')
        bfile.read(4)
        #--------------
    
        if path.basename(fname).startswith('con_tt'):
            p_nd = sim.domain.ld*sim.ny*sim.nz_tot*sim.n_con
            pcon = unpack('d'*p_nd, bfile.read(8*p_nd))
            pcon = np.array(pcon).reshape((sim.domain.ld, sim.ny, sim.nz_tot, sim.n_con), order='F')
            return pcon
    
        elif path.basename(fname).startswith('vel_sc'):
            u_nd = domain.ld*domain.ny*domain.nz_tot
            u = np.array(unpack('d'*u_nd, bfile.read(8*u_nd))).reshape((domain.ld, domain.ny, domain.nz_tot), order='F')
            v = np.array(unpack('d'*u_nd, bfile.read(8*u_nd))).reshape((domain.ld, domain.ny, domain.nz_tot), order='F')
            w = np.array(unpack('d'*u_nd, bfile.read(8*u_nd))).reshape((domain.ld, domain.ny, domain.nz_tot), order='F')

            if simulation.s_flag and simulation.pcon_flag:
                p_nd = sim.domain.ld*sim.ny*sim.nz_tot*sim.n_con
                T = np.array(unpack('d'*u_nd, bfile.read(8*u_nd))).reshape((sim.domain.ld, sim.ny, sim.nz_tot), order='F')
                pcon = np.array(unpack('d'*p_nd, bfile.read(8*p_nd))).reshape((sim.domain.ld, sim.ny, sim.nz_tot, sim.n_con), order='F')
                return u,v,w,T,pcon

            elif simulation.s_flag and (not simulation.pcon_flag):
                T = np.array(unpack('d'*u_nd, bfile.read(8*u_nd))).reshape((sim.domain.ld, sim.ny, sim.nz_tot), order='F')
                return u,v,w,T
    
            elif (not simulation.s_flag) and (not simulation.pcon_flag):
                return u,v,w

            elif (not simulation.s_flag) and simulation.pcon_flag:
                p_nd = sim.domain.ld*sim.ny*sim.nz_tot*sim.n_con
                pcon = np.array(unpack('d'*p_nd, bfile.read(8*p_nd))).reshape((sim.domain.ld, sim.ny, sim.nz_tot, sim.n_con), order='F')
                return u,v,w,pcon
        return
    #---------
    
    #---------
    # Reading with fortran is less flexible but 20 times faster
    elif engine=='fortran':
        from . import read_instant3 as read

        if path.basename(fname).startswith('con_tt'):
            u,v,w,T,pcon = read.read_binary(fname, sim.nx, sim.ny, sim.nz_tot, sim.n_con, sim.s_flag, sim.pcon_flag, sim.flag_endless, 'con_tt')
            return pcon

        elif path.basename(fname).startswith('vel_sc'):
            u,v,w,T,pcon = read.read_binary(fname, sim.nx, sim.ny, sim.nz_tot, sim.n_con, sim.s_flag, sim.pcon_flag, sim.flag_endless, 'vel_sc')
            if simulation.s_flag and simulation.pcon_flag:
                return u,v,w,T,pcon
            if simulation.s_flag and (not simulation.pcon_flag):
                return u,v,w,T

        return
    #---------

