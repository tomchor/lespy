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
        simulation = simClass.simulation(path.join(model_outputdir, 'codebkp/param.nml'))

    nz = simulation.domain.Nz
    Lz = simulation.domain.Lz
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


def readBinary(fname, simulation=None, engine='fortran'):
    """Reads a binary file according to the simulation object passed
    """
    from os import path

    if isinstance(simulation, str):
        from ..simClass import simulation as Sim
        simulation = Sim(simulation)

    sim = simulation

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
            p_nd = sim.domain.Ld*sim.ny*sim.nz_tot*sim.n_con
            pcon = unpack('d'*p_nd, bfile.read(8*p_nd))
            pcon = np.array(pcon).reshape((sim.domain.Ld, sim.ny, sim.nz_tot, sim.n_con), order='F')
            return pcon
    
        elif path.basename(fname).startswith('vel_sc'):
            u_nd = sim.domain.Ld*sim.ny*sim.nz_tot
            u = np.array(unpack('d'*u_nd, bfile.read(8*u_nd))).reshape((sim.domain.Ld, sim.ny, sim.nz_tot), order='F')
            v = np.array(unpack('d'*u_nd, bfile.read(8*u_nd))).reshape((sim.domain.Ld, sim.ny, sim.nz_tot), order='F')
            w = np.array(unpack('d'*u_nd, bfile.read(8*u_nd))).reshape((sim.domain.Ld, sim.ny, sim.nz_tot), order='F')

            if simulation.s_flag and simulation.pcon_flag:
                p_nd = sim.domain.Ld*sim.ny*sim.nz_tot*sim.n_con
                T = np.array(unpack('d'*u_nd, bfile.read(8*u_nd))).reshape((sim.domain.Ld, sim.ny, sim.nz_tot), order='F')
                pcon = np.array(unpack('d'*p_nd, bfile.read(8*p_nd))).reshape((sim.domain.Ld, sim.ny, sim.nz_tot, sim.n_con), order='F')
                return u,v,w,T,pcon

            elif simulation.s_flag and (not simulation.pcon_flag):
                T = np.array(unpack('d'*u_nd, bfile.read(8*u_nd))).reshape((sim.domain.Ld, sim.ny, sim.nz_tot), order='F')
                return u,v,w,T
    
            elif (not simulation.s_flag) and (not simulation.pcon_flag):
                return u,v,w

            elif (not simulation.s_flag) and simulation.pcon_flag:
                p_nd = sim.domain.Ld*sim.ny*sim.nz_tot*sim.n_con
                pcon = np.array(unpack('d'*p_nd, bfile.read(8*p_nd))).reshape((sim.domain.Ld, sim.ny, sim.nz_tot, sim.n_con), order='F')
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

