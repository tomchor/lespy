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

def postProcessAvgs(outputdir, t_ini=55000, t_end=None, simulation=None, postpfile='postp.out'):
    """
    Postprocess average results from LES output
    Compiled for python from the postproc-vel2.f90 routine
    """
    from .postmod import postproc
    if not simulation:
        from os import path
        params = paramParser(path.join(outputdir, 'codebkp/params.nml.bkp'))

    nz = simulation.domain.Nz
    Lz = simulation.domain.Lz
    u_star = simulation.u_scale
    dt_dim = simulation.dt
    if not t_end:
        t_end = simulation.timelength
    T_scale = simulation.T_scale
    nt = int(simulation.timelength/simulation.avglength)
    postproc(outputdir, postpfile, nz, Lz, u_star, dt_dim, nt, t_ini, t_end, T_scale)
    return


def postProcess2D(outputdir, t_ini=55000, t_end=None, simulation=None, postpfile='aux.txt'):
    """
    Postprocess average results from LES output
    Compiled for python from the postproc-vel2.f90 routine
    """
    from .fpostproc import postproc
    import numpy as np
    if not simulation:
        from os import path
        params = paramParser(path.join(outputdir, 'codebkp/params.nml.bkp'))

    nz = simulation.domain.Nz
    Lz = simulation.domain.Lz
    u_star = simulation.u_scale
    dt_dim = simulation.dt
    if not t_end:
        t_end = simulation.timelength
    T_scale = simulation.T_scale
    nt = int(simulation.timelength/simulation.avglength)
    aux = postproc(outputdir, nz, Lz, u_star, dt_dim, nt, t_ini, t_end, T_scale)

    return aux


