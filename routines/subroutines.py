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


    
def readBinary2(fname, simulation=None, domain=None, read_pcon=True, n_con=None, only_pcon=False, 
        trim=True, as_dataarray=False):
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
    u_nd2 = domain.ld*domain.ny
    if read_pcon or only_pcon:
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
    # Spectra
    elif path.basename(fname).startswith('spec_uvwT'):
        ndz=4*(domain.nx//2)*(domain.nz_tot-1)
        u,v,w,T = np.fromfile(bfile, dtype=np.float32, count=ndz).reshape((4,domain.nx//2, domain.nz_tot-1), order='F')
        return u,v,w,T
    #---------


    #---------
    # Straightforward
    elif path.basename(fname).startswith('div_z0_t'):
        #-----
        # Note that we will multiply by u_star at the bottom of the function, so we only vidide by z_i here
        u = np.fromfile(bfile, dtype=np.float64, count=u_nd2).reshape((domain.ld, domain.ny), order='F')/sim.inversion_depth
        v = np.fromfile(bfile, dtype=np.float64, count=u_nd2).reshape((domain.ld, domain.ny), order='F')/sim.inversion_depth
        w =-np.fromfile(bfile, dtype=np.float64, count=u_nd2).reshape((domain.ld, domain.ny), order='F')/sim.inversion_depth
        #-----
    #---------

    elif path.basename(fname).startswith('vel_srf'):
        u = np.fromfile(bfile, dtype=np.float64, count=u_nd2).reshape((domain.ld, domain.ny), order='F')*sim.u_scale
        v = np.fromfile(bfile, dtype=np.float64, count=u_nd2).reshape((domain.ld, domain.ny), order='F')*sim.u_scale
        return u, v

 
    elif path.basename(fname).startswith('dvel_srf'):
        ux = np.fromfile(bfile, dtype=np.float64, count=u_nd2).reshape((domain.ld, domain.ny), order='F')*sim.u_scale/sim.inversion_depth
        uy = np.fromfile(bfile, dtype=np.float64, count=u_nd2).reshape((domain.ld, domain.ny), order='F')*sim.u_scale/sim.inversion_depth
        vx = np.fromfile(bfile, dtype=np.float64, count=u_nd2).reshape((domain.ld, domain.ny), order='F')*sim.u_scale/sim.inversion_depth
        vy = np.fromfile(bfile, dtype=np.float64, count=u_nd2).reshape((domain.ld, domain.ny), order='F')*sim.u_scale/sim.inversion_depth
        wz =-np.fromfile(bfile, dtype=np.float64, count=u_nd2).reshape((domain.ld, domain.ny), order='F')*sim.u_scale/sim.inversion_depth
        return ux, uy, vx, vy, wz

    #---------
    # Here lies the problem! The same file name for a bunch of formats! (should consider change)
    elif path.basename(fname).startswith('vel_sc'):
        #---------
        # If you want just the concentration, this will skip everything else and will read faster
        if only_pcon:
            if sim.s_flag:
                bfile.read(8*u_nd*4)
            else:
                bfile.read(8*u_nd*3)
            pcon = np.fromfile(bfile, dtype=np.float64, count=p_nd).reshape((domain.ld, domain.ny, domain.nz_tot, sim.n_con), order='F')
            if trim:
                pcon=pcon[:sim.nx]
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
                if n_con!=None and sim.pcon_flag:
                    print("Trying to read pcon failed. Set n_con=0 or read_pcon=0 to specify you don't want to read it")
                    raise e
                elif str(e).startswith('cannot reshape array of size 0'):
                    pass
                else:
                    print("Couldn't read pcon value at the end. Life goes on.")
        #---------
    #---------

    #---------
    # Now we set up the output and re-scale one by one of they exist
    outlist = []
    if isinstance(u, np.ndarray):
        u = u*sim.u_scale
        v = v*sim.u_scale
        w =-w*sim.u_scale
        outlist+=[u,v,w]
    if isinstance(T, np.ndarray):
        T = 2.*sim.t_init - T*sim.t_scale
        outlist.append(T)
    if isinstance(pcon, np.ndarray):
        pcon *= sim.pcon_scale
        outlist.append(pcon)
    #---------

    if trim:
        for i in range(len(outlist)):
            outlist[i] = outlist[i][:sim.nx]

    if as_dataarray:
        for i in range(len(outlist)):
            outlist[i] = sim.DataArray(outlist[i])

    #---------
    # We simplify if there's only one output
    if len(outlist)==1:
        return outlist[0]
    else:
        return outlist
    #---------
    

def interp_w(w, axis=-1):
    """ Interpolate w to uv_nodes keeping the same shape """
    import numpy as np
    w=np.asarray(w)
    w_int=(w.take(np.arange(0, w.shape[axis]-1), axis=axis) + w.take(np.arange(1, w.shape[axis]), axis=axis))/2
    pads=[ (0,0) for ax in range(w.ndim) ]
    pads[axis]=(0,1)
    w_int=np.pad(w_int, pads, 'constant', constant_values=0)
    return w_int


def monitor_stats(path, simulation=None, block=50, nblocks=6, Nstart=0, outname=None, use_deardorff=False, **kwargs):
    """
    Monitor important statistics and profiles of the simulation to see if it has converged or not
    """
    import numpy as np
    from matplotlib import pyplot as plt

    sim=simulation

    w_star=sim.w_star
    u_scale=sim.u_scale
    if use_deardorff:
        uw_scale=sim.u_scale/w_star
        t_conv=sim.inv_depth/w_star
    else:
        uw_scale=sim.u_scale/u_scale
        t_conv=np.inf
    t_star=sim.inv_depth/sim.u_scale
    count=sim.p_count
    Z=sim.domain.z_w[:-1]

    def fromtxt(fname, **kw_args):
        kw_args.update(kwargs)
        return np.genfromtxt(fname, **kw_args)
    
    #-----
    # Read wT
    print('Opening', path+'/aver_sgs_t3.out')
    sgs_t3=fromtxt(path+'/aver_sgs_t3.out', skip_header=Nstart//count)
    ndtimes=sgs_t3[:,0]
    if use_deardorff:
        sgs_t3=sgs_t3[:,1:]*sim.t_scale*sim.u_scale/sim.wt_s
    else:
        sgs_t3=sgs_t3[:,1:]
    
    print('Using', nblocks*block, 'lines of', len(sgs_t3), 'lines read from files.')

    print('Opening', path+'/aver_wt.out')
    if use_deardorff:
        res_t3=fromtxt(path+'/aver_wt.out', skip_header=Nstart//count)*sim.t_scale*sim.u_scale/sim.wt_s
    else:
        res_t3=fromtxt(path+'/aver_wt.out', skip_header=Nstart//count)
    res_t3=res_t3[:,1:]
    
    wT=res_t3+sgs_t3
    wT_mean=wT[-nblocks*block:].reshape(block,-1,len(Z)).mean(0)
    #zi_list = lp.physics.get_zi(wT, simulation=sim)
    #-----
    
    #-----
    # Read u2
    print('Opening', path+'/aver_u2.out')
    res_u2=fromtxt(path+'/aver_u2.out', skip_header=Nstart//count)*uw_scale**2
    res_u2=res_u2[:,1:]
    res_u2_mean=res_u2[-nblocks*block:].reshape(block,-1,len(Z)).mean(0)
    #-----
    
    #-----
    # Read w2
    print('Opening', path+'/aver_w2.out')
    res_w2=fromtxt(path+'/aver_w2.out', skip_header=Nstart//count)*uw_scale**2
    res_w2=res_w2[:,1:]
    res_w2_mean=res_w2[-nblocks*block:].reshape(block,-1,len(Z)).mean(0)
    #-----
    
    #-----
    # Read w3
    print('Opening', path+'/aver_w3.out')
    res_w3=-fromtxt(path+'/aver_w3.out', skip_header=Nstart//count)*uw_scale**3
    res_w3=res_w3[:,1:]
    res_w3_mean=res_w3[-nblocks*block:].reshape(block,-1,len(Z)).mean(0)
    #-----
    
    #-----
    # Read T2
    print('Opening', path+'/aver_var_t.out')
    res_T2=fromtxt(path+'/aver_var_t.out', skip_header=Nstart//count)*sim.t_scale*w_star/(sim.wt_s)
    res_T2=res_T2[:,1:]
    res_T2_mean=res_T2[-nblocks*block:].reshape(block,-1,len(Z)).mean(0)
    #-----
    
    if use_deardorff: 
        ndtimes=ndtimes*(t_star/t_conv)
    #itimes=ndtimes*t_star/sim.dt
    
    def adjust_fig(figure):
        for ax in figure.axes:
            ax.grid()
            ax.set_xticks(ax.get_xticks()[::2])
        figure.tight_layout()
    
    #-----
    # u^2
    fig, axes = plt.subplots(2,5, gridspec_kw=dict(height_ratios=[3,1]), figsize=(19,9))
    axes[0,0].set_title('$u^2$')
    axes[1,0].set_title('MAX$(u^2)$')
    #axes[0,0].plot(res_u2[-30:].mean(0), Z)
    for i, arr_i in enumerate(res_u2_mean):
        axes[0,0].plot(arr_i, Z, label='t{}'.format(i))
    axes[0,0].legend()
    axes[1,0].plot(ndtimes, res_u2.max(1))
    #-----
    
    if 1:
        #-----
        # SGS_t
        sgs_t3_mean=sgs_t3[-nblocks*block:].reshape(block,-1,len(Z)).mean(0)
        axes[0,1].set_title('SGS $wT$')
        axes[1,1].set_title('MAX(SGS $T^2)$')
        for i, arr_i in enumerate(sgs_t3_mean):
            axes[0,1].plot(arr_i, Z, label='t{}'.format(i))
        axes[0,1].legend()
        axes[1,1].plot(ndtimes, sgs_t3.max(1))
        #-----
    else:
        #-----
        # T^2
        axes[0,1].set_title('$T^2$')
        axes[1,1].set_title('MAX$(T^2)$')
        #axes[0,1].plot(res_T2[-30:].mean(0), Z)
        for i, arr_i in enumerate(res_T2_mean):
            axes[0,1].plot(arr_i, Z, label='t{}'.format(i))
        axes[0,1].legend()
        axes[1,1].plot(ndtimes, res_T2.max(1))
        #-----
    
    #-----
    # w^2
    axes[0,2].set_title('$w^2$')
    axes[1,2].set_title('MAX$(w^2)$')
    #axes[0,2].plot(res_w2[-30:].mean(0), Z)
    for i, arr_i in enumerate(res_w2_mean):
        axes[0,2].plot(arr_i, Z, label='t{}'.format(i))
    axes[0,2].legend()
    axes[1,2].plot(ndtimes, res_w2.max(1))
    #-----
    
    #-----
    # w^3
    axes[0,3].set_title('$w^3$')
    axes[1,3].set_title('MIN$(w^3)$')
    #axes[0,3].plot(res_w3[-30:].mean(0), Z)
    for i, arr_i in enumerate(res_w3_mean):
        axes[0,3].plot(arr_i, Z, label='t{}'.format(i))
    axes[0,3].legend()
    axes[1,3].plot(ndtimes, res_w3.min(1))
    #-----
    
    #-----
    # w*T
    axes[0,4].set_title('$wT$')
    axes[1,4].set_title('MIN$(wT)$')
    #axes[0,4].plot(wT[-30:].mean(0), Z)
    for i, arr_i in enumerate(wT_mean):
        axes[0,4].plot(arr_i, Z, label='t{}'.format(i))
    axes[0,4].legend()
    axes[1,4].plot(ndtimes, wT.min(1))
    adjust_fig(fig)
    fig.savefig(outname)
    #-----
    return axes


def combine_DAs(DAs):
    """
    Combines DataArrays based on their coordinates. Usefull for endless patches
    """
    out = DAs[0]
    for DA in DAs[1:]:
        out = out.combine_first(DA)
    return out

