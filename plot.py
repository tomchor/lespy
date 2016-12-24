from matplotlib import pyplot as plt
from matplotlib import animation as anim

def check_avgs(outputs, t_ini, t_end, savefigs=False, return_df=True,
        normalize=False, plot_kwargs={}, means=True, variances=True, covariances=False):
    """
    Plots important averages from the 2D outputs

    Parameters
    ----------
    outputs: string
        output directory where the averages are.
    t_ini: int
        timestep to start the average.
    t_end: int
        timestep to end the average.
    savefigs: bool
        whether or not to save the figs. If false, figs are just shown.
    """
    from . import postProcess2D
    from . import Simulation

    sim = Simulation(outputs)
    dat = postProcess2D(outputs, t_ini=t_ini, t_end=t_end, simulation=sim, return_df=return_df)

    if normalize:
        dat.iloc[:, :3] = -dat.iloc[:, :3]/sim.inversion_depth
    else:
        dat.iloc[:, :3] = -dat.iloc[:, :3]

    if means:
        xlim=[1.2*dat['<U>/u*'].min(), 1.2*dat['<U>/u*'].max()]
        dat.plot(x='<U>/u*', y='z_uv', grid=True, xlim=xlim, **plot_kwargs)
        if savefigs: plt.savefig('U_{!s}-{!s}.png'.format(t_ini, t_end))
    
        xlim=[1.2*dat['<V>/u*'].min(), 1.2*dat['<V>/u*'].max()]
        dat.plot(x='<V>/u*', y='z_uv', grid=True, xlim=xlim, **plot_kwargs)
        if savefigs: plt.savefig('V_{!s}-{!s}.png'.format(t_ini, t_end))
    
    if variances:
        xlim=[1.2*dat['<u^2>/u*^2'].min(), 1.2*dat['<u^2>/u*^2'].max()]
        dat.plot(x='<u^2>/u*^2', y='z_uv', xlim=xlim, grid=True, **plot_kwargs)
        if savefigs: plt.savefig('uu_{!s}-{!s}.png'.format(t_ini, t_end))
    
        xlim=[1.2*dat['<v^2>/u*^2'].min(), 1.2*dat['<v^2>/u*^2'].max()]
        dat.plot(x='<v^2>/u*^2', y='z_uv', xlim=xlim, grid=True, **plot_kwargs)
        if savefigs: plt.savefig('vv_{!s}-{!s}.png'.format(t_ini, t_end))
        
        xlim=[1.2*dat['<w^2>/u*^2'].min(), 1.2*dat['<w^2>/u*^2'].max()]
        dat.plot(x='<w^2>/u*^2', y='z_uv', xlim=xlim, grid=True, **plot_kwargs)
        if savefigs: plt.savefig('ww_{!s}-{!s}.png'.format(t_ini, t_end))
    
    if covariances:
        xlim=[1.2*dat['<uw>/u*^2'].min(), 1.2*dat['<uw>/u*^2'].max()]
        dat.plot(x='<uw>/u*^2', y='z_uv', xlim=xlim, grid=True, **plot_kwargs)
        plt.savefig('uw_{!s}-{!s}.png'.format(t_ini, t_end))

        xlim=[1.2*dat['<vw>/u*^2'].min(), 1.2*dat['<vw>/u*^2'].max()]
        dat.plot(x='<vw>/u*^2', y='z_uv', xlim=xlim, grid=True, **plot_kwargs)
        plt.savefig('vw_{!s}-{!s}.png'.format(t_ini, t_end))

        xlim=[1.2*dat['<wT>/u*T*'].min(), 1.2*dat['<wT>/u*T*'].max()]
        dat.plot(x='<wT>/u*T*', y='z_uv', xlim=xlim, grid=True, **plot_kwargs)
        plt.savefig('wT_{!s}-{!s}.png'.format(t_ini, t_end))

    if not savefigs: plt.show()
    return dat



def pcon_2D_animation(results, outname=None, which='xy', simulation=None, title='',
        clim=[], logscale=True, axes=[], levels=None, timelist=None, aspect='equal', **kwargs):
    """
    Prints 2D animations from binary data for oil concentrations

    Each frame is done with the plane() function, and changes to that can be
    passed to plane() with the plane_kw keyword

    Parameters
    ----------
    resuts: np.array
        3D array where first axis is time
    outname: string
        path for animation to be saved. If None, animation is returned.
    which: string
        either xz or xy. If it's anything else no label will be put and x and y axis will be nodes.
    simulation: lespy.Simulation
        simulation from where to retrieve information
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as anim
    from .utils import get_ticks

    fig = plt.figure()
    snaps = []
    sim = simulation

    #--------
    # This sets the labels for each frame
    if timelist==None:
        timelist = range(0, results.shape[0])
    #--------

    #------------
    # Sets levels and ticks to use on the contour
    ticklabels, formatting, levels_con = get_ticks(results, levels=levels, logscale=logscale, clim=clim, nbins=6)
    #------------

    #-------
    # To take care of slightly negative or zero concentrations for plotting purposes
    results[ results==0 ] += 1.e-20
    results = abs(results)
    #-------

    for i, planecut in zip(timelist, results):
        print(i)
        aux = plane(planecut, which=which, axes=axes, outname=None, levels=levels_con, 
                simulation=simulation, logscale=logscale, set_cbar=False, **kwargs)
        timelabel = aux.ax.annotate(i, (.8,0.05), annotation_clip=False, xycoords='figure fraction')
        aux.collections.append(timelabel)
        snaps.append(aux.collections)
    fig.axes[0].set_title(title)

    #--------
    # Doesnt waste space on the plot and sets aspect ratio to be equal
    plt.gca().set_aspect(aspect, 'datalim')
    plt.tight_layout()
    #--------

    #--------
    # Animate and put the colorbar
    animated = anim.ArtistAnimation(fig, snaps, interval=50, blit=True)
    cbar = plt.colorbar(aux, ticks = levels_con, extend=None)
    cbar.ax.set_yticklabels(formatting(ticklabels))
    #--------

    if outname:
        animated.save(outname, bitrate=-1, dpi=120)
    else:
        return animated



def plane(plane, outname=None, which='xy', simulation=None, 
        axes=[], levels=None, logscale=True, clim=[],
        set_cbar=True, cmap=None, xlim=[], ylim=[]):
    """
    Prints 2D animations from binary data for oil concentrations
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from .utils import get_ticks

    sim = simulation
    if cmap:
        cmap_con = cm.get_cmap(cmap)
    else:
        cmap_con = cm.get_cmap("jet")

    #------------
    # Sets levels and ticks to use on the contour
    ticklabels, formatting, levels_con = get_ticks(plane, levels=levels, logscale=logscale, clim=clim, nbins=6)
    #------------

    #------------
    # Sets x y grid
    if axes:
        ax1, ax2 = axes
    else:
        if simulation:
            if which == 'xy':
                ax1 = np.arange(0, plane.shape[0]*sim.domain.dx, sim.domain.dx)
                ax2 = np.arange(0, plane.shape[1]*sim.domain.dy, sim.domain.dy)
            elif which == 'xz':
                ax1 = np.arange(0, plane.shape[0]*sim.domain.dx, sim.domain.dx)
                ax2 = -np.arange(0, plane.shape[1]*sim.domain.dz, sim.domain.dz)
        else:
            ax1 = np.arange(0, plane.shape[0])
            ax2 = np.arange(0, plane.shape[1])
    X2, X1 = np.meshgrid(ax2, ax1)
    #------------

    #----------
    # Sets labels latex style (no label is set if which is not given
    if which=='xz':
        plt.xlabel('$\mathbf{x(m)}$')
        plt.ylabel('$\mathbf{z(m)}$')
    elif which=='xy':
        plt.xlabel('$\mathbf{x(m)}$')
        plt.ylabel('$\mathbf{y(m)}$')
    #----------

    #----------
    # To plot logscale, we must do it manually
    if logscale:
        plane_log = np.log10(plane)
    else:
        plane_log = plane
    #----------

    #-------
    # Actual plotting is done here
    aux = plt.contourf(X1, X2, plane_log, levels_con, cmap=cmap_con, extend='both')
    aux.ax.set_xlim(*xlim)
    aux.ax.set_ylim(*ylim)
    #-------

    #-------
    # If this is going to animation, cbar shoudnt be set here
    if set_cbar:
        cbar = plt.colorbar(aux, ticks = ticklabels, extend='both')
        cbar.ax.set_yticklabels(formatting(ticklabels))
    #-------

    if outname:
        return plt.savefig(outname, bbox_inches='tight')
    else:
        return aux



def pcon_side_animation(bins, outname=None, simulation=None,
        n_pcon=0, xz_func=lambda x: x.mean(axis=1),
        xy_func=lambda x: x.mean(axis=2), trim_x=True):
    """
    Prints 2D animations from binary data for oil concentrations
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as anim
    from matplotlib import colors
    from . import routines

    #fig = plt.figure()
    fig, axes = plt.subplots(2, sharex=True, figsize=(6,12))
    fig.tight_layout()
    axes[0].set_aspect('equal')
    axes[1].set_aspect('equal')
    snaps = []

    sim = simulation
    if trim_x:
        xar = np.arange(0,sim.domain.Nx)*sim.domain.dx
    else:
        xar = np.arange(0,sim.domain.Ld)*sim.domain.dx
    yar = (np.arange(0,sim.ny)*sim.domain.dy).reshape(-1,1)
    zar = (np.arange(0,sim.nz_tot)*sim.domain.dz).reshape(-1,1)


    for fname in bins:
        print(fname)
        u,v,w,T,pcon = routines.readBinary(fname, simulation=sim)
        if trim_x:
            pcon = pcon[:-2,:,:, n_pcon]
        else:
            pcon = pcon[:,:,:, n_pcon]

        #-------
        # To take care of slightly negative or zero concentrations for plotting purposes
        pcon[ pcon==0 ] += 1.e-20
        pcon = abs(pcon)
        #-------

        flat1 = xz_func(pcon)
        axes[1].set_xlabel('x')
        axes[1].set_ylabel('z')
        a2 = -zar
        a1 = xar
        im1 = axes[1].pcolormesh(a1, a2, flat1.T, norm=colors.LogNorm(vmin=1.e-4,vmax=1.e-1))

        flat2 = xy_func(pcon)
        axes[0].set_xlabel('x')
        axes[0].set_ylabel('y')
        a1 = xar
        a2 = yar
        im2 = axes[0].pcolormesh(a1, a2, flat2.T, norm=colors.LogNorm(vmin=1.e-4,vmax=1.e-1))
        snaps.append([im1, im2])

    animated = anim.ArtistAnimation(fig, snaps, interval=50, blit=True, repeat_delay=300)

    if outname:
        animated.save(outname)
    else:
        return animated
