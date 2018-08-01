
def check_avgs(outputs, t_ini, t_end, savefigs=False, return_df=True, simulation=None, 
        normalize=False, plot_kwargs={}, means=True, variances=True, covariances=False, theta=False):
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
        whether or not to save the figs. If false, figs are just shown and plt.gcf() is returned instead.
    normalize: boolean
        whether to normalize height with inversion height.
    """
    from matplotlib import pyplot as plt
    from . import postProcess2D
    from . import Simulation
    from .utils import get_lims

    if simulation:
        sim = simulation
    else:
        sim = Simulation(outputs)
    t_ini = int(t_ini)
    t_end = int(t_end)
    dat = postProcess2D(outputs, t_ini=t_ini, t_end=t_end, simulation=sim, return_df=return_df)

    #------
    # The theta that comes out is actually more like rho/rho0, to we transform to absolute temperature
    dat.loc[:, '<Theta>'] = 2*sim.t_init - sim.t_scale*dat['<Theta>']
    #------

    #------
    # Inverts and maybe normalizes z axis
    if normalize:
        dat.iloc[:, :3] = -dat.iloc[:, :3]/sim.inversion_depth
    else:
        dat.iloc[:, :3] = -dat.iloc[:, :3]
    #------

    #------
    # plots the means
    if means:
        xlim = get_lims(dat['<U>/u*'])
        dat.plot(x='<U>/u*', y='z_uv', grid=True, xlim=xlim, **plot_kwargs)
        if savefigs: plt.savefig('U_{!s}-{!s}.png'.format(t_ini, t_end))
    
        xlim = get_lims(dat['<V>/u*'])
        dat.plot(x='<V>/u*', y='z_uv', grid=True, xlim=xlim, **plot_kwargs)
        if savefigs: plt.savefig('V_{!s}-{!s}.png'.format(t_ini, t_end))

        xlim = get_lims(dat['<W>/u*'])
        dat.plot(x='<W>/u*', y='z_uv', grid=True, xlim=xlim, **plot_kwargs)
        if savefigs: plt.savefig('W_{!s}-{!s}.png'.format(t_ini, t_end))
    #------
    
    #------
    # plots the variances
    if variances:
        xlim = get_lims(dat['<u^2>/u*^2'])
        dat.plot(x='<u^2>/u*^2', y='z_uv', xlim=xlim, grid=True, **plot_kwargs)
        if savefigs: plt.savefig('uu_{!s}-{!s}.png'.format(t_ini, t_end))
    
        xlim = get_lims(dat['<v^2>/u*^2'])
        dat.plot(x='<v^2>/u*^2', y='z_uv', xlim=xlim, grid=True, **plot_kwargs)
        if savefigs: plt.savefig('vv_{!s}-{!s}.png'.format(t_ini, t_end))
        
        xlim = get_lims(dat['<w^2>/u*^2'])
        dat.plot(x='<w^2>/u*^2', y='z_uv', xlim=xlim, grid=True, **plot_kwargs)
        if savefigs: plt.savefig('ww_{!s}-{!s}.png'.format(t_ini, t_end))
    #------
    
    #------
    # plots the covariances
    if covariances:
        xlim = get_lims(dat['<uw>/u*^2'])
        dat.plot(x='<uw>/u*^2', y='z_uv', xlim=xlim, grid=True, **plot_kwargs)
        if savefigs: plt.savefig('uw_{!s}-{!s}.png'.format(t_ini, t_end))

        xlim = get_lims(dat['<vw>/u*^2'])
        dat.plot(x='<vw>/u*^2', y='z_uv', xlim=xlim, grid=True, **plot_kwargs)
        if savefigs: plt.savefig('vw_{!s}-{!s}.png'.format(t_ini, t_end))

        xlim = get_lims(dat['<wT>/u*T*'])
        dat.plot(x='<wT>/u*T*', y='z_uv', xlim=xlim, grid=True, **plot_kwargs)
        if savefigs: plt.savefig('wT_{!s}-{!s}.png'.format(t_ini, t_end))
    #------

    #------
    # plots theta
    if theta:
        xlim = get_lims(dat['<Theta>'])
        dat.plot(x='<Theta>', y='z_uv', grid=True, xlim=xlim, **plot_kwargs)
        if savefigs: plt.savefig('Theta_{!s}-{!s}.png'.format(t_ini, t_end))
    #------
 
    if not savefigs:
        figs = plt.gcf()
        plt.show()
        plt.close()
        return figs
    else:
        return dat



def animate(results, outname=None, which='xy', simulation=None, title='', grid=True,
        clim=[], logscale=False, axes=[], levels=None, nbins=6, interpolation=None, interval=50,
        timelist=None, aspect='equal', cmap='viridis', clabel='', verbose=True, cbar_format=None, dpi=120):
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
    fig, ax = plt.subplots()
    sim = simulation

    results = np.asarray(results)
    nt = results.shape[0]

    #--------
    # This sets the labels for each frame
    if type(timelist)==type(None):
        timelist = range(0, nt)
    #--------

    #------------
    # Sets levels and ticks to use on the contour
    ticklabels, formatting, levels_con = get_ticks(results, levels=levels, logscale=logscale, clim=clim, nbins=nbins)
    xmin=sim.domain.x.min()
    xmax=sim.domain.x.max()
    #------------

    if which=='xz':
        results=results.transpose((0,2,1))[:,::-1]
        ymin=sim.domain.z.min()
        ymax=sim.domain.z.max()
    else:
        results = results.transpose(0,2,1)
        ymin=sim.domain.y.min()
        ymax=sim.domain.y.max()

    #-------
    # To take care of slightly negative or zero concentrations for plotting purposes
    global im, tx
    if logscale:
        from matplotlib.colors import LogNorm
        results = abs(results)
        results[ results==0 ] += 1.e-20
        im = plt.imshow(results[0], interpolation=interpolation, animated=True, cmap=cmap, origin='lower',
                norm=LogNorm(vmin=levels_con.min(), vmax=levels_con.max()), extent=[xmin, xmax, ymin, ymax])
    else:
        im = plt.imshow(results[0], interpolation=interpolation, animated=True, cmap=cmap, origin='lower',
                vmin=levels_con.min(), vmax=levels_con.max(), extent=[xmin, xmax, ymin, ymax])
    tx = im.axes.annotate(timelist[0], (.8,.03), annotation_clip=False, xycoords='figure fraction')
    im.axes.set_title(title)
    #-------

    def updatefig(it):
        global im, tx
        if verbose: print(timelist[it])
        im.set_array(results[it])
        tx.set_text(timelist[it])
        return im, tx

    #--------
    # Doesnt waste space on the plot and sets aspect ratio to be equal
    if which=='xz':
        plt.xlabel('$\mathbf{x(m)}$')
        plt.ylabel('$\mathbf{z(m)}$')
    elif which=='xy':
        plt.xlabel('$\mathbf{x(m)}$')
        plt.ylabel('$\mathbf{y(m)}$')
    cbar = plt.colorbar(label=clabel, format=cbar_format)
    plt.gca().set_aspect(aspect)
    plt.tight_layout()
    #--------

    #--------
    # Animate and put the colorbar
    animated = anim.FuncAnimation(fig, updatefig, frames=range(1,nt), interval=interval, blit=True)
    #--------

    if outname:
        animated.save(outname, bitrate=-1, dpi=dpi, writer='ffmpeg', codec='libx264', 
                extra_args=['-pix_fmt', 'yuv420p', '-g', '30', '-r', '30', '-profile:v', 'high'])
        #animated.save(outname) 
    else:
        print('returning')
        return animated



def pcon_2D_animation(results, outname=None, which='xy', simulation=None, title='',
        clim=[], logscale=True, axes=[], levels=None, nbins=6,
        timelist=None, aspect='equal', **kwargs):
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

    try:
        results = results.values
    except:
        pass

    #--------
    # This sets the labels for each frame
    if type(timelist)==type(None):
        timelist = range(0, results.shape[0])
    #--------

    #-------
    # To take care of slightly negative or zero concentrations for plotting purposes
    results[ results==0 ] += 1.e-20
    results = abs(results)
    #-------

    #------------
    # Sets levels and ticks to use on the contour
    ticklabels, formatting, levels_con = get_ticks(results, levels=levels, logscale=logscale, clim=clim, nbins=nbins)
    #------------

    for i, planecut in zip(timelist, results):
        print(i)
        aux = plane(planecut, which=which, axes=axes, outname=None, levels=levels_con, 
                simulation=simulation, logscale=logscale, set_cbar=False, **kwargs)
        timelabel = aux.ax.annotate(i, (.8,.05), annotation_clip=False, xycoords='figure fraction')
        aux.collections.append(timelabel)
        aux.grid(grid)
        snaps.append(aux.collections)
    fig.axes[0].set_title(title)

    #--------
    # Doesnt waste space on the plot and sets aspect ratio to be equal
    plt.gca().set_aspect(aspect)
    plt.tight_layout()
    #--------

    #--------
    # Animate and put the colorbar
    animated = anim.ArtistAnimation(fig, snaps, interval=50, blit=True)
    cbar = plt.colorbar(aux, ticks = levels_con, extend=None)
    cbar.ax.set_yticklabels(formatting(ticklabels))
    #--------

    if outname:
        animated.save(outname, bitrate=-1, dpi=150, writer='ffmpeg')
    else:
        return animated



def plane(plane, outname=None, which='xy', simulation=None, 
        axes=[], levels=None, logscale=True, clim=[], title='',
        set_cbar=True, cmap=None, xlim=[], ylim=[], aspect='equal', nbins=6, clabel=None):
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
        cmap_con = cm.get_cmap("viridis")

    #------------
    # Sets levels and ticks to use on the contour
    ticklabels, formatting, levels_con = get_ticks(plane, levels=levels, logscale=logscale, clim=clim, nbins=nbins)
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
    aux.ax.grid()
    #-------

    #-------
    # If this is going to animation, cbar shoudnt be set here
    if set_cbar:
        cbar = plt.colorbar(aux, ticks = levels_con, extend='both', label=clabel)
        cbar.ax.set_yticklabels(formatting(ticklabels))
        plt.gca().set_aspect(aspect)
        plt.tight_layout()
    #-------

    if outname:
        print('saving figure...', end='')
        plt.title(title)
        plt.savefig(outname, bbox_inches='tight')
        plt.close()
        print('done.')
        return
    else:
        return aux



def plane2(plane, outname=None, which='xy', simulation=None, interpolation=None,
        axes=[], levels=None, logscale=False, clim=[], title='',
        set_cbar=True, cmap='viridis', xlim=[], ylim=[], aspect='equal', 
        nbins=6, clabel='', axis=None, label_xy=True, norm_axes=False):
    """
    Prints 2D animations from binary data for oil concentrations
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from .utils import get_ticks

    sim = simulation
    if axis: plt.sca(axis)

    #------------
    # Sets levels and ticks to use on the contour
    ticklabels, formatting, levels_con = get_ticks(plane, levels=levels, logscale=logscale, clim=clim, nbins=nbins)
    vmin=levels_con.min()
    vmax=levels_con.max()
    #------------

    #------------
    # Sets x y grid
    if axes:
        ax1, ax2 = axes
    else:
        if simulation:
            if which == 'xy':
                ax1 = sim.domain.x
                ax2 = sim.domain.y
            elif which == 'xz':
                ax1 = sim.domain.x
                ax2 = sim.domain.z
            else:
                ax1 = np.arange(0, plane.shape[0])
                ax2 = np.arange(0, plane.shape[1])
            if norm_axes:
                ax1/=sim.inv_depth
                ax2/=sim.inv_depth
        else:
            ax1 = np.arange(0, plane.shape[0])
            ax2 = np.arange(0, plane.shape[1])
    xmin=ax1.min()
    ymin=ax2.min()
    xmax=ax1.max()
    ymax=ax2.max()
    #------------

    #----------
    # Sets labels latex style (no label is set if which is not given
    if label_xy:
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
        from matplotlib.colors import LogNorm
        im = plt.imshow(plane.T, interpolation=interpolation, animated=True, cmap=cmap, origin='lower',
                                norm=LogNorm(vmin=vmin, vmax=vmax), extent=[xmin, xmax, ymin, ymax])
    else:
        im = plt.imshow(plane.T, interpolation=interpolation, animated=True, cmap=cmap, origin='lower',
                                vmin=vmin, vmax=vmax, extent=[xmin, xmax, ymin, ymax])
    #----------

    #-------
    # Actual plotting is done here
    im.axes.set_xlim(*xlim)
    im.axes.set_ylim(*ylim)
    im.axes.grid()
    #-------

    #-------
    # If this is going to animation, cbar shoudnt be set here
    if set_cbar:
        cbar = plt.colorbar(label=clabel)
    plt.gca().set_aspect(aspect)
    plt.title(title)
    plt.tight_layout()
    #-------

    if outname:
        print('saving figure...', end='')
        plt.savefig(outname, bbox_inches='tight')
        plt.close()
        print('done.')
        return im
    else:
        return im



def plot_edls(lsc, ax=None, extent=None, logscale=True, clim=[None, None], aspect=None, interp=None):
    """Plot list of endless patches in their correct place.
    The patches must be xarray.DataArrays
    """
    import numpy as np
    import matplotlib.pyplot as plt
    options = dict(vmin=clim[0], vmax=clim[1], add_colorbar=True, x='x', y='y', interpolation=None, aspect=aspect)
    options2 = dict(vmin=-0.02, vmax=0.02, add_colorbar=True, x='x', y='y', interpolation=None, aspect=aspect, cmap="seismic")

    #------
    # Determine extent
    if type(extent)==type(None):
        minx, maxx, miny, maxy = [], [], [], []
        for arr3 in lsc:
            minx.append(arr3.coords['x'].values.min())
            maxx.append(arr3.coords['x'].values.max())
            miny.append(arr3.coords['y'].values.min())
            maxy.append(arr3.coords['y'].values.max())
        extent =  min(minx), max(maxx), min(miny), max(maxy)
    #------

    #------
    # Determine scale
    if logscale:
        from matplotlib.colors import LogNorm
        options['norm'] = LogNorm()
    print(options)
    #------

    if type(ax)==type(None):
        fig, axes = plt.subplots(ncols=1, figsize=(8,8), sharex=True, sharey=True, squeeze=False)
        axes=axes.flatten()
    else:
        axes=ax
    print(axes)

    for i, patch in enumerate(lsc):
        patch=abs(patch)
        patch.plot.imshow(ax=axes[0], **options)
        if i==0:
            options['add_colorbar']=False
            options2['add_colorbar']=False

    for ax in axes:
        ax.set_xlim(*extent[:2])
        ax.set_ylim(*extent[2:])
        ax.set_xticks(np.arange(extent[0],extent[1], 500), minor=True)
        ax.set_yticks(np.arange(extent[2],extent[3], 500), minor=True)
        ax.set_xticks(np.arange(extent[0],extent[1], 1000), minor=False)
        ax.set_yticks(np.arange(extent[2],extent[3], 1000), minor=False)
        ax.set_aspect(1)
    axes[0].grid(True, which="both", color="w")
    return ax

