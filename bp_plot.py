from matplotlib import pyplot as plt
from matplotlib import animation as anim

def check_avgs(outputs, t_ini, t_end, savefigs=False, return_df=True,
        normalize=False, plot_kwargs={}, means=True, variances=True, covariances=False):
    """Plots important averages from the 2D outputs
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



def pcon_2D_animation(results, outname=None, which='xy', simulation=None, ax1=None, ax2=None, levels=None,
        logscale=True, timelist=None):
    """
    Prints 2D animations from binary data for oil concentrations
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as anim
    from matplotlib import colors
    from matplotlib import cm
    from . import routines
    from . import utils

    fig = plt.figure()
    snaps = []
    sim = simulation
    cmap_con = cm.get_cmap("jet")

    if not timelist:
        timelist = results.shape[0]

    if levels==None:
        if logscale:
            level_bd=(-3,1)
            levels_con = np.arange(level_bd[0],level_bd[1],0.4)
            ticklabels = np.power(10.0, levels_con)
            formatting_fun = np.vectorize(lambda f: format(f, '4.0e'))

    if simulation:
        if which == 'xy':
            ax1p = np.arange(0, results.shape[1]*sim.domain.dx, sim.domain.dx)
            ax2p = np.arange(0, results.shape[2]*sim.domain.dy, sim.domain.dy)
        elif which == 'xz':
            ax1p = np.arange(0, results.shape[1]*sim.domain.dx, sim.domain.dx)
            ax2p = np.arange(0, results.shape[2]*sim.domain.dz, sim.domain.dz)
    else:
        ax1p = np.arange(0, results.shape[1])
        ax2p = np.arange(0, results.shape[2])

    if ax1 is None:
        ax1 = ax1p
    if ax2 is None:
        ax2 = ax2p
    X2, X1 = np.meshgrid(ax2, ax1)

    #-------
    # To take care of slightly negative or zero concentrations for plotting purposes
    results[ results==0 ] += 1.e-20
    results = abs(results)
    for i, plane in zip(timelist, results):
        print(i)

        #-------

        if which=='xz':
            plt.xlabel('$\mathbf{x(m)}$')
            plt.ylabel('$\mathbf{z(m)}$')
        elif which=='xy':
            plt.xlabel('$\mathbf{x(m)}$')
            plt.ylabel('$\mathbf{y(m)}$')

        if logscale:
            plane_log = np.log10(plane)
        else:
            plane_log = plane
        plt.tight_layout()
        aux = plt.contourf(X1, X2, plane_log, levels_con, cmap=cmap_con, extend='both')
        title = aux.ax.annotate(i, (1000,-30), annotation_clip=False)

        aux.collections.append(title)
        snaps.append(aux.collections)

    animated = anim.ArtistAnimation(fig, snaps, interval=50, blit=True)
    cbar = plt.colorbar(aux, ticks = levels_con, extend='both')
    cbar.ax.set_yticklabels(formatting_fun(ticklabels))

    if outname:
        animated.save(outname, bitrate=-1, dpi=120)
    else:
        return animated



def plane(plane, outname=None, which='xy', simulation=None, ax1=None, ax2=None, levels=None, logscale=True):
    """
    Prints 2D animations from binary data for oil concentrations
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as anim
    from matplotlib import colors
    from matplotlib import cm
    from . import routines
    from . import utils

    fig = plt.figure()
    snaps = []
    sim = simulation
    cmap_con = cm.get_cmap("jet")


    if levels==None:
        if logscale:
            level_bd=(-3,1)
            levels_con = np.arange(level_bd[0],level_bd[1],0.4)
            ticklabels = np.power(10.0, levels_con)
            formatting_fun = np.vectorize(lambda f: format(f, '4.0e'))
        else:
            level_bd=(1e-6,1.)
            levels_con = np.arange(level_bd[0],level_bd[1],0.3)
#        levels_con = np.logspace(level_bd[0], level_bd[1], 30)

    if simulation:
        if which == 'xy':
            ax1p = np.arange(0, plane.shape[0]*sim.domain.dx, sim.domain.dx)
            ax2p = np.arange(0, plane.shape[1]*sim.domain.dy, sim.domain.dy)
        elif which == 'xz':
            ax1p = np.arange(0, plane.shape[0]*sim.domain.dx, sim.domain.dx)
            ax2p = np.arange(0, plane.shape[1]*sim.domain.dz, sim.domain.dz)
    else:
        ax1p = np.arange(0, plane.shape[0])
        ax2p = np.arange(0, plane.shape[1])

    if ax1 is None:
        ax1 = ax1p
    if ax2 is None:
        ax2 = ax2p
    X2, X1 = np.meshgrid(ax2, ax1)

    plane = np.ma.array(plane, mask=plane==np.nan)
    #-------
    # To take care of slightly negative or zero concentrations for plotting purposes
    plane[ plane==0 ] += 1.e-20
    plane = abs(plane)
    #-------

    if which=='xz':
        plt.xlabel('$\mathbf{x(m)}$')
        plt.ylabel('$\mathbf{z(m)}$')
    elif which=='xy':
        plt.xlabel('$\mathbf{x(m)}$')
        plt.ylabel('$\mathbf{y(m)}$')

    if logscale:
        plane_log = np.log10(plane)
    else:
        plane_log = plane

    aux = plt.contourf(X1, X2, plane_log, levels_con, cmap=cmap_con, extend='both')
    cbar = plt.colorbar(aux, ticks = levels_con, extend='both')
    cbar.ax.set_yticklabels(formatting_fun(ticklabels))

    if outname:
        return plt.savefig(outname, bbox_inches='tight')
    else:
        return aux




def pcon_2D_animation_old(bins, outname=None, which='xz', simulation=None,
        n_pcon=0, func=lambda x: x.mean(axis=1), trim_x=True):
    """
    Prints 2D animations from binary data for oil concentrations
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as anim
    from matplotlib import colors
    from . import routines
    from . import utils

    fig = plt.figure()
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
        ndtime = utils.nameParser(fname)

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

        flat = func(pcon)
        if which=='xz':
            plt.xlabel('x')
            plt.ylabel('z')
            a2 = -zar
            a1 = xar
        elif which=='xy':
            plt.xlabel('x')
            plt.ylabel('y')
            a1 = xar
            a2 = yar
        else:
            a1=range(flat.shape[0])
            a2=range(flat.shape[1])
        aux = plt.pcolormesh(a1, a2, flat.T, norm=colors.LogNorm(vmin=1.e-4,vmax=1.e-1))

        snaps.append([aux])


    animated = anim.ArtistAnimation(fig, snaps, interval=50, blit=True, repeat_delay=300)

    if outname:
        animated.save(outname)
    else:
        return animated


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
