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


def pcon_2D_animation(bins, outname=None, which='xz', simulation=None,
        n_pcon=0, func=lambda x: x.mean(axis=1), trim_x=True):
    """
    Prints 2D animations from binary data for oil concentrations
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as anim
    from matplotlib import colors
    from . import routines

    fig = plt.figure()
    snaps = []

    sim = simulation
    if trim_x:
        xar = np.arange(0,sim.domain.Ld)*sim.domain.dx
    else:
        xar = np.arange(0,sim.domain.Nx)*sim.domain.dx
    yar = (np.arange(0,sim.ny)*sim.domain.dy).reshape(-1,1)
    zar = (np.arange(0,sim.nz_tot)*sim.domain.dz).reshape(-1,1)


    for fname in bins:
        print(fname)
        u,v,w,T,pcon = routines.readBinary(fname, simulation=sim)
        pcon = pcon[:,:,:, n_pcon]

        #-------
        # To take care of slightly negative or zero concentrations for plotting purposes
        pcon[ pcon==0 ] += 1.e-20
        pcon = abs(pcon)
        #-------


        plt.xlabel('x')
        if which=='xz' or which=='zx':
            flat = pcon.mean(axis=1)[:-2]
            plt.ylabel('z')
            snaps.append(( plt.pcolormesh(xar, -zar, flat.T, norm=colors.LogNorm(vmin=1.e-4,vmax=1.e-1)), ))
        elif which=='xy':
            flat = pcon.mean(axis=2)[:-2]
            plt.ylabel('y')
            snaps.append(( plt.pcolormesh(flat.T, norm=colors.LogNorm(vmin=1.e-4,vmax=1.e-1)), ))
        else:
            raise Exception

    animated = anim.ArtistAnimation(fig, snaps, interval=50, blit=True, repeat_delay=300)

    if outname:
        animated.save(outname)
    else:
        return animated
