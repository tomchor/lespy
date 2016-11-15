from matplotlib import pyplot as plt
from matplotlib import animation as anim

def check_avgs(outputs, t_ini, t_end, savefigs=False, return_df=True,
        normalize=False, plot_kwargs={}, means=True, variances=True, covariances=False):
    """Plots important averages from the 2D outputs
    """
    from . import postProcess2D
    from . import simulation

    sim = simulation(outputs)
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


def xz_animation(arrays, outname=None):
    fig = plt.figure()
    snaps = []
    for tstamp, array in enumerate(arrays):
        xz_surf = array.mean(axis=1).T
        snaps.append(plt.imshow(xz_surf, xlabel="x nodes", ylabel="y nodes"))
    animated = anim.ArtistAnimation(fig, snaps, interval=50, blit=True, repeat_delay=300)

    if outname:
        animated.save(outname)
    else:
        return animated
