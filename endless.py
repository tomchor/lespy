import numpy as _np
def get_edls_minmax(arr):
    arr=np.array(arr)
    minx, maxx, miny, maxy = [], [], [], []
    for arr1 in arr:
        for arr2 in arr1:
            for arr3 in arr2:
                minx.append(arr3.coords['x'].values.min())
                maxx.append(arr3.coords['x'].values.max())
                miny.append(arr3.coords['y'].values.min())
                maxy.append(arr3.coords['y'].values.max())
    return min(minx), max(maxx), min(miny), max(maxy)


def plot_edls(lsc, ax=None, extent=None, logscale=True, clim=[None, None], aspect=None, grid=True, interp=None):
    """Plot list of endless patches in their correct place.
    The patches must be xarray.DataArrays
    """
    options = dict(vmin=clim[0], vmax=clim[1], add_colorbar=True, x='x', y='y', interpolation=None, aspect=aspect)

    if type(extent)==type(None):
        minmaxs = get_edls_minmax(lsc)
    #------
    # Determine extent
#    if type(extent)==type(None):
#        for arr3 in lsc:
#            minx.append(arr3.coords['x'].values.min())
#            maxx.append(arr3.coords['x'].values.max())
#            miny.append(arr3.coords['y'].values.min())
#            maxy.append(arr3.coords['y'].values.max())
#        extent =  min(minx), max(maxx), min(miny), max(maxy)
    #------

    #------
    # Determine scale
    if logscale:
        from matplotlib.colors import LogNorm
        options['norm'] = LogNorm()
    print(options)
    #------

    if type(ax)==type(None):
        fig, ax = plt.subplots()

    for i, patch in enumerate(lsc):
        patch=abs(patch)
        patch.plot.imshow(ax=ax, **options)
        if i==0: options['add_colorbar']=False
    ax.set_xlim(*extent[:2])
    ax.set_ylim(*extent[2:])
    ax.grid(grid)
    return ax



