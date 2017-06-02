def paramParser(nmlpath):
    """Function that parses parameters from param.nml namelist files
    """
    from .nml import read
    from os.path import basename, isfile, exists, join, abspath
    from os import remove

    buffername = nmlpath.replace('/', '_') + '.buffer'
    if exists(nmlpath):
        #----------
        # If nmlpath is a file, it's probably just the param.nml file
        if isfile(nmlpath):
            namelist=open(nmlpath, 'rt')
        #----------

        #----------
        # If it's a dir, we first look the param.nml in the dir and then in the codebkp dir
        else:
            try:
                namelist = open(join(nmlpath, 'param.nml'), 'rt')
                print('Found param.nml at {}'.format(abspath('param.nml')))
            except:
                namelist = open(join(nmlpath, 'codebkp/param.nml'), 'rt')
                print('Found param.nml at {}'.format(abspath('codebkp/param.nml')))
        #----------

    else:
        raise ValueError('Path {} doesnt exist'.format(nmlpath))

    #----------
    # lines with % are not supported (type character resolver)
    # rls_src has variables in it, so python can't understand it
    nml = namelist.readlines()
    nml = [ line for line in nml if '%' not in line ]
    nml = [ line for line in nml if 'rls_src' not in line ]
    aux = open(buffername, 'wt')
    aux.writelines(nml)
    aux.close()
    #----------

    groups = read(buffername)
    remove(buffername)
    params = {}
    for key in groups.keys():
        params.update(groups[key])
    return params


def nameParser(fname):
    """
    Parser that gets name of output file and returns time, endless patch location and etc
    """
    from os import path
    fname = path.basename(fname)

    if 'con_tt' in fname:
        import re
        numbers = re.split('[a-z. +_]',fname)
        ndtime, pcon_n, row, col = [ int(el) for el in numbers if el is not '' ]
        return int(ndtime), pcon_n, row, col

    if 'vel_sc' in fname:
        start='vel_sc'

    if 'vel_t' in fname:
        start='vel_t'

    if 'temp_t' in fname:
        start='temp_t'

    if 'div_z0_t' in fname:
        start='div_z0_t'


    ndtime = fname.strip(start).strip('.out')
    return int(ndtime)


def get_ticks(array, levels=None, logscale=None, clim=[], nbins=6):
    """
    Auxiliar function to get plitting limits for plane and pcon_2D_animation
    """
    import numpy as np

    #-------
    nseps = nbins+1
    if type(levels)==type(None):
        if logscale:
            if not clim:
                logarray = np.log10(np.abs(array))
                clim = (np.nanpercentile(logarray, 20), np.nanpercentile(logarray, 95))
            levels_con = np.linspace(clim[0], clim[1], nseps)
            ticklabels = np.power(10.0, levels_con)
        else:
            if not clim:
                #clim=(float(array.min()), float(array.max()))
                #clim=get_lims(array, method='median', increase=0.)
                clim = (np.nanpercentile(array, 20), np.nanpercentile(array, 95))
            levels_con = np.linspace(clim[0], clim[1], nseps)
            ticklabels = levels_con
    #-------

    #-------
    # If levels are given, result is straightfoward
    else:
        levels_con = levels
        if logscale:
            ticklabels = np.power(10.0, levels_con)
        else:
            ticklabels = levels_con
    #-------
    
    #-------
    # The formatting can be done differently
    if logscale:
        formatting = np.vectorize(lambda f: format(f, '4.0e'))
    else:
        formatting = np.vectorize(lambda f: format(f, 'G'))
    #-------

    return ticklabels, formatting, levels_con



def get_lims(data, increase=.15, method='bounds'):
    """
    Quick and simple function that gets nice limits to plot the data in
    """
    import numpy as np
    if method=='bounds':
        totdelta = abs(float(data.max()) - float(data.min()))
        incr = increase*totdelta
        botlim = data.min() - incr
        toplim = data.max() + incr
    elif method=='median':
        med = float(np.median(data))
        maxval = float(np.max(data))
        minval = float(np.min(data))
        delta = (1.-increase)*(abs(abs(med)-abs(maxval)) + abs(abs(med)-abs(minval)))/2.
        botlim = med - delta
        toplim = med + delta

    return botlim, toplim



def find_in_tree(name, path, type='f'):
    """Find name inside path and subdirectories"""
    import os

    if os.path.basename(path)==name:
        result = [os.path.abspath(path)]
    else:
        result = []

    if type=='f':
        for root, dirs, files in os.walk(path):
            if name in files:
                result.append(os.path.abspath(os.path.join(root, name)))
    if type=='d':
        for root, dirs, files in os.walk(path):
            if name in dirs:
                result.append(os.path.abspath(os.path.join(root, name)))

    return result



def np2vtr(arrays, outname):
    """
    Writes an xarray DataArray in vtr format
    The input MUST be a dict with xarray DataArrays:
    dataarray = xr.DataArray(np_array, dims=['time', 'x', 'y'], coords={'time':timestamps', 'x':x_array, 'y':y_array})
    """
    from .pyevtk.hl import gridToVTK
    coords = list(arrays.values())[0].coords
    x = coords['x'].values
    y = coords['y'].values
    z = coords['z'].values

    if 'time' not in coords.dims:
        try:
            gridToVTK(outname, x,y,z, pointData = arrays)
        except AssertionError:
            arrays = { key:val.values for key,val in arrays.items() }
            gridToVTK(outname, x,y,z, pointData = arrays)
    else:
        timestamps = coords['time']
        for tstep in timestamps:
            tstep = int(tstep)
            print('Writing t=',tstep,'to vtr')
            points={}
            for key, val in arrays.items():
                points.update({ key:val.sel(time=tstep).values })
            gridToVTK(outname.format(tstep), x,y,z, pointData = points)

    return



def get_dataarray(pcons, simulation=None, with_time=False):
    """
    Gets a dataarray from pcons
    
    pcons: list or np.array
        If it's a list the domains can be different
    with_time: list, array
        list that will serve as the time index
    """
    from collections import OrderedDict as ODic
    sim = simulation

    if with_time:
        time=1
    else:
        time=0

    #--------
    # This is just to find out the dimensions
    if isinstance(pcons, list):
        pcon = pcons[0]
    else:
        pcon = pcons

    if with_time:
        coords = ODic({'time':with_time})
    else:
        coords = ODic({})
    #--------

    #--------
    # If there is more than one dataset (list), each can have a different domain
    if isinstance(pcons, list):
        pcons_da = []
        for pcon in pcons:
            x,y,z = sim.domain.makeAxes(pcon[0])
            if (len(pcon.shape)-time) >= 1:
                coords.update({'x':x})
            if (len(pcon.shape)-time) >= 2:
                coords.update({'y':y})
            if (len(pcon.shape)-time) == 3:
                coords.update({'z':z})
            pcons_da.append(sim.DataArray(pcon, coords=coords))
        del pcons, pcon
        return pcons_da
    #--------

    #----------
    # Else, we do it just for one
    else:
        print(pcons.shape)
        if (len(pcon.shape)-time)>=1:
            coords.update({'x':sim.domain.x})
        if (len(pcon.shape)-time) >= 2:
            coords.update({'y':sim.domain.y})
        #-------
        # If pcons in not list, then the last dim is always the size
        if (len(pcon.shape)-time) >= 3:
            import numpy as np
            if (len(pcon.shape) - time) == 4:
                coords.update({'z':sim.domain.z})
            coords.update({'size':np.arange(pcon.shape[-1])})
        #-------
        if (len(pcon.shape)-time) == 5:
            raise ValueError("Too many dimensions in array")
        return sim.DataArray(pcons, coords=coords)
    #----------


def radial_prof3D(data, r=None, simulation=None, func=None):
    """
    Gets a radial profile around the center of `data`.

    Parameters
    ----------
    data: np.array
        3d array, shape = (nt, nx, ny)
    r: np.array
        matrix of distances from the point. Same shape as `data`.
    func: function
        defaults to np.mean
    """
    import numpy as np
    if type(func)==type(None):
        func=np.mean
    if type(r)==type(None):
        sim=simulation
        nt, Lx1, Ly1 = (np.array(data.shape)).astype(int)
        Lx, Ly = (np.array([Lx1, Ly1])/2).astype(int)
        x = np.arange(-Lx,-Lx+Lx1,1)*sim.domain.dx
        y = np.arange(-Ly,-Ly+Ly1,1)*sim.domain.dy
        xx, yy = np.meshgrid(x, y)
        r = np.sqrt(xx**2. + yy**2.)
    uniq = np.unique(r)
    prof = np.array([ np.mean(data[:, r==un ], axis=1) for un in uniq ]).T
    return uniq, prof



def radial_prof(data, r=None, simulation=None, func=None, axes=(0,1)):
    """
    Gets a radial profile around the center of `data`.

    Parameters
    ----------
    data: np.array
        2d array
    r: np.array
        matrix of distances from the point. Same shape as `data`.
    func: function
        defaults to np.mean
    """
    import numpy as np
    axes = list(axes)
    if type(func)==type(None):
        func=np.mean
    if type(r)==type(None):
        sim=simulation
        Lx1, Ly1 = (np.array(data.shape)[axes]).astype(int)
        Lx, Ly = (np.array(data.shape)[axes]/2).astype(int)
        x = np.arange(-Lx,-Lx+Lx1,1)*sim.domain.dx
        y = np.arange(-Ly,-Ly+Ly1,1)*sim.domain.dy
        xx, yy = np.meshgrid(x, y)
        r = np.sqrt(xx**2. + yy**2.)
    uniq = np.unique(r)
    prof = np.array([ func(data[ r==un ]) for un in uniq ])
    return uniq, prof


import numpy as _np
def classbin(x, y, bins_number=100, function=_np.mean, xfunction=_np.mean, logscale=True, clean_nans=False):
    """
    Separates x and y inputs into bins based on the x array.
    x and y do not have to be ordered.

    Parameters
    -----------
    x: np.array
        independent variable
    y: np.array
        dependent variable
    bins_number: int
        number of classes (or bins) desired
    function: callable
        funtion to be applied to both x and y-bins in order to smooth the data
    logscale: boolean
        whether or not to use a log-spaced scale to set the bins

    Returns
    -------
    np.array:
        x binned
    np.array:
        y binned
    """
    import warnings
    import numpy as np

    xmin=np.min(x)
    xmax=np.max(x)
    if logscale:
        #-----------
        # The following if statement gets rid of negative or zero values in the x array, since we are using log-scale
        if (x<=0).any():
            y=np.array([ yy for yy, xx in zip(y,x) if xx > 0 ])
            x=np.array([ el for el in x if el > 0])
            xmin=np.min(x)
            xmax=np.max(x)
        #-----------
        bins=np.logspace(np.log(xmin), np.log(xmax), bins_number+1, base=np.e)
    else:
        bins=np.linspace(xmin, xmax, bins_number+1)
    xsm = np.zeros(bins_number)
    ysm = np.zeros(bins_number)



    #-----------
    # The following process is what actually bins the data using numpy
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for i in range(bins_number):
            if i == bins_number - 1:
                sel = bins[i] <= x
            else:
                sel = (bins[i] <= x) & (x < bins[i+1])
            xsm[i] = xfunction(x[sel])
            ysm[i] = function(y[sel])
    #-----------

    if clean_nans:
        xsm = xsm[ np.isfinite(ysm) ]
        ysm = ysm[ np.isfinite(ysm) ]

    return xsm, ysm


def nearest(array, values, return_idx=False):
    """ Searches for the nearest instances of values inside array """
    import numpy as np
    values=np.asarray(values)
    array2, values2 = np.meshgrid(array, values)
    idx = (np.abs(array2-values2)).argmin(axis=1)
    if return_idx:
        return idx
    else:
        return array[idx]



from matplotlib.colors import LinearSegmentedColormap as _lin_cmap

cdict1 = {'red':   ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 0.1),
                   (1.0, 1.0, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 1.0),
                   (0.5, 0.1, 0.0),
                   (1.0, 0.0, 0.0))
        }

cdict2 = {'red':   ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 1.0),
                   (1.0, 0.1, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 0.1),
                   (0.5, 1.0, 0.0),
                   (1.0, 0.0, 0.0))
        }

blue_red1 = _lin_cmap('BlueRed1', cdict1)
blue_red2 = _lin_cmap('BlueRed2', cdict2)

