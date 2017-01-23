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

    if 'vel_sc' in fname:
        ndtime = fname.strip('vel_sc').strip('.out')
        return int(ndtime)

    if 'vel_t' in fname:
        ndtime = fname.strip('vel_t').strip('.out')
        return int(ndtime)

    if 'temp_t' in fname:
        ndtime = fname.strip('temp_t').strip('.out')
        return int(ndtime)

    if 'con_tt' in fname:
        import re
        numbers = re.split('[a-z. +_]',fname)
        ndtime, pcon_n, row, col = [ int(el) for el in numbers if el is not '' ]
        return int(ndtime), pcon_n, row, col


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
                #clim = (np.nanmin(logarray), np.nanmax(logarray))
                clim = (np.nanpercentile(logarray, 50), np.nanpercentile(logarray, 95))
            levels_con = np.linspace(clim[0], clim[1], nseps)
            ticklabels = np.power(10.0, levels_con)
        else:
            if not clim:
                clim=(array.min(), array.max())
            levels_con = np.linspace(clim[0],clim[1], nseps)
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



def get_lims(data, increase=.15):
    """
    Quick and simple function that gets nice limits to plot the data in
    """
    totdelta = abs(data.max() - data.min())
    incr = increase*totdelta
    botlim = data.min() - incr
    toplim = data.max() + incr

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
            for key, val in arrays.items():
                points = { key:val.sel(time=tstep).values }
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


