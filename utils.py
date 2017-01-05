#from .decorators import autoassign

def paramParser(nmlpath):
    """Function that parses parameters from param.nml namelist files
    """
    from .f90nml import read
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
    if levels==None:
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


def np2vtr(array, outname):
    """Writes a numpy array in vtr format"""
    from .pyevtk.hl import gridToVTK
    return
