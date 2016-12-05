#from .decorators import autoassign

def paramParser(nmlpath):
    """Function that parses parameters from param.nml namelist files
    """
    from .f90nml import read
    from os.path import basename, isfile, exists, join
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
                print('Found param.nml at {}'.format(join(nmlpath, 'param.nml')))
            except:
                namelist = open(join(nmlpath, 'codebkp/param.nml'), 'rt')
                print('Found param.nml at {}'.format(join(nmlpath, 'codebkp/param.nml')))
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
        #ndtime = fname.split('vel_sc')[1].split('.out')[0]
        ndtime = fname.strip('vel_sc').strip('.out')
        return int(ndtime)

    if 'con_tt' in fname:
        import re
        numbers = re.split('[a-z. +_]',fname)
        ndtime, pcon_n, row, col = [ int(el) for el in numbers if el is not '' ]
        return ndtime, pcon_n, row, col


