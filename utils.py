#from .decorators import autoassign

def paramParser(nmlpath):
    """Function that parses parameters from param.nml namelist files
    """
    from .f90nml import read
    from os.path import basename

    buffername = nmlpath.replace('/', '_') + '.buffer'
    if isinstance(nmlpath, str):
        namelist=open(nmlpath,'rt')

    nml = namelist.readlines()
    nml = [ line for line in nml if '%' not in line ]
    aux = open(buffername, 'wt')
    aux.writelines(nml)
    aux.close()
    groups = read(buffername)
    params = {}
    for key in groups.keys():
        params.update(groups[key])
    return params
