#from .decorators import autoassign

def paramParser(namelist):
    """Function that parses parameters from param.nml namelist files
    """
    from .f90nml import read
    if isinstance(namelist, str):
        namelist=open(namelist,'rt')

    nml = namelist.readlines()
    nml = [ line for line in nml if '%' not in line ]
    aux = open('param_buffer.nml', 'wt')
    aux.writelines(nml)
    aux.close()
    groups = read('param_buffer.nml')
    params = {}
    for key in groups.keys():
        params.update(groups[key])
    return params
