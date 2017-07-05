__version__ = '1.0.1'
print('Initiating lespy version %s'%__version__, end='...')
from .core import *
from .routines import *
from . import core, utils, routines, plot, physics
from . import vector, numerical, stats, spectral#, nml
#from . import convection
#from .fortran import numfor
#import xarray
#from .pyevtk.hl import gridToVTK
#__all__ = ["dmClass", "concentration", "langmuir"]
print('done')
