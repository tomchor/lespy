__version__ = '0.0'
print('Initiating lespy version %s'%__version__, end='...')
from .core import *
from .routines import *
from . import core, utils, routines, plot, physics
from . import vector, numerical, nml
from . import convection
#import xarray
#from .pyevtk.hl import gridToVTK
#__all__ = ["dmClass", "concentration", "langmuir"]
print('done')
