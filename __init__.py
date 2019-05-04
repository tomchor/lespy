__version__ = '2.0.0'
print('Initiating lespy version %s'%__version__, end='...')
from .core import *
from .routines import *
from . import core, utils, routines, plot, plot2, physics
from . import vector, numerical, stats, spectral#, nml
from .io import *
print('done')
