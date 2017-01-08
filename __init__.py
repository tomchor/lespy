from .core import *
from .routines import *
#from .vector import *
#from .utils import *
#from .plot import *
from . import core, utils, routines, plot, physics, vector, nml
from .pyevtk.hl import gridToVTK
import xarray
__all__ = ["dmClass", "concentration", "langmuir"]
__version__ = '0.0'
print('Initiating lespy version %s'%__version__)
