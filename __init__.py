from .core import *
from .routines import *
#from .vector import *
#from .utils import *
#from .plot import *
from . import core, utils, routines, plot, physics, vector
#from .pyevtk.hl import GridToVTK
import xarray
__all__ = ["dmClass", "concentration", "langmuir"]
__version__ = '0.0'
print('Initiating lespy version %s'%__version__)
