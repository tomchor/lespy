from .dmClass import domain
from .simClass import simulation
from .routines import *
from . import dmClass, simClass, concentration, langmuir, utils, routines
__all__ = ["dmClass", "concentration", "langmuir"]
__version__ = '0.0'
print('Initiating lespy version %s'%__version__)
