import numpy as _np
from .. import physics as phs

class driftVel:
    """class for Stokes drift velocity"""
    def __init__(self, simulation=None, wl=None, wn=None, amp=None, sigma=None, us0=None, g=phs.g):
        if simulation:
            self.wl = simulation.lambda_w
            self.amp = simulation.amp_w
        else:
            self.wl = wl
            self.wn = wn
            self.amp = amp
            self.sigma = sigma
        self.g = g

        if (self.wl == None):
            self.wn2wl()
        if (self.wn == None):
            self.wl2wn()
        if (self.sigma == None):
            self.sigma = self.get_sigma()
        if (self.amp !=None) and (self.wn != None):
            self.Us_function = phs.get_Ustokes(self.amp, self.wn, g=self.g)
        self.Us0 = self.get_Us0()
        self.efolding_z = self.get_efolding()


    def wl2wn(self):
        if (self.wl != None):
            self.wn = 2 * _np.pi / self.wl 

    def wn2wl(self):
        if (self.wn != None):
            self.wl = 2 * _np.pi / self.wn

    def get_sigma(self):
        """Gets the angular frequency of the waves sqrt(g*wn)"""
        if (self.wn != None):
            return _np.sqrt(self.g * self.wn)

    def get_Us0(self):
        if (self.wn != None and self.amp != None and self.sigma != None):
            return self.sigma * self.wn * self.amp ** 2

    def get_efolding(self):
        """Gets e-folding depth for convenience"""
        if (self.wn != None):
            return 1./(2.*self.wn)

    def Us(self, z, angle=None):
        """Returns the drift velocity at depth z"""
        if (self.wn != None and self.us0 != None):
            z = _np.asarray(z)
            us = self.us0 * _np.exp(2 * self.wn * z)
            return (us*_np.cos(angle), us*_np.sin(angle))

    def dUsdz(self, z, angle=None):
        """Retuns the derivative of the stokes drift at height z (meters)"""
        if (self.wn != None and self.us0 != None):
            z = _np.asarray(z)
            dusdz = 2 * self.wn * self.us0 * _np.exp(2 * self.wn * z)
            if angle:
                angle = _np.radians(angle)
                return (dusdz*_np.cos(angle), dusdz*_np.sin(angle))
            else:
                return dusdz

    def show(self):
        print('#' * 10,' Stokes Drift Parameters ','#' * 10)
        print("The gravitational acceleration is ",self.g," m/s^2")
        print("The wavelength is ",self.wl, " m")
        print("The wavenumber is ",self.wn, " /m")
        print("The frequency is ",self.sigma, " /s")
        print("The velocity of Stokes drift at the surface is ",self.us0,
                " m/s")
        print('#' * 10,' END Stokes Drift Parameters ','#' * 10)

