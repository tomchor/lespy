#!/usr/bin/env python3

### Import module ###
import numpy as np
### END Import module ###

### Class definition ###
class stokesDrift:
  """class for Stokes drift"""
  def __init__(self, wl=-1, wn=-1, amp=-1, sigma=-1, us0=-1, g=9.81):
    self.wl = wl
    self.wn = wn
    self.amp = amp
    self.sigma = sigma
    self.us0 = us0
    self.g = g
    if (self.wl == -1):
      self.wn2wl()
    if (self.wn == -1):
      self.wl2wn()
    if (self.sigma == -1):
      self.getSigma()
    if (self.us0 == -1):
      self.getUs0()

  def show(self):
    print('#' * 10,' Stokes Drift Parameters ','#' * 10)
    print("The gravitational acceleration is ",self.g," m/s^2")
    print("The wavelength is ",self.wl, " m")
    print("The wavenumber is ",self.wn, " /m")
    print("The frequency is ",self.sigma, " /s")
    print("The velocity of Stokes drift at the surface is ",self.us0,
        " m/s")
    print('#' * 10,' END Stokes Drift Parameters ','#' * 10)

  def wl2wn(self):
    if (self.wl != -1):
      self.wn = 2 * np.pi / self.wl 

  def wn2wl(self):
    if (self.wn != -1):
      self.wl = 2 * np.pi / self.wn

  def getSigma(self):
    if (self.wn != -1):
      self.sigma = np.sqrt(self.g * self.wn)

  def getUs0(self):
    if (self.wn != -1 and self.amp != -1 and self.sigma != -1):
      self.us0 = self.sigma * self.wn * self.amp ** 2

  def getUs(self, z, ag_rad):
    if (self.wn != -1 and self.us0 != -1):
      z = np.asarray(z)
      us = self.us0 * np.exp(2 * self.wn * z)
      return (us*np.cos(ag_rad), us*np.sin(ag_rad))

  def getDusdz(self, z, ag_rad):
    if (self.wn != -1 and self.us0 != -1):
      z = np.asarray(z)
      dusdz = 2 * self.wn * self.us0 * np.exp(2 * self.wn * z)
      return (dusdz*np.cos(ag_rad), dusdz*np.sin(ag_rad))
### END Class definition ###
