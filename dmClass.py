#!/usr/bin/env python3.4

class domain:
  """class for domain parameters"""
  def __init__(self, nx=0, ny=0, nz=0, dx=None, dy=None, dz=None, lx=0, ly=0, lz=0):
    self.Nx = nx
    self.Ld = self.Nx + 2
    self.Ny = ny
    self.Nz = nz
    self.Lx = lx
    self.Ly = ly
    self.Lz = lz

    if dx!=None:
        self.Dx = dx
    else:
        try:
            self.Dx = self.Lx/self.Nx
        except:
            self.Dx = dx
    if dy!=None:
        self.Dy = dy
    else:
        try:
            self.Dy = self.Ly/self.Ny
        except:
            self.Dy = dy
    if dz!=None:
        self.Dz = dz
    else:
        try:
            self.Dz = self.Lz/self.Nz
        except:
            self.Dz = dz

  def __str__(self):
    print('#' * 10,' Domain Parameters ','#' * 10)
    print("The grid numbers in 3D are ",self.Nx, self.Ny, self.Nz)
    print("The grid sizes in 3D are ",self.Dx, self.Dy, self.Dz, " m")
    print("The domain sizes in 3D are ",self.Lx, self.Ly, self.Lz, " m")
    print('#' * 10,' END Domain Parameters ','#' * 10)

  __repr__ = __str__
