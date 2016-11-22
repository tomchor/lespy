
class domain:
  """Class for domain parameters"""
  def __init__(self, nx=None, ny=None, nz=None, nz_tot=None,
          dx=None, dy=None, dz=None, lx=None, ly=None, lz=None):
    self.Nx = nx
    self.Ld = 2*((self.Nx//2)+1)
    self.Ny = ny
    self.Nz = nz
    self.Nz_tot = nz_tot
    self.Lx = lx
    self.Ly = ly
    self.Lz = lz

    if dx!=None:
        self.dx = dx
    else:
        try:
            self.dx = self.Lx/self.Nx
        except:
            self.dx = dx
    if dy!=None:
        self.dy = dy
    else:
        try:
            self.dy = self.Ly/self.Ny
        except:
            self.dy = dy
    if dz!=None:
        self.dz = dz
    else:
        try:
            self.dz = self.Lz/self.Nz
        except:
            self.dz = dz

  def __str__(self):
      buff = ''
      buff+="\nx, y, z grid numbers: {}, {}, {}\n".format(self.Nx, self.Ny, self.Nz)
      buff+="dx, dy, dz: {}, {}, {} m\n".format(self.dx, self.dy, self.dz)
      buff+="Lx, Ly, Lz: {}, {}, {} m\n".format(self.Lx, self.Ly, self.Lz)
      return buff

  __repr__ = __str__
