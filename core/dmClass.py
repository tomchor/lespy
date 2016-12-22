
class domain:
  """Class for domain parameters"""
  def __init__(self, nx=None, ny=None, nz=None, nz_tot=None,
          dx=None, dy=None, dz=None, lx=None, ly=None, lz=None):
    self.nx = nx
    self.ld = 2*((self.nx//2)+1)
    self.ny = ny
    self.nz = nz
    self.nz_tot = nz_tot
    self.lx = lx
    self.ly = ly
    self.lz = lz

    if dx!=None:
        self.dx = dx
    else:
        try:
            self.dx = self.lx/self.nx
        except:
            self.dx = dx
    if dy!=None:
        self.dy = dy
    else:
        try:
            self.dy = self.ly/self.ny
        except:
            self.dy = dy
    if dz!=None:
        self.dz = dz
    else:
        try:
            self.dz = self.lz/self.nz
        except:
            self.dz = dz

  def __str__(self):
      buff = ''
      buff+="\nx, y, z grid numbers: {}, {}, {}\n".format(self.nx, self.ny, self.nz)
      buff+="dx, dy, dz: {}, {}, {} m\n".format(self.dx, self.dy, self.dz)
      buff+="lx, ly, lz: {}, {}, {} m\n".format(self.lx, self.ly, self.lz)
      return buff

  __repr__ = __str__
