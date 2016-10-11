
class domain:
  """class for domain parameters"""
  def __init__(self, nx=0, ny=0, nz=0, nz_tot=None, dx=None, dy=None, dz=None, lx=0, ly=0, lz=0):
    self.Nx = nx
    self.Ld = self.Nx + 2
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
      buff+='#'*10 + ' Domain Parameters '+ '#'* 10
      buff+="\nThe grid numbers in 3D are {}, {}, {}\n".format(self.Nx, self.Ny, self.Nz)
      buff+="The grid sizes in 3D are {}, {}, {} m\n".format(self.dx, self.dy, self.dz)
      buff+="The domain sizes in 3D are {}, {}, {} m\n".format(self.Lx, self.Ly, self.Lz)
      buff+='#' * 10+' END Domain Parameters ' + '#'*10
      return buff

  __repr__ = __str__
