
class Domain(object):
    """Class for domain parameters"""
    def __init__(self, nx=None, ny=None, nz=None, nz_tot=None, 
                 dx=None, dy=None, dz=None, lx=None, ly=None, lz=None, 
                 origin=(0,0,0), ocean_flag=False):
        """Initialize the class. dx, dy, and dz can be calculated"""
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.nz_tot = nz_tot
        self.ld = 2*((self.nx//2)+1)
        self.origin_node = origin

        self.lx = lx
        self.ly = ly
        self.lz = lz

        self.dx = dx
        self.dy = dy
        self.dz = dz

        self.ocean_flag = ocean_flag

        self.points = self.nx*self.ny*self.nz
        self.get_resolution()
        self.get_delta()
        self.makeAxes()

    def get_resolution(self):
        """gets the x, y, z resolution if it isn't given"""
        if self.dx==None:
            try:
                self.dx = self.lx/self.nx
            except:
                pass
        if self.dy==None:
            try:
                self.dy = self.ly/self.ny
            except:
                pass
        if self.dz==None:
            try:
                self.dz = self.lz/self.nz
            except:
                pass
        return
    
    def get_delta(self):
        """Gets the characteristic volume box delta"""
        self.delta = (self.dx*self.dy*self.dz)**(1./3.)
        return
    

    def makeAxes(self, array=None, origin=None):
        """
        Creates the x, y and z axes based on dx and nx
        
        Parameters
        ----------
        array: numpy.ndarray
            array from which to get the number of points (for endless)
        origin: tuple, list
            3D coord of node that is to be taken as origin. Default is
            whatever is in the domain object. Ex.: origin=(50,50,50)
        """
        import numpy as np

        if origin:
            origin_node = origin
        else:
            origin_node = self.origin_node

        if type(array) != type(None):
            if len(array.shape)==3:
                nx, ny, nz_tot = array.shape
            elif len(array.shape)==2:
                nx, ny = array.shape
                nz_tot = None
            else:
                raise ValueError("makeAxes() can only work with (x,y)- or (x,y,z)-shaped arrays")
        else:
            nx = self.nx
            ny = self.ny
            nz_tot = self.nz_tot

        x = np.arange(0, nx)*self.dx - origin_node[0]*self.dx
        y = np.arange(0, ny)*self.dy - origin_node[1]*self.dy
        if self.ocean_flag:
            z_w = -(np.arange(0, nz_tot) - origin_node[2])*self.dz
            z_u = -(np.arange(0, nz_tot) - origin_node[2] + 1/2)*self.dz
        else:
            z_w = (np.arange(0, nz_tot) - origin_node[2])*self.dz
            z_u = (np.arange(0, nz_tot) - origin_node[2] + 1/2)*self.dz

        if type(array) != type(None):
            return x,y,z
        else:
            self.x = x+self.dx/2
            self.y = y+self.dy/2
            self.z = z_u
            self.z_u = z_u
            self.z_w = z_w
            return


    def __str__(self):
        buff = """<lespy.domain object>
nx, ny, nz: {self.nx} x {self.ny} x {self.nz} = {self.points} points
dx, dy, dz: {self.dx:.2f} x {self.dy:.2f} x {self.dz:.2f}
lx, ly, lz: {self.lx} x {self.ly} x {self.lz}""".format(**locals())
        return buff

    def __repr__(self):
        buff = """<lespy.domain with {}x{}x{} points>""".format(self.nx, self.ny, self.nz)
        return buff



class hDomain(object):
    """Class for domain parameters"""
    def __init__(self, nx=None, ny=None, nz=None, 
                 dx=None, dy=None, dz=None, lx=None, ly=None, lz=None, 
                 origin=(0,0,0), environment=False):
        """Initialize the class. dx, dy, and dz can be calculated"""
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.nz_tot = nz
        self.ld = 2*((self.nx//2)+1)
        self.origin_node = origin

        self.lx= lx
        self.ly= ly
        self.lz= lz

        self.dx = dx
        self.dy = dy
        self.dz = dz

        self.environment = environment

        self.points = self.nx*self.ny*self.nz
        self.get_resolution()
        self.get_delta()
        self.makeAxes()

    def get_resolution(self):
        """gets the x, y, z resolution if it isn't given"""
        if self.dx==None:
            try:
                self.dx = self.lx/self.nx
            except:
                pass
        if self.dy==None:
            try:
                self.dy = self.ly/self.ny
            except:
                pass
        if self.dz==None:
            try:
                self.dz = self.lz/self.nz
            except:
                pass
        return
    
    def get_delta(self):
        """Gets the characteristic volume box delta"""
        self.delta = (self.dx*self.dy*self.dz)**(1./3.)
        return
    

    def makeAxes(self, array=None, origin=None):
        """
        Creates the x, y and z axes based on dx and nx
        
        Parameters
        ----------
        array: numpy.ndarray
            array from which to get the number of points (for endless)
        origin: tuple, list
            3D coord of node that is to be taken as origin. Default is
            whatever is in the domain object. Ex.: origin=(50,50,50)
        """
        import numpy as np

        if origin:
            origin_node = origin
        else:
            origin_node = self.origin_node

        if type(array) != type(None):
            if len(array.shape)==3:
                nx, ny, nz_tot = array.shape
            elif len(array.shape)==2:
                nx, ny = array.shape
                nz_tot = None
            else:
                raise ValueError("makeAxes() can only work with (x,y)- or (x,y,z)-shaped arrays")
        else:
            nx = self.nx
            ny = self.ny
            nz_tot = self.nz_tot

        x = np.arange(0, nx)*self.dx - origin_node[0]*self.dx
        y = np.arange(0, ny)*self.dy - origin_node[1]*self.dy
        z_w = (np.arange(0, nz_tot) - origin_node[2])*self.dz
        z_u = (np.arange(0, nz_tot) - origin_node[2] + 1/2)*self.dz
        if self.environment=="ocean":
            z_w -= z_w
            z_u -= z_u

        if type(array) != type(None):
            return x,y,z
        else:
            self.x = x+self.dx/2
            self.y = y+self.dy/2
            self.z = z_u
            self.z_u = z_u
            self.z_w = z_w
            return


    def __str__(self):
        buff = """<lespy.domain object>
nx, ny, nz: {self.nx} x {self.ny} x {self.nz} = {self.points} points
dx, dy, dz: {self.dx:.2f} x {self.dy:.2f} x {self.dz:.2f}
lx, ly, lz: {self.lx} x {self.ly} x {self.lz}""".format(**locals())
        return buff

    def __repr__(self):
        buff = """<lespy.domain with {}x{}x{} points>""".format(self.nx, self.ny, self.nz)
        return buff


