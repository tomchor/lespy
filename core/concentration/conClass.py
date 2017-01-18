
class twoDHrzCon:
    """class for 2D horizontal concentration field"""
    #def __init__(self, con=np.zeros(1)):
    def __init__(self, con=None):
        import numpy as np
        if (type(con) is np.ndarray):
            self.con = con
            self.nd = con.shape
        else:
            import sys
            sys.exit("The data type for initializing class \
                    twoDHztCon is not numpy.ndarray!")

    def getcenlin(self, dix=0, diy=0, dx=-1, dy=-1, transpose=False):
        import numpy as np
        """calculate the center line of plume"""
        # the plume must have 2 dimensions
        if (len(self.nd) != 2):
            sys.exit("Cannot calculate center line without 2D field")
        else:
            ix_cmax, iy_cmax = np.where(self.con == self.con.max())
        # the dixection of plume
        if (dix == 0 and diy == 0):
            if (ix_cmax < np.floor(self.nd[0]/2)):
                span_x = ix_cmax
            else:
                span_x = self.nd[0] - 1 - ix_cmax

            if (iy_cmax < np.floor(self.nd[1]/2)):
                span_y = iy_cmax
            else:
                span_y = self.nd[1] - 1 - iy_cmax

            span = np.amin([span_x, span_y, 15])
            ix_check = np.array([ix_cmax-span, ix_cmax+span])
            iy_check = np.array([iy_cmax-span, iy_cmax+span])

            con_dix = self.con[ix_check[0], iy_check[0]]
            for ix in ix_check:
                if (self.con[ix, iy_cmax] > con_dix):
                    con_dix = self.con[ix, iy_cmax]
                    ix_dix = ix
                    iy_diy = iy_cmax
            for iy in iy_check:
                if (self.con[ix_cmax, iy] > con_dix):
                    con_dix = self.con[ix_cmax, iy]
                    ix_dix = ix_cmax
                    iy_diy = iy

            if ix_dix > ix_cmax:
                dix = 1
            elif ix_dix < ix_cmax:
                dix = -1

            if iy_diy > iy_cmax:
                diy = 1
            elif iy_diy < iy_cmax:
                diy = -1

        # search the center line
        if (dix and diy) != 0:
            sys.exit("dix and diy can not be non-zero at the same time!!")
        ix = np.array(ix_cmax)
        iy = np.array(iy_cmax)
        ix_cen = np.array(ix_cmax)
        iy_cen = np.array(iy_cmax)
        con_cen = np.asarray(self.con[ix, iy])
        while True:
            if dix != 0:
                ix += dix
                if ix < 0 or ix >= self.nd[0]:
                    break
                else:
                    ix_cen = np.append(ix_cen, ix)
                    iy_cen = np.append(iy_cen, self.con[ix,:].argmax())
                    con_cen = np.append(con_cen, self.con[ix_cen[-1], iy_cen[-1]])
                    if iy_cen[-1] == self.nd[1]-1:
                        break
            elif diy !=0:
                iy += diy
                if iy < 0 or iy >= self.nd[1]:
                    break
                else:
                    ix_cen = np.append(ix_cen, self.con[:,iy].argmax())
                    iy_cen = np.append(iy_cen, iy)
                    con_cen = np.append(con_cen, self.con[ix_cen[-1], iy_cen[-1]])
                    if ix_cen[-1] == self.nd[0]-1:
                        break

        # save
        self.cenlin = dict()
        if transpose == False:
            self.cenlin.update({'ix':ix_cen})
            self.cenlin.update({'iy':iy_cen})
            self.cenlin.update({'con':con_cen})
            if dx > 0:
                self.cenlin.update({'x':ix_cen*dx})
            if dy > 0:
                self.cenlin.update({'y':iy_cen*dy})
        elif transpose == True:
            self.cenlin.update({'ix':iy_cen})
            self.cenlin.update({'iy':ix_cen})
            self.cenlin.update({'con':con_cen})
            if dy > 0:
                self.cenlin.update({'x':iy_cen*dy})
            if dx > 0:
                self.cenlin.update({'y':ix_cen*dx})

class twoDVtlCon:
    """class for 2D vertiyal concentration field"""
    #def __init__(self, con=np.zeros(1)):
    def __init__(self, con=None):
        if (type(con) is np.ndarray):
            self.con = con
            self.nd = con.shape
        else:
            sys.exit("The data type for initializing class \
                    twoDVtlCon is not numpy.ndarray!")

class threeDCon:
    """class for 3D concentration field"""
    #def __init__(self, con=np.zeros(1)):
    def __init__(self, con=None):
        if (type(con) is np.ndarray):
            self.con = con
            self.nd = con.shape
        else:
            sys.exit("The data type for initializing class \
                    threeDCon is not numpy.ndarray!")

