
def radial_dist(data, cond0=lambda x, y: x<=np.percentile(y,5), condr=lambda x, y: x>=np.percentile(y,95), simulation=None, bins=None):
    """
    data: np.ndarray
        indices are:
            0: time
            1: x
            2: y
    """
    import numpy as np
    sim=simulation
    nt, Lx1, Ly1 = data.shape
    Lx, Ly = (np.array([Lx1, Ly1])/2).astype(int)
    x = np.arange(-Lx,-Lx+Lx1,1)*sim.domain.dx
    y = np.arange(-Ly,-Ly+Ly1,1)*sim.domain.dy
    xx, yy = np.meshgrid(x, y)
    r = np.sqrt(xx**2. + yy**2.)
    print(r.argmin())
    exit()

    if type(bins)==type(None):
        bins = np.arange(0, 700, 10)

    x = np.arange(-Lx,-Lx+Lx1,1)
    y = np.arange(-Ly,-Ly+Ly1,1)

    full_hist = np.zeros((nt, Lx1, Ly1, len(bins)-1))
    for it in range(nt):
        origins = np.where(cond0(data[it], data[it]))
        for ix,iy in zip(*origins):
            rolled = np.roll(data[it], -x[ix], axis=0)
            rolled = np.roll(rolled,-y[iy],axis=1)
            high_r = r[ condr(rolled, rolled) ]
            full_hist[it,ix,iy,:] = np.histogram(high_r, bins=bins)[0]
    hist = full_hist.mean(axis=(1,2))
    summ = hist.sum(axis=(1), keepdims=True)
    summ[ summ==0 ] = 1.
    hist = hist/summ
    norm = np.histogram(r, bins=bins)[0]
    hist = hist/norm
    centers = (bins[:-1]+bins[1:])/2
    return hist, centers


def fromRadial(Hist, bins, window=None):
    """ Gets estimate if L from radial distr function """
    from . import vector
    import numpy as np
    maxima=[]
    if window:
        Hist = vector.moving_average(Hist, axis=1, window=window)
        bins = bins[:(1-window)]
    for hist0 in Hist:
        aux = vector.detect_local_minima(-hist0)
        maxima.append(aux[0][0])
    return bins[np.array(maxima)]


def radial_cdensity(Vars, simulation=None, nc=None, func=None):
    """ Calculates the normalized conditional density as a function of radius """
    from . import vector, utils
    import numpy as np
    sim=simulation
    timelength=Vars.shape[1]

    if type(func)==type(None):
        func=vector.condnorm2d_fft
    if type(nc)==type(None):
        nc=sim.nx//2

    try:
        x, y, condC = func(Vars, simulation=sim, nxc=nc, nyc=nc)
    except TypeError:
        x, y, condC = func(Vars, simulation=sim)
    print(condC.shape)
    nv, nt, nnx, nny = condC.shape

    print('Calculating from phi(x,y) condnorm to phi(r) ... ')
    Rs = utils.radial_prof(condC[0,0], simulation=sim, axes=(0,1))[0]
    rCond = np.zeros((nv, nt, Rs.shape[0]), dtype=np.float64)
    for iv in range(nv):
        print('For variable {} of {}'.format(iv+1,nv))
        rCond[iv,:,:] = utils.radial_prof3D(condC[iv], simulation=sim)[1]
#    for iv in range(nv):
#        print('For variable {} of {}'.format(iv+1,nv))
#        for it in range(nt):
#            rCond[iv, it, :] = utils.radial_prof(condC[iv,it,:,:], simulation=sim, axes=(0,1))[1]
    print('done')
    return np.array(Rs), np.array(rCond)

def power_law(r, rc, gamma):
    """ Theoretical power law for the shape of the normalized conditional density """
    out = (r/rc)**gamma
    out[ out<=1. ] = 1.
    return out

def fromCond4(Rs, rNorm, p0=(100., -1)):
    """
    Must be a 3D array with index 2 being radian conditional density and
    index 0 and 1 being whatever
    """
    import numpy as np
    from scipy.optimize import curve_fit
    Rs=Rs[1:]
    rNorm=rNorm[:,:,1:]
    Lcs = []
    Gcs = []
    for iv, rnm0 in enumerate(rNorm):
        Lcs.append([])
        Gcs.append([])
        for it, rnm1 in enumerate(rnm0):
            rc, gamma = curve_fit(power_law, Rs, rnm1, p0=p0, maxfev=1000)[0]
            Lcs[iv].append(rc)
            Gcs[iv].append(gamma)
            #------
            # In case we want to see what's happening
            if 0:
                from matplotlib import pyplot as plt
                plt.close()
                plt.loglog(Rs, rnm0)
                plt.loglog(Rs, power_law(Rs, rc, gamma))
                plt.show()
            #print(rc, gamma)
            #------
    return np.array(Lcs), np.array(Gcs)


def fromCond3(Rs, rnorm, p0=(10., -1)):
    """
    Must be a 2D array with indexes 1 being radian conditional density and
    index 0 being whatever
    """
    from scipy.optimize import curve_fit
    from matplotlib import pyplot as plt
    Rs=Rs[1:]
    rnorm=rnorm[:,1:]
    rcs = []
    gammas = []
    for rnm0 in rnorm:
        rc, gamma = curve_fit(power_law, Rs, rnm0, p0=p0, maxfev=1000)[0]
        rcs.append(rc)
        gammas.append(gamma)
        #------
        # In case we want to see what's happening
        if 0:
            plt.close()
            plt.loglog(Rs, rnm0)
            plt.loglog(Rs, power_law(Rs, rc, gamma))
            plt.show()
        #print(rc, gamma)
    return rcs, gammas


def get_L(Vars, simulation=None, func=None, p0=(10, -1)):
    """
    Vars should ve 4D, with x, y being the last two dimensions
    """
    import numpy as np
    from . import vector
    sim=simulation
    if func==None:
        func=vector.condnorm2d_fft
    Rs, Phi = radial_cdensity(Vars, simulation=sim, func=func)

#    Gammas=[]
#    Ls=[]
#    for iv in range(Phi.shape[0]):
#        Ls_rc, gammas = fromCond3(Rs, Phi[iv])
#        Gammas.append(gammas)
#        Ls.append(Ls_rc)
#    Gammas=np.array(Gammas)
#    Ls=np.array(Ls)

    Ls, Gammas = fromCond4(Rs, Phi, p0=p0)

    return Ls, Gammas

def from_2Dfft_old(R1, R2, P, cut=None):
    import numpy as np
    PP = P.flatten()
    C=int(PP.shape[0]/2)
    PP=PP[C:]
    r2=R2.flatten()[C:]
    r1=R1.flatten()[C:]
    ID=PP.argmax()
    u0=r1[ID]
    v0=r2[ID]
    Lx=abs(2./(np.sqrt(3.)*u0))
    Ly=abs(2./(3.*v0))
    return np.mean([Lx, Ly])


def from_2Dfft(R1, R2, P, cut=250):
    import numpy as np
    R=np.sqrt(R1**2+R2**2)
    K=np.average(R, weights=P)
    return 1./K


def L_from_fft2(data, simulation=None, what='L'):
    """ """
    import numpy.fft as FFT
    import numpy as np
    from . import utils
    sim=simulation
    S = FFT.fftshift

    nv, nt, nx, ny = data.shape
    data = data-data.mean(axis=(-2,-1), keepdims=True)
    fx=S(np.fft.fftfreq(nx, d=sim.domain.dx))
    fy=S(np.fft.fftfreq(ny, d=sim.domain.dy))
    R2,R1=np.meshgrid(fy,fx)
    R0=np.sqrt(R1**2+R2**2)

    fftv = FFT.fft2(data, axes=(-2,-1))
    pspec = R0*S((fftv*fftv.conj()).real, axes=(-2,-1))

    if what=='radial':
        radR, radP = [[]]*nv, [[]]*nv
        for iv in range(nv):
            for it in range(nt):
                radr, radp = utils.radial_prof(pspec[iv, it], r=R0)
                radR[iv].append(radr)
                radP[iv].append(radp)
        return np.array(radR), np.array(radP)

    if what=='bi':
        return R1, R2, pspec

    if what=='integral':
        R = np.tile(R0, (nv, nt, 1, 1))
        R = np.ma.masked_where(R==0., R)
        L = np.average(1./R, weights=pspec, axis=(-2,-1))
        return L

    if what=='max':
        peak = pspec.reshape(nv, nt, -1).argmax(axis=-1)
        K1 = R1.flatten()[peak]
        K2 = R1.flatten()[peak]
        return 1./np.sqrt(K1**2 + K2**2)

    else:
        R = np.tile(R0, (nv, nt, 1, 1))
        K = np.average(R, weights=pspec, axis=(-2,-1))
        return 1./K

def L_from_fft(data, simulation=None, what='L'):
    """ """
    import numpy.fft as FFT
    import numpy as np
    from . import utils
    sim=simulation
    S = FFT.fftshift

    nv, nt, nx, ny = data.shape
    data = data-data.mean(axis=(-2,-1), keepdims=True)
    fx=S(np.fft.fftfreq(nx, d=sim.domain.dx))
    fy=S(np.fft.fftfreq(ny, d=sim.domain.dy))
    R2,R1=np.meshgrid(fy,fx)
    R0=np.sqrt(R1**2+R2**2)

    fftv = FFT.fft2(data, axes=(-2,-1))
    pspec = S((fftv*fftv.conj()).real, axes=(-2,-1))

    if what=='radial':
        radR, radP = [[]]*nv, [[]]*nv
        for iv in range(nv):
            for it in range(nt):
                radr, radp = utils.radial_prof(pspec[iv, it], r=R0)
                radR[iv].append(radr)
                radP[iv].append(radp)
        return np.array(radR), np.array(radP)

    if what=='bi':
        return R1, R2, pspec

    if what=='integral':
        R = np.tile(R0, (nv, nt, 1, 1))
        R = np.ma.masked_where(R==0., R)
        L = np.average(1./R, weights=pspec, axis=(-2,-1))
        return L

    if what=='max':
        peak = pspec.reshape(nv, nt, -1).argmax(axis=-1)
        K1 = R1.flatten()[peak]
        K2 = R1.flatten()[peak]
        return 1./np.sqrt(K1**2 + K2**2)

    else:
        R = np.tile(R0, (nv, nt, 1, 1))
        K = np.average(R, weights=pspec, axis=(-2,-1))
        return 1./K

