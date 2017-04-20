
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

    print('Calculating from phi(x,y) condnorm to phi(r) ... ', end='')
    rCond = [ 
    [ utils.radial_prof(condC[iv,it,:,:], simulation=sim, axes=(0,1))[1] for it in range(timelength) ]
       for iv in range(condC.shape[0]) ]
    rCond = np.asarray(rCond)
    Rs = utils.radial_prof(condC[0,0], simulation=sim, axes=(0,1))[0]
    print('done')
    return Rs, rCond

def power_law(r, rc, gamma):
    """ Theoretical power law for the shape of the normalized conditional density """
    out = (r/rc)**gamma
    out[ out<=1. ] = 1.
    return out

def fromCond3(Rs, rnorm):
    from scipy.optimize import curve_fit
    from matplotlib import pyplot as plt
    Rs=Rs[1:]
    rnorm=rnorm[:,1:]
    rcs = []
    gammas = []
    for rnm0 in rnorm:
        rc, gamma = curve_fit(power_law, Rs, rnm0, p0=(100., -1.), maxfev=1000)[0]
        rcs.append(rc)
        gammas.append(gamma)
        if 0:
            plt.close()
            plt.loglog(Rs, rnm0)
            plt.loglog(Rs, power_law(Rs, rc, gamma))
            plt.show()
        #print(rc, gamma)
    return rcs, gammas


