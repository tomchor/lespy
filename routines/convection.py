
def radial_homogFunction(Vars, simulation=None, nc=None, func=None):
    """ Calculates the normalized conditional density as a function of radius """
    from .. import stats, utils
    import numpy as np
    sim=simulation
    timelength=Vars.shape[1]

    if type(func)==type(None):
        func=stats.condnorm2d_fft
    if type(nc)==type(None):
        nc=sim.nx//2

    x, y, condC = func(Vars, simulation=sim)
    nv, nt, nnx, nny = condC.shape

    print('Calculating phi(r) from phi(x,y) ... ')
    Rs = utils.radial_prof(condC[0,0], simulation=sim, axes=(0,1))[0]
    rCond = np.zeros((nv, nt, Rs.shape[0]), dtype=np.float64)
    for iv in range(nv):
        print('For variable {} of {}'.format(iv+1,nv))
        rCond[iv,:,:] = utils.radial_prof3D(condC[iv], simulation=sim)[1]
    print('done')
    return np.array(Rs), np.array(rCond)



def power_law(r, rc, gamma):
    """ Theoretical power law for the shape of the normalized conditional density """
    import numpy as np
    import warnings
    warnings.filterwarnings('error')
    try:
        return np.piecewise(r, [r < rc, r >= rc], [lambda x: (x/rc)**-gamma, lambda x: 1])
    except Exception as e:
        print(rc, gamma)
        print(e)
        raise e




def from_Phi(Rs, rNorm):
    """
    Must be a 3D array with index 2 being radian conditional density and
    index 0 and 1 being whatever
    """
    import numpy as np
    from scipy.optimize import curve_fit
    import warnings
    warnings.filterwarnings('error')
    Rs=Rs[1:]
    lX=np.log(Rs)
    rNorm=rNorm[:,:,1:]
    Lcs = []
    Gcs = []
    #-----
    # Iterate in variable
    for iv, rnm0 in enumerate(rNorm):
        print('Curve-fitting for variable {} of {}.'.format(iv+1, len(rNorm)))

        lY=np.log(rnm0)
        A = np.vstack([lX, np.ones(len(lX))]).T
        gamma, c = np.linalg.lstsq(A, lY.T)[0]
        L=np.exp(-c/gamma)

        Lcs.append([L])
        Gcs.append([gamma])
    #-----
    return np.array(Lcs), np.array(Gcs)


def fromCond5(Rs, rNorm, p0=(30., 1e-1), maxG=10, upperlim=100):
    """
    Must be a 3D array with index 2 being radian conditional density and
    index 0 and 1 being whatever
    """
    import numpy as np
    from scipy.optimize import curve_fit
    from ..utils import nearest
    import warnings
    warnings.filterwarnings('error')

    #-------
    # Remove point where it's zero
    Rs=Rs[1:]
    rNorm=rNorm[:,:,1:]
    #-------

    #--------
    # Optimization:
    xmax=nearest(Rs, upperlim, return_idx=True)[0]
#    stride = int(max(1, xmax / 10))

    X_sparse=Rs[0:xmax]
    Y_sparse=rNorm[:,:,0:xmax]
    #--------

    Lcs = []
    Gcs = []
    #-----
    # Iterate in variable
    for iv, rnm0 in enumerate(rNorm):
        print('Curve-fitting for variable {} of {}.'.format(iv+1, len(rNorm)))
        Lcs.append([])
        Gcs.append([])
        #-----
        # Iterate in time
        for it, rnm1 in enumerate(rnm0):
            p1 = curve_fit(power_law, X_sparse, Y_sparse[iv,it], p0=p0, check_finite=False, 
                    bounds=([1e-5,0], [np.inf,maxG]))[0]

            rc, gamma = curve_fit(power_law, Rs, rnm1, p0=p1, check_finite=False, 
                    bounds=([1e-5,0], [np.inf,maxG]))[0]
            Lcs[iv].append(rc)
            Gcs[iv].append(gamma)
            #------
        #-----
    #-----
    return np.array(Lcs), np.array(Gcs)



def fromCond4(Rs, rNorm, p0=(30., 1e-1), maxG=10):
    """
    Must be a 3D array with index 2 being radian conditional density and
    index 0 and 1 being whatever
    """
    import numpy as np
    from scipy.optimize import curve_fit
    import warnings
    warnings.filterwarnings('error')
    Rs=Rs[1:]
    rNorm=rNorm[:,:,1:]
    Lcs = []
    Gcs = []
    #-----
    # Iterate in variable
    for iv, rnm0 in enumerate(rNorm):
        print('Curve-fitting for variable {} of {}.'.format(iv+1, len(rNorm)))
        Lcs.append([])
        Gcs.append([])
        #-----
        # Iterate in time
        for it, rnm1 in enumerate(rnm0):
            rc, gamma = curve_fit(power_law, Rs, rnm1, p0=p0, check_finite=False, 
                    bounds=([1e-5,0], [np.inf,maxG]))[0]
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
        #-----
    #-----
    return np.array(Lcs), np.array(Gcs)


def normalizeVars(Vars):
    """
    Normalize Vars to be used for estimating L
    Vars should be 4D, with x, y being the last two dimensions
    """
    import numpy as np
    STD=Vars.std(axis=(2,3), keepdims=True)

    #------
    # Prepare for the case where STD is zero
    mask=(STD==0)
    if mask.any():
        STD=np.ma.array(STD, mask=mask)
    #------

    stVars=(Vars-Vars.min(axis=(2,3), keepdims=True))/STD
    return stVars


def get_L(Vars, simulation=None, p0=(30, 1e-1), maxG=10, 
        return_phi=False, pre_process=True, method=4):
    """
    Vars should be 4D, with x, y being the last two dimensions

    If std at some time for some variable is zero, it is "fixed" by creating a mask
    over those occasions. After the calculations are performed the masked places
    receive a NaN value.

    Vars: ndarray
        array from which to get L
    simulation: lespy.Simulation
        context
    p0: list, array, float
        for curve_fit
    maxG: float
        upper bound for Gamma in the fit (prevents some errors).
    return_phi: bool
        whether to return Phi(r) or return L, Gamma
    pre_process: bool
        whether to first pre-process the data by X=(x-x.min())/x.std()
    """
    import numpy as np
    from .. import stats
    sim=simulation

    #-------
    # We normalize the data as pre-proc to make them range between [0:~10]
    if pre_process:
        Vars=normalizeVars(Vars)
    #-------

    #-------
    # We obtain the homogeneity function using the fastest method
    Rs, Phi = radial_homogFunction(Vars, simulation=sim)
    #-------

    #-------
    if np.ma.is_masked(Vars):
        mask2D=Vars.mask.any(axis=(2,3))
    #-------

    #-------
    # We can return Phi for checking or return the scales (default)
    if return_phi:
        #------
        # We mask Phi if Vars is masked because the creation of Phi ignores masks
        if np.ma.is_masked(Vars):
            mask3D=np.tile(mask2D[:,:,None], (1,1,len(Rs)))
            Phi=np.ma.array(Phi, mask=mask3D)
        #------

        return Rs, Phi

    else:
        if method==4:
            Ls, Gammas = fromCond4(Rs, Phi, p0=p0, maxG=maxG)
        elif method==5:
            Ls, Gammas = fromCond5(Rs, Phi, p0=p0, maxG=maxG)

        #------
        # We put nans where there's a mask in Vars because the curve_fit ignores masks
        if np.ma.is_masked(Vars):    
            Ls[mask2D]=np.nan
            Gammas[mask2D]=np.nan
        #------

        return 2.*Ls, Gammas
    #-------


def _fromRadial(Hist, bins, window=None):
    """ Gets estimate if L from radial distr function """
    from .. import stats, utils
    import numpy as np
    maxima=[]
    if window:
        Hist = stats.moving_average(Hist, axis=1, window=window)
        bins = bins[:(1-window)]
    for hist0 in Hist:
        aux = utils.detect_local_minima(-hist0)
        maxima.append(aux[0][0])
    return bins[np.array(maxima)]



def _from_2Dfft(R1, R2, P, cut=250):
    import numpy as np
    R=np.sqrt(R1**2+R2**2)
    K=np.average(R, weights=P)
    return 1./K


def _L_from_fft2(data, simulation=None, what='L'):
    """ """
    import numpy.fft as FFT
    import numpy as np
    from .. import utils
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

def _power_law_old(r, rc, gamma):
    """ Theoretical power law for the shape of the normalized conditional density """
    import warnings
    warnings.filterwarnings('error')
    try:
        out = (r/rc)**(-gamma)
        out[ out<=1. ] = 1.
    except Exception as e:
        print(rc, gamma)
        print(e)
        raise e
    return out


