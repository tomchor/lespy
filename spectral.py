
def spectra_2d(data, simulation=None, what='radial'):
    """
    Calculates 2D spectra on data. Uses last two axes to calculate it

    Parameters
    ----------
    data: ndarray
    simulation: lespy.Simulation
    what: str
        {radial, 2d, integral} for radial average, 2D spectra or integral scale
    """
    import numpy.fft as FFT
    import numpy as np
    from . import utils
    sim=simulation
    S = FFT.fftshift

    nv, nt, nx, ny = data.shape
    data = data-data.mean(axis=(-2,-1), keepdims=True)
    fx=S(np.fft.fftfreq(nx, d=sim.domain.dx))
    fy=S(np.fft.fftfreq(ny, d=sim.domain.dy))
    K1,K2=np.meshgrid(fy,fx)
    Kabs=np.sqrt(K2**2+K1**2)

    fftv = FFT.fft2(data, axes=(-2,-1))
    pspec = Kabs*S((fftv*fftv.conj()).real, axes=(-2,-1))

    if what=='radial':
        radP = []
        for iv in range(nv):
            radK, radp = utils.radial_prof3D(pspec[iv], r=Kabs)
            radP.append(radp)
        return radK, 2.*np.pi*radK[None,None]*np.array(radP)

    elif what=='2d':
        return K1, K2, pspec

    elif what=='integral':
        R = np.tile(Kabs, (nv, nt, 1, 1))
        R = np.ma.masked_where(R==0., R)
        L = np.average(1./R, weights=pspec, axis=(-2,-1))
        return L

    else:
        raise ValueError("'what' keyword wasn't correct. Should be 'radial', '2d' or 'integral'")


