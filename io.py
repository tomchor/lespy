import numpy as _np

def empty_array(simulation, n_con=False):
    """ Creates an empty array with the shape dictated by simulation """
    sim=simulation
    if n_con:
        blank=_np.full((sim.domain.ld, sim.ny, sim.nz_tot, sim.n_con), _np.nan, dtype=_np.float64)
    else:
        blank=_np.full((sim.domain.ld, sim.ny, sim.nz_tot), _np.nan, dtype=_np.float64)
    return blank


def write_to_les(array, fname, simulation=None, **kwargs):
    """
    Writes array into a file fname in a format that LES can easily understand
    """
    sim=simulation
    array[sim.nx:] = 0.
    array.T.tofile(fname, **kwargs)
    return


def read_aver(fname, simulation, squeeze=True, return_times=False, **kwargs):
    """Reads aver_* files from LES"""
    sim=simulation
    aver=_np.loadtxt(fname, **kwargs)
    if 'pcon' in fname.lower():
        aver=aver.reshape(-1, sim.n_con, aver.shape[-1], order='C').transpose(0,2,1)
    ndtimes = aver[:,0]
    aver = aver[:,1:]
    if squeeze:
        aver = _np.squeeze(aver)

    if return_times:
        return ndtimes, aver
    else:
        return aver

