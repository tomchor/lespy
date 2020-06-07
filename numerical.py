import numpy as np

def integrate(arr, axis=0, dx=1., chunks=1):
    """
    Re-implementation of numpy.trapz but saving RAM memory

    axis can't be -1, as opposed to numpy.trapz
    """
    arr = arr.values
    if chunks==1:
        return dx*(0.5*(arr.take(0, axis=axis) + arr.take(-1, axis=axis)) + arr.take(range(1,arr.shape[axis]-1), axis=axis).sum(axis=axis))
    else:
        import numpy as np
        newshape = [ el for index, el in enumerate(arr.shape) if index != axis ]
        integ = np.zeros(newshape)
        limits = np.linspace(0, arr.shape[axis], chunks+1, dtype=int)
        for ini, end in zip(limits,limits[1:]):
            print('Integrating axis {} from {} to {}'.format(axis,ini,end))
            arr2 = arr.take(np.arange(ini, end), axis=axis, mode='clip')
            integ += dx*(0.5*(arr2.take(0, axis=axis) + arr2.take(-1, axis=axis)) + arr2.take(range(1,arr2.shape[axis]-1), axis=axis).sum(axis=axis))
        del arr2
        return integ


def diff_fft_xr(da, n=1, dim=None, shift=False, real=True, reshape=False, **kwargs):
    """
    Calculate derivative of a DataArray using fft

    Notes:
    - Changing the zero frequency doesn't change the output as far
        as the derivative goes.
    - Normalizing the fourier coefficients by N, yields derivatives
        that are 1/N the expected value.
    - The nyquist frequency is generally pretty small so setting it to
        zero generally doesn't do anything.
    """
    import numpy as np
    import xrft
    import xarray as xr

    kdim = 'freq_'+dim

    #----
    # Perform forward DFT
    if real:
        wave = xrft.dft(da, dim=dim, shift=shift, real=dim, **kwargs)
    else:
        wave = xrft.dft(da, dim=dim, shift=shift, real=None, **kwargs)
    freq = wave.coords[kdim]
    #----

    #----
    # Multiply amplitudes by 2πk
    ikF = (pow(1.j*2.*np.pi, n)*wave*freq)
    #----

    #----
    # Perform inverse DFT
    if real:
        deriv = np.fft.irfft(ikF, axis=wave.dims.index(kdim))
    else:
        deriv = np.fft.ifft(ikF, axis=wave.dims.index(kdim)).real
    #----

    #----
    # Repackage into dataarray
    newdims = [ dim if dim!=kdim else kdim.replace("freq_", "") for dim in wave.dims ]
    deriv = xr.DataArray(deriv, coords=da.coords, dims=newdims)
    return deriv.transpose(*da.dims)
    #----



def diff_fft(array, n=1, axis=0, dx=1., zero_nyquist=False):
    """
    Notes:
    - Changing the zero frequency doesn't change the output as far
        as the derivative goes.
    - Normalizing the fourier coefficients by N, yields derivatives
        that are 1/N the expected value.
    - The nyquist frequency is generally pretty small so setting it to
        zero generally doesn't do anything.
    """
    import numpy as np
    wave = np.fft.rfft(array, axis=axis)
    freq = np.fft.rfftfreq(array.shape[axis], d=dx)
    newshape = tuple( 1 if ax!=axis else -1 for ax in range(len(array.shape)) )
    freq = pow(1.j*2.*np.pi, n)*freq.reshape(newshape)

#    if zero_nyquist:
#        wave[-1]=0
    ikF = (wave*freq)
    return np.fft.irfft(ikF, axis=axis)


def ax_replace(arr, indices, val, axis):
    s = [slice(None)]*arr.ndim
    s[axis] = indices
    arr[s] = val


def recover_cont(F, axes=(0,1), coeff=1):
    import numpy as np
    f_c = np.fft.fft2(F, axes=axes)
    for ax in axes:
        nx=F.shape[ax]
        f_c.swapaxes(0, ax)[nx//2] = 0
    #f_c[nx//2]=0.
    #f_c[:,ny//2]=0.
    return np.real(np.fft.ifft2(coeff*f_c, axes=axes))


def _correlate(in1, in2, mode="full", axis=0):
    """
    Convolve two N-dimensional arrays using FFT.
    Convolve `in1` and `in2` using the fast Fourier transform method, with
    the output size determined by the `mode` argument.
    This is generally much faster than `convolve` for large arrays (n > ~500),
    but can be slower when only a few output values are needed, and can only
    output float arrays (int or object array inputs will be cast to float).
    Parameters
    ----------
    in1 : array_like
        First input.
    in2 : array_like
        Second input. Should have the same number of dimensions as `in1`.
        If operating in 'valid' mode, either `in1` or `in2` must be
        at least as large as the other in every dimension.
    mode : str {'full', 'valid', 'same'}, optional
        A string indicating the size of the output:
        ``full``
           The output is the full discrete linear convolution
           of the inputs. (Default)
        ``valid``
           The output consists only of those elements that do not
           rely on the zero-padding.
        ``same``
           The output is the same size as `in1`, centered
           with respect to the 'full' output.
    Returns
    -------
    out : array
        An N-dimensional array containing a subset of the discrete linear
        convolution of `in1` with `in2`.
    """
    from scipy import fftpack
    import numpy as np
    s1 = np.array(in1.shape)
    s2 = np.array(in2.shape)
    in1 = np.asarray(in1)
    #in2 = np.asarray(in2.take(range(in2.shape[axis])[::-1], axis=axis))
    if isinstance(axis,int):
        in2 = np.asarray(np.take(in2, range(s2[axis])[::-1], axis=axis))
    else:
        for ax in axis:
            in2 = np.asarray(np.take(in2, range(s2[ax])[::-1], axis=ax))
    print(in2.shape)

    if in1.ndim == in2.ndim == 0:  # scalar inputs
        return in1 * in2
    elif not in1.ndim == in2.ndim:
        raise ValueError("in1 and in2 should have the same dimensionality")
    elif in1.size == 0 or in2.size == 0:  # empty arrays
        return array([])

    complex_result = (np.issubdtype(in1.dtype, complex) or
                      np.issubdtype(in2.dtype, complex))
    shape = s1
    shape[axis] += s1[axis] - 1
    print('in shape:',s1)
    print('out shape:',shape)

    # Speed up FFT by padding to optimal size for FFTPACK
    fshape = np.array([fftpack.helper.next_fast_len(int(d)) for d in shape])
    fslice = tuple([slice(0, int(sz)) for sz in shape])
    # Pre-1.9 NumPy FFT routines are not threadsafe. For older NumPys, make
    # sure we only call rfftn/irfftn from one thread at a time.
    if axis!=None:
        axis=(axis,)
    print(fshape)
    print(axis)
    #sp1 = np.fft.rfftn(in1, axes=axis)
    sp1 = np.fft.rfftn(in1, s=[fshape[axis]], axes=axis)
    #sp2 = np.fft.rfftn(in2, axes=axis)
    sp2 = np.fft.rfftn(in2, s=[fshape[axis]], axes=axis)
    ret = (np.fft.irfftn(sp1 * sp2, axes=axis)[fslice].copy())
    #ret = (np.fft.irfftn(sp1 * sp2, s=fshape, axes=axis)[fslice].copy())

    if mode == "full":
        return ret
    elif mode == "same":
        return _centered(ret, s1)
    elif mode == "valid":
        return _centered(ret, s1 - s2 + 1)
    else:
        raise ValueError("Acceptable mode flags are 'valid',"
                         " 'same', or 'full'.")


def diff_z(da, n=1, z_final=None):
    """ n-th order derivative according to LES rules """
    #----
    # Get z-dependentr variables
    z = da.z.values
    Δz = abs(np.median(np.diff(z)))
    #----

    #----
    # Get final z recursively
    da_z = da.diff("z", n) / Δz**2
    for ni in range(1,n+1):
        z = (z[1:] + z[:-1])/2
    #----

    #----
    # Re-sets the z corrd if it's good and return deriv
    if z_final is not None:
        if np.allclose(z, z_final):
            da_z = da_z.assign_coords(z=z_final)
        else:
            raise AssertionError
    else:
        da_z = da_z.assign_coords(z=z)
    return da_z
    #----


