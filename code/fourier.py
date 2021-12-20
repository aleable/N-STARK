"""
N-STARK (Non-STAtionary loads Routing on networKs) -- https://github.com/aleable/N-STARK
Contributors:
    Alessandro Lonardi
    Caterina De Bacco
"""

import numpy as np
from scipy.sparse import csr_matrix
import warnings

warnings.filterwarnings('ignore')


def fs_coeff(f, N, return_complex):
    """Calculates the first 2*N+1 Fourier series coeff. of a periodic function.
    Taken from: https://stackoverflow.com/questions/4258106/how-to-calculate-a-fourier-series-in-numpy

    f(t) ~= a0/2+ sum_{k=1}^{N} ( a_k*cos(2*pi*k*t/T) + b_k*sin(2*pi*k*t/T) )

    If return_complex is set to True, it returns instead the coefficients
    {c0,c1,c2,...}
    such that:

    f(t) ~= sum_{k=-N}^{N} c_k * exp(i*2*pi*k*t/T)

    where we define c_{-n} = complex_conjugate(c_{n})

    Parameters
    ----------
    f : the periodic function, a callable like f(t)
    N : the function will return the first N_max + 1 Fourier coeff.

    Returns
    -------
    if return_complex == False, the function returns:
        a0 : float
        a,b : numpy float arrays describing respectively the cosine and sine coeff.
    if return_complex == True, the function returns:
        c : numpy 1-dimensional complex-valued array of size N+1
    """

    T = 1       # period assumed to be 1
    f_sample = 2 * N
    t, dt = np.linspace(0, T, f_sample + 2, endpoint=False, retstep=True)

    y = np.fft.rfft(f(t)) / t.size

    if return_complex:
        return y
    else:
        y *= 2
        return y[0].real, y[1:-1].real, -y[1:-1].imag


def coefficient_matrix(forcings, N_chosen, return_complex):
    """Compute the Fourier coefficients matrix

    Parameters:
        forcings: function, time vaying forcings
        T: float, period
        N_chosen: int, truncation of Fourier series
        return_complex: bool, complex or real values in coeff matrix

    Returns:
        coeff_matrix: numpy.array(complex), matrix C
    """

    coeff_ = fs_coeff(forcings, N_chosen, return_complex)

    coeff_extended = []
    for row in coeff_:
        ext_row = list(row) + list(row[1:].conj())
        coeff_extended.append(ext_row)
    coeff_extended = np.array(coeff_extended)

    nnode = coeff_.shape[0]
    coeff_matrix = np.zeros((nnode, nnode), dtype=float)
    for i in range(nnode):
        for j in range(nnode):
            coeff_matrix[i][j] = np.dot(coeff_extended[i], coeff_extended[j].conj())

    coeff_matrix = csr_matrix(coeff_matrix)

    return coeff_matrix
