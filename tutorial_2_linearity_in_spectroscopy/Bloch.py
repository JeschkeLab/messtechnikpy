import numpy as np
from numpy import pi

def Bloch_stationary(frq_axis,frq,T1,T2,nu1 = 0.002):
    """
    Computes spectrum from specification of a single spectral line and an
    optional frequency axis based on stationary solutions of the Bloch
    equations
    an appropriate frequency axis is determined, if none is provided

    Input parameters:
        frq    frequency offset from carrier (resonance offset)
        T1     longitudinal relaxation time
        T2     transverse relaxation time
    Optional parameters:
        nu1    irradiation amplitude (default=0.002, linear regime)
    Output parameters:
    frq_axis  frequency axis
    spectrum  complex spectrum, real part absorption, imaginary part dispersion

    Luis Fabregas, 2020 adapted from G. Jeschke, 2011, for lecture course Messtechnik
    """ 

    ndefault = 1024 # default number of points of frequency axis
    verbose = True  # warnings are output, if true and suppressed, if false

    # create frequency axis, if necessary
    linewidth = 1/(pi*T2)

    # determine if irradiation amplitude leads to saturation
    S = (2*pi*nu1)**2*T1*T2

    # determine range of significant line intensity
    minfrq = frq - 10*linewidth
    maxfrq = frq + 10*linewidth

    if min(frq_axis)>minfrq and verbose:
        print('Warning: Spectral line extends beyond minimum of frequency axis.\n')

    if max(frq_axis)<maxfrq and verbose:
        print('Warning: Spectral line extends beyond maximum of frequency axis.\n')


    arg = (2*pi*(frq_axis - frq)*T2)**2
    denom = arg + 1 + S
    dispersion = -2*pi*nu1*(2*pi*(frq_axis - frq)*T2**2)/denom
    absorption = 2*pi*nu1*T2*np.ones_like(frq_axis)/denom

    spectrum = absorption + 1j*dispersion

    return spectrum