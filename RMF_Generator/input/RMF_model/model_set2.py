import numpy as np
import scipy.integrate
from gaussian_model import gaussian


def mormf_set2(E, l, Bt, sig, tau, norm, f0, mu, sig_esc, wt_esc, Ein):

    """
    Model response function for X-ray RMF set_2 [E= 1.8425 - 5.7875 keV].

    This model simulates the detector response to an incident X-ray photon energy having primary peak and escape peak,
    incorporating attenuation and charge loss effects across two regions transitioning at length l from func1 to func2
    with added gaussian for the escape peak.

    Parameters
    ----------
        E        : channel Energy array (keV) over which the response is evaluated.
        sig      : Gaussian sigma representing energy resolution or Fano noise.
        norm     : Normalization constant to scale the final model output.
        l        : Characteristic length scale (μm)        
        f0       : Fraction of charge collected at x = 0 (initial collection efficiency)?
        Bt       : Beta parameter.
        mu       : Linear attenuation coefficient (1/μm).
        Ein      : Incident photon energy (keV).
        sig_esc  : Standard deviation for the escape peak
        wt_esc   : Weights (intensities) for the escape peak
        Eesc     : Escape photon energy (keV)

    Returns
    -------
        model    : Detector response N(x) evaluated at each energy E.
    """


    # Calculate Al and Gm constants
    Al = (l * (1 - f0)) / (l + (Bt * tau))
    Gm = ((1 - f0) / (l + (Bt * tau))) * Bt * tau

    # Define integrand functions
    def func1(x, E, l, f0, Bt, tau, mu, Ein, sig, Al):
        return np.exp(-mu * x) * np.exp(-((E - (Ein * (f0 + Al * (x / l) ** Bt))) ** 2) / (2 * (sig ** 2)))

    def func2(x, E, l, f0, Bt, tau, mu, Ein, sig, Gm):
        return np.exp(-mu * x) * np.exp(-((E - (Ein * (1 - Gm * (np.exp(-(x - l) / tau))))) ** 2) / (2 * (sig ** 2)))

    
    # Integrate with vectorized input over two regions
    intgl1, _ = scipy.integrate.quad_vec(func1, 0, l, args=(E, l, f0, Bt, tau, mu, Ein, sig, Al))
    intgl2, _ = scipy.integrate.quad_vec(func2, l, 300, args=(E, l, f0, Bt, tau, mu, Ein, sig, Gm))

    # Use external Gaussian function for escape peak
    Eesc = Ein - 1.7475
    gauss1 = gaussian(E, Eesc, sig_esc, wt_esc)
    
    # Model output
    model = (intgl1 + intgl2 + gauss1) * (1 / norm)
    return model






 


