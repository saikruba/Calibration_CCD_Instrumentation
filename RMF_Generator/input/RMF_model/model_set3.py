import numpy as np
import scipy.integrate
from gaussian_model import gaussian  # Import general Gaussian
   
    
def mormf_set3(E, l, Bt, sig, tau, norm,f0, mu, sig_esc, wt_esc, Ein, sig_fl, wt_fl, sig_noise1, wt_noise1, sig_noise2, wt_noise2):

    """
    Model response function for X-ray RMF set_3 [E= 5.7925 - 10.2475 keV].

    This model simulates the detector response to an incident X-ray photon energy having primary peak, Escape peak, Fluorescence
    peak and two noise peaks. It incorporates attenuation and charge loss effects across two regions transitioning at length l
    from func1 to func2 with added gaussians.
    
    Parameters
    ----------
        E         : channel Energy array (keV) over which the response is evaluated.
        sig       : Gaussian sigma representing energy resolution or Fano noise.
        norm      : Normalization constant to scale the final model output.
        l         : Characteristic length scale (μm)        
        f0        : Fraction of charge collected at x = 0 (initial collection efficiency)?
        Bt        : Beta parameter.
        mu        : Linear attenuation coefficient (1/μm).
        Ein       : Incident photon energy (keV).
        sig_esc   : Standard deviation for the escape peak
        wt_esc    : Weights (intensities) for the escape peak
        Eesc      : Escape photon energy (keV)
        sig_fl    : Standard deviation for the fluorescence peak
        wt_fl     : Weights (intensities) for the fluorescence peak 
        Efl       : Energy of the fluorescence line (keV)
       sig_noise1, 
       sig_noise2 : Standard deviations of the two noise Gaussians.
       wt_noise1,
       wt_noise2  : Weights (intensities) for the two noise Gaussians.        
       Eno1, Eno2 : Energies for low-energy noise components (keV)        

    Returns
    -------
        model    : Detector response N(x) evaluated at each energy E.
    """


    Al = (l * (1 - f0)) / (l + (Bt * tau))
    Gm = ((1 - f0) / (l + (Bt * tau))) * Bt * tau

    def func1(x, E, l, f0, Bt, tau, mu, Ein, sig, Al):
        return np.exp(-mu * x) * np.exp(-((E - (Ein * (f0 + Al * (x / l) ** Bt))) ** 2) / (2 * sig ** 2))

    def func2(x, E, l, f0, Bt, tau, mu, Ein, sig, Gm):
        return np.exp(-mu * x) * np.exp(-((E - (Ein * (1 - Gm * (np.exp(-(x - l) / tau))))) ** 2) / (2 * sig ** 2))


    intgl1, _ = scipy.integrate.quad_vec(func1, 0, l, args=(E, l, f0, Bt, tau, mu, Ein, sig, Al))
    intgl2, _ = scipy.integrate.quad_vec(func2, l, 300, args=(E, l, f0, Bt, tau, mu, Ein, sig, Gm))


    Eesc = Ein - 1.7475
    gauss1 = gaussian(E, Eesc, sig_esc, wt_esc)           # Escape peak
    Efl = 1.7475
    gauss2 = gaussian(E, Efl, sig_fl, wt_fl)              # Fluorescence
    Eno1= 0.17
    gauss3 = gaussian(E, Eno1, sig_noise1, wt_noise1)     # Noise 1
    Eno2= 0.39
    gauss4 = gaussian(E, Eno2, sig_noise2, wt_noise2)     # Noise 2

    # Final model
    model = (intgl1 + intgl2 + gauss1 + gauss2 + gauss3 + gauss4) * (1 / norm)
    return model
    
    


 


