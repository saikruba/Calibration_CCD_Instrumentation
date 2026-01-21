import numpy as np

"""
    General Gaussian model component.

    Parameters:
        E       : Energy array (keV)
        E_center  : Center of the Gaussian peak (keV)
        sigma   : Standard deviation of the peak (keV)
        weight  : Amplitude or weight of the peak

    Returns:
        Gaussian function evaluated over E
"""
    
    
def gaussian(E, E_center, sigma, weight):
    return (weight * np.exp(-((E - E_center) ** 2) / (2 * (sigma ** 2))))


