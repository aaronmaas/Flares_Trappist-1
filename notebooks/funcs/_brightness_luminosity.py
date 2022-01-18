import numpy as np


def brightness_star(spectrum_star,spectrumerror_star, T_star, Terror_star, limit, limiterror, model = "None"):
    '''
    Apparent Brightness calculater for star.

    Parameters
    -------------
    spectrum_star : 1-d array
        The spectrum of the star given as a fluxdensity on top of Earths atmosphere without 
        atmospheric extinction [Erg/s/cm^2/A]
    spectrumerror_star: 1-d array
        The error of the spectrum of the star given as error of the fluxdensity on top of the atmsophere without 
        atmospheric extinction [Erg/s/cm^2/A]
    T_star: float
        Blackbody temperature of the star [K]
    Terror_star: float
        Blackbody temperature uncertainty of the star [K]
    limit: float
        Indices for limits in the spectrum_star array (for slicing in the passbands)
    limiterror: float
        Indices for the limits in the spectrumerror_star array
    model: string
        deafault "None"; 

        
    Return 
    ----------
    1-d array
        brightness of flare on top of Earth atmosphere
    1-d array
        error on brightness of flare on top of Earth atmosphere
    '''
    if model == "blackbody":
        #sum model
        brightness_star = np.asarray([np.nansum(brightness_othick_(wavelength_TRAP*1e-10, 0,0, T_star)[limit[i][0]:limit[i][1]]) for i in range(len(dfs)) ])
        T_e = T_error + T_star
        brightnesserror_star = np.asarray([np.nansum(brightness_othick_(wavelength_TRAP*1e-10, 0,0, T_e)[limit[i][0]:limit[i][1]]) for i in range(len(dfs)) ])
        brightnesserror_star = brightnesserror_star - brightness_star
    else:
        #brighntess out of spectrum
        brightness_star = np.asarray([np.nansum(flux_TRAP[limit[i][0]:limit[i][1]]) for i in range(len(dfs)) ] ) * 1e7 #erg Angstrom cm^-2  conversion constant 
        brightnesserror_star = np.asarray([np.nansum(flux_TRAP_error[limiterror[i][0]:limiterror[i][1]]) for i in range(len(dfs)) ] ) * 1e7
    return brightness_star, brightnesserror_star
