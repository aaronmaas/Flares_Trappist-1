import numpy as np
import scipy.constants as c

def brightness_mod_continous(wavelength, T_flare, a, fluxdensity_tot, T_star, R_star, dist_star, model = "blackbody"):
    pre_factor = 1/(dist_star**2) * 2*c.h*c.c**2 /(wavelength**5) * np.pi
    flare_contribution = (a*R_star)**2/(np.exp((c.h * c.c)/(T_flare*c.k*wavelength)) -1)
    
    if model == "blackbody": 
        star_contribution = -(a*R_star)**2/(np.exp((c.h * c.c)/(T_star*c.k*wavelength)) -1) #if blackbody here be carful with R_star its not correct at the moment
    else: 
        star_contribution = -(a)**2* fluxdensity_tot  
        #fluxdensity exaclty at the points from griz otherwise it will not work
    return pre_factor * (flare_contribution)  + star_contribution 
	#das mit dem dist star kommt mir sehr verdächtig vor überlege noch mal genau, wie du das einbinden musst
	
def brightness_mod(wavelength_tot, limit, T_flare, a, fluxdensity_tot, T_star, R_star, dist_star, model = "thick"):
    
    
    brightness = []
    for i in range(len(limit)):
        wavelength = wavelength_tot[limit[i][0]:limit[i][1]]
        pre_factor = 1/(dist_star**2) * 2*c.h*c.c**2 /(wavelength**5) * np.pi
        flare_contribution = (a*R_star)**2/(np.exp((c.h * c.c)/(T_flare*c.k*wavelength)) -1)
        if model == "thick":
            #total star in given passband #das ist super klein, deswegen hat es kaum einen Einfluss
            star_contribution =  -(a*R_star/dist_star)**2 * np.nansum(fluxdensity_tot[i]) * 1e7
        else:
            #thin case, a sign
            star_contribution = (a*R_star)**2 * np.nansum(fluxdensity_tot[i]) * 1e7
        brightness.append(np.nansum(pre_factor * (flare_contribution))  + star_contribution) 
        
    return np.asarray(brightness)


def _brightness_mod(wavelength_tot, limit, T_flare, a, fluxdensity_tot, T_star, R_star, dist_star, model = "thick"):
    
    '''
    Modified brightness calculator for every passband. It summs the brightness in the given limits. 
    So for Muscat we get four summed points. The model uses an optical thick approach. 

    Parameters
    -------------
    wavelength_tot : 1-d array
        total wavelength array to be sliced into passband contributions by limit
    limit: mx2-array
        m is the number of passbands, 2 because two boarders 
    T_flare: float
        Flare temperature [K]
    a: float
        Fraction of Radius of the star asociate with flaring area. 
    fluxdensity_tot: 1-d array
        total fluxdensity of the star given [W/m^3] 
    Stellar Paramer 
    --------------
    T_star: float 
        Temperature of star [K]
    R_star: float
        Radius of star [m]
    dist_star
        distance of star [m]
        
    Return 
    ----------
    mx1 array
        Apparent brightness of flare in the passbands at Earth
    
    '''
    
    
    brightness = []
    for i in range(len(limit)):
        wavelength = wavelength_tot[limit[i][0]:limit[i][1]]
        pre_factor = 1/(dist_star**2) * 2*c.h*c.c**2 /(wavelength**5) * np.pi
        flare_contribution = (a*R_star)**2/(np.exp((c.h * c.c)/(T_flare*c.k*wavelength)) -1)
        if model == "thick":
            #total star in given passband #das ist super klein, deswegen hat es kaum einen Einfluss
            star_contribution = -(a**2) * fluxdensity_tot[limit[i][0]:limit[i][1]]   #R_star/dist_star weg #es macht Sinn nur den faktro a da zu haben, wegen l = sum flux l_a = a sum flux
        else: 
            #blackbody optical thick atM
            star_contribution = -(a**2*R_star**2)/(np.exp((c.h * c.c)/(T_star*c.k*wavelength)) -1) 
        brightness.append(np.trapz(pre_factor * (flare_contribution) + star_contribution) )
        
    return np.asarray(brightness)
