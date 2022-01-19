
import numpy as np
from funcs.blackbody_model import brightness_mod

def log_likelihood(theta, x, y, yerr, limit, fluxdensity_star, T_star, R_star, dist_star):
    T, a, = theta
    model = brightness_mod(x, limit, T, a , fluxdensity_star, T_star, R_star, dist_star) 
    return -0.5 * (np.sum(((y - model)/yerr) ** 2 ))
                                               
def log_prior_global_uniform(theta):
    T, a  = theta
    if 3000 < T < 15000 and 0.0 < a < 1.:
        return 0.0
    return -np.inf
    
def log_prior_peak_uniform(theta):
    T, a  = theta
    if 3000 < T < 18000 and 0.0 < a < 1.:
        return 0.0
    return -np.inf

def log_probability(theta, x, y, yerr, limit, fluxdensity_star, T_star,R_star, dist_star, prior):
    lp = prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr, limit, fluxdensity_star, T_star, R_star, dist_star)
