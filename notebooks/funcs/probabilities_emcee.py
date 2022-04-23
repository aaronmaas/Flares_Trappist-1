
import numpy as np
import matplotlib.pyplot as plt
from funcs.blackbody_model import _brightness_mod, brightness_mod_continous
from IPython.display import display, Math
import corner


def log_likelihood(theta, x, y, yerr, limit, fluxdensity_star, T_star, R_star, dist_star):
    T, a, = theta
    model = _brightness_mod(x, limit, T, a , fluxdensity_star, T_star, R_star, dist_star) 
    return -0.5 * (np.sum(((y - model)/yerr) ** 2 ))
                                               
def log_loglikelihood(theta, x, y, yerr, limit, fluxdensity_star, T_star, R_star, dist_star):
    logT, loga = theta
    model = _brightness_mod(x, limit, 10**(logT), 10**(loga) , fluxdensity_star, T_star, R_star, dist_star) 
    return -0.5 * (np.sum(((y - model)/yerr) ** 2 ))


def log_loglikelihood_simoul(theta, x, y, yerr, x1, y1, yerr1, limit, fluxdensity_star, T_star, R_star, dist_star):
    logT, loga, loga1 = theta
    model = _brightness_mod(x, limit, 10**(logT), 10**(loga) , fluxdensity_star, T_star, R_star, dist_star)
    model1 = _brightness_mod(x1, limit, 10**(logT), 10**(loga1) , fluxdensity_star, T_star, R_star, dist_star)
    return -0.5 * (np.sum(((y - model)/yerr) ** 2 )) + -0.5 * (np.sum(((y1 - model1)/yerr1) ** 2 ))


#combined likelihood for two models 
                                                                                              
                                               
def log_prior_global_uniform(theta):
    T, a  = theta
    if 2000 < T < 15000 and 0.0 < a < 1.:
        return 0.0
    return -np.inf
    
def log_logprior_global_uniform(theta):
    logT, loga  = theta
    if np.log10(2000) < logT < np.log10(25000) and -10 < loga < 0:
        return 0.0
    return -np.inf

def log_logprior_peak_uniform(theta):
    logT, loga  = theta
    if np.log10(2000) < logT < np.log10(45000) and -20 < loga < 0:
        return 0.0
    return -np.inf

def log_logprior_global_uniform_simoul(theta):
    logT, loga, loga1  = theta
    if np.log10(2000) < logT < np.log10(25000) and -10 < loga < 0 and -10 < loga1 < 0:
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

def log_logprobability(theta, x, y, yerr, limit, fluxdensity_star, T_star,R_star, dist_star, logprior):
    lp = logprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_loglikelihood(theta, x, y, yerr, limit, fluxdensity_star, T_star, R_star, dist_star)

def log_logprobability_simoul(theta, x, y, yerr, x1, y1, yerr1, limit, fluxdensity_star, T_star,R_star, dist_star, logprior_simoul):
    lp = logprior_simoul(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_loglikelihood_simoul(theta, x, y, yerr, x1, y1, yerr1, limit, fluxdensity_star, T_star, R_star, dist_star)
    

def plot_walker_emcee(samples,labels = ["T", "a"], log = True): 
    if log == "True":
        labels = ["logT", "loga"]
    #Analysis Plot Walker
    for j in range(len(samples)):
        fig, axes = plt.subplots(2, figsize=(10, 7), sharex=True)
        for i in range(len(labels)):
            
            ax = axes[i]
            ax.plot(samples[j][:, :, i], "k", alpha=0.3)
            ax.set_xlim(0, len(samples[j]))
            #ax.set_ylabel(labels[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)
            axes[-1].set_xlabel("step number");
    return 
    
def plot_walker_emcee_simoul(samples,labels = ["T", "a", "a1"]): 
    #Analysis Plot Walker
    fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
    for i in range(len(labels)):
            
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        #ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)
        axes[-1].set_xlabel("step number");
    return 


def plot_corner_emcee(samples_flat, bins = 100, labels = ["T","a"], log = True):
    if log == "True":
        labels = ["logT","loga"]
    for j in range(len(samples_flat)):
        fig = corner.corner(
        samples_flat[j], labels=labels, bins = bins
        );
    return 
    
def plot_corner_emcee_simoul(samples_flat, bins = 100, labels = ["logT","loga", "loga1"]):
    
    fig = corner.corner(samples_flat, labels=labels, bins = bins);
    return 


def display_median_from_chain(samplesflat, labels = ["T", "a"]):
    
    T, Terror, a, aerror = [], [], [], []
    
    for j in range(len(samplesflat)):
        for i in range(2):
            mcmc = np.percentile(samplesflat[j][:, i], [16, 50, 84])
            q = np.diff(mcmc)
            txt = "\mathrm{{{3}}} = {0:.5f}_{{-{1:.5f}}}^{{{2:.5f}}}"
            txt = txt.format(mcmc[1], q[0], q[1], labels[i])
            display(Math(txt))
            if i == 0:
                T_med = mcmc[1]
                T.append(T_med)
                Terror.append(q)
            if i == 1:
                a_med = mcmc[1]
                a.append(a_med) 
                aerror.append(q) 
                
    return T, Terror, a, aerror  


#This is a diagnostik plot, so it does not have to be perfect
def plot_fit_with_MCMCsamples(samplesflat_total, brightness_flares, brightnesserror_flares, T_Med, a_Med, \
                              T_star, R_star, Dist_star, Flux_SED, Wavelength_M2, Wavelength_SED, Limit, \
                              Limit_sum, continous = True):
    fig, axes = plt.subplots(2, figsize=(10, 7), sharex=True)
    plt.xlabel("Wavelength [m]")
    plt.xlim(0.4e-6,1e-6)
    #plt.ylim(0)
    plt.ylabel(r"Brightness [$W m^{-2}$]")
    for i in range(len(brightness_flares)):
        ax = axes[i]
        inds = np.random.randint(len(samplesflat_total[i]), size=100) #
        if continous != True:
            for ind in inds:
                sample = samplesflat_total[i][ind]
                ax.plot(Wavelength_M2, _brightness_mod(Wavelength_SED,Limit,10**sample[0],10**sample[1], \
                                            Flux_SED,T_star,R_star, Dist_star), alpha=0.1)
            ax.errorbar(Wavelength_M2, brightness_flares[i], yerr=brightnesserror_flares[i], fmt=".k", capsize=0)
            ax.plot(Wavelength_M2, _brightness_mod(Wavelength_SED,Limit,T_Med[i],a_Med[i], Flux_SED,T_star, \
                                            R_star, Dist_star), "green", alpha=1)
        else:
            ax.errorbar(Wavelength_M2, brightness_flares[i]/Limit_sum, yerr=brightnesserror_flares[i]/Limit_sum, fmt=".k", capsize=0)
            ax.plot(Wavelength_SED, brightness_mod_continous(Wavelength_SED,T_Med[i], a_Med[i], Flux_SED, T_star, \
                                                             R_star, Dist_star, model = None))
