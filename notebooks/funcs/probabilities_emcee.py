
import numpy as np
import matplotlib.pyplot as plt
from funcs.blackbody_model import _brightness_mod
from IPython.display import display, Math
import corner


def log_likelihood(theta, x, y, yerr, limit, fluxdensity_star, T_star, R_star, dist_star):
    T, a, = theta
    model = _brightness_mod(x, limit, T, a , fluxdensity_star, T_star, R_star, dist_star) 
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
    

def plot_walker_emcee(samples,labels = ["T", "a"]): 
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

def plot_corner_emcee(samples_flat):
    labels = ["T","a"]
    for j in range(len(samples_flat)):
        fig = corner.corner(
        samples_flat[j], labels=labels
        );
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
