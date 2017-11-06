import numpy as np
import scipy.optimize as opt
from scipy.stats import norm, binom, weibull_min

def pmf_loss(eval_at, n_corr, n_trials, pmf_model, pmf_args):
    
    pmf_vals = pmf_model(eval_at, *pmf_args);
    
    loss_eval = binom.pmf(n_corr, n_trials, pmf_vals);
    
    return sum(-np.log(loss_eval));
    
def opt_pmf(eval_at, n_corr, n_trials, pmf_model):
        
    obj = lambda params: pmf_loss(eval_at, n_corr, n_trials, pmf_model, params) 
    
    init_params = [np.median(eval_at), 1]; # location at median of all values, scale = 1
    methodStr = 'L-BFGS-B';
    
    rapper = opt.minimize(obj, init_params, method=methodStr);
    
    return rapper