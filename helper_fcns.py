import numpy as np
import scipy.optimize as opt
from scipy.stats import norm, binom, weibull_min

def pmf_loss(eval_at, n_corr, n_trials, pmf_model, pmf_args, lapse, ceiling):

    pmf_vals = lapse + (1-lapse-ceiling)*pmf_model(eval_at, *pmf_args);
    
    loss_eval = binom.pmf(n_corr, n_trials, pmf_vals);
    
    return sum(-np.log(np.maximum(1e-6, loss_eval))); # ensure no value is zero...can't take log of zero, yadig?
    
def opt_pmf(eval_at, n_corr, n_trials, pmf_model, nFits):
        
    obj = lambda params: pmf_loss(eval_at, n_corr, n_trials, pmf_model, params[0:2], params[2], params[3]) 
    
    init_params = [np.median(eval_at), 1, 0.05, 0.02]; # location at median of all values, scale = 1
    methodStr = 'L-BFGS-B';
    sfConstr = (eval_at[0], eval_at[-1]);
    slopeConstr = (0.25, 6)
    lapseConstr = (0, 0.2)
    ceilingConstr = (0, 0.2)

    boundsAll = (sfConstr, slopeConstr, lapseConstr, ceilingConstr);
    
    bestOpt = opt.minimize(obj, init_params, method=methodStr, bounds=boundsAll);
    
    for i in range(nFits-1):
        init_params = [np.median(eval_at), slopeConstr[0] + np.diff(slopeConstr)*np.random.random(), 
                       lapseConstr[0]+np.diff(lapseConstr)*np.random.random(), 
                       ceilingConstr[0]+np.diff(ceilingConstr)*np.random.random()]; # location at median of all values, scale = 1

        if np.mod(i, 2) == 1:
            methodStr = 'L-BFGS-B';
        else:
            methodStr = 'SLSQP';
        
        currOpt = opt.minimize(obj, init_params, method=methodStr, bounds=boundsAll);

        if currOpt['fun'] < bestOpt['fun']:
            bestOpt = currOpt;
    
    return bestOpt