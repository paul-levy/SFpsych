import numpy as np
import scipy.optimize as opt
from scipy.stats import norm, binom, weibull_min

import pdb

def pmf_loss(eval_at, n_testResp, n_trials, pmf_model, pmf_means, pmf_slope, lapse, ceiling):

    # n_testResp and n_trials are nCons x nSF, so we loop over all contrasts separately (each gets its own mean param)
    loss_by_con = np.zeros((n_testResp.shape[0], 1));
    for c in range(n_testResp.shape[0]):
        if n_testResp.shape[0] == 1: # only one pmf_means value
            pmf_args = [pmf_means, pmf_slope]; 
        else: # it's an array - grab the right one
            pmf_args = [pmf_means[c], pmf_slope];
        pmf_vals = lapse + (1-lapse-ceiling)*pmf_model(eval_at, *pmf_args);
        loss_eval = binom.pmf(n_testResp[c], n_trials[c], pmf_vals);
        loss_by_con[c] = sum(-np.log(np.maximum(1e-6, loss_eval))); # ensure no value is zero...can't take log of zero, yadig?
    
    return sum(loss_by_con);
    
def opt_pmf(eval_at, n_testResp, n_trials, pmf_model, nFits):
        
        # params[0-nCon] will be mean (if CDF)
        # params[nCon]   will be slope (sigma)
        # params[nCon+1,nCon+2] will be lapse (floor/ceiling)
    nCon = n_testResp.shape[0]; # number of contrast conditions
    obj = lambda params: pmf_loss(eval_at, n_testResp, n_trials, pmf_model, params[0:4], params[4], params[5], params[6]) 
    
    distr_arg1 = np.repeat(np.median(eval_at), nCon);
    init_params = np.hstack((distr_arg1, [1, 0.05, 0.02])); # scale = 1, floor = 0.05, ceil = 0.02
    methodStr = 'L-BFGS-B';
    z = [eval_at[0], eval_at[-1]]; # ugly python
    sfConstr = [tuple(x) for x in np.broadcast_to(z, (nCon, 2))]
    slopeConstr = (0.25, 6)
    lapseConstr = (0, 0.2)
    ceilingConstr = (0, 0.2)
    
    boundsAll = np.vstack((sfConstr, [slopeConstr, lapseConstr, ceilingConstr]));
    boundsAll = [tuple(x) for x in boundsAll]; # turn the (innter) arrays into tuples...
    
    bestOpt = opt.minimize(obj, init_params.tolist(), method=methodStr, bounds=tuple(boundsAll)); # turn the array into a tuple
    
    # doing more fits?
    for i in range(nFits-1):
        distr_arg1 = np.repeat(np.median(eval_at), nCon);
        init_slope = slopeConstr[0] + np.diff(slopeConstr)*np.random.random();
        init_floor = lapseConstr[0]+ np.diff(lapseConstr)*np.random.random();
        init_ceil  = ceilingConstr[0]+ np.diff(ceilingConstr)*np.random.random(); # location at median of all values, scale = 1

        #pdb.set_trace();
        
        init_params = np.hstack((distr_arg1, np.hstack((init_slope, init_floor, init_ceil))));
                           
        if np.mod(i, 2) == 1:
            methodStr = 'L-BFGS-B';
        else:
            methodStr = 'SLSQP';
        
        currOpt = opt.minimize(obj, init_params.tolist(), method=methodStr, bounds=boundsAll);

        if currOpt['fun'] < bestOpt['fun']:
            bestOpt = currOpt;
    
    return bestOpt