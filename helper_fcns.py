import numpy as np
import scipy.optimize as opt
from scipy.stats import norm, binom, weibull_min
import statsmodels.api as sm
from scipy.special import logit

import pdb

def geo_mean(nparr):
    
    n = len(nparr);
    return np.power(np.prod(nparr), 1.0/n);

def find_pse(true_pmf, params):
    # 0.5 is PSE; i.e., we believe the test and reference SF are equivalent
    pse_obj = lambda x: np.square(0.5 - true_pmf(params, x));
    pse_opt = opt.minimize(pse_obj, 0.5);
    return pse_opt['x'], pse_opt;

def fit_pmf_logit(obs, sfs):
    '''
    obs: binary vector of responses (test SF higher [1] or lower [0] than reference SF?)

    sfs: spatial frequency of the test grating; corresponds to obs
    '''
    regressors = sm.add_constant(sfs);    
    glm_binom = sm.GLM(obs, regressors, family=sm.families.Binomial(sm.families.links.logit));
    return glm_binom.fit();

def get_logit_results(glm_fit, eval_sfs = np.logspace(np.log10(1), np.log10(10), 31), pse = 0.5):
    '''
    Pass in a glm.fit() structure
    
    This will return the PSE, standard error of the estimate, a curve (given eval_sfs), etc...
    '''
    glm_params = glm_fit.params;
    stand_errors = glm_fit.bse;
    
    pse = (logit(pse)-glm_params[0])/glm_params[1];
    se_pse = abs(stand_errors[0]/glm_params[1]);
    
    predict_this = sm.add_constant(eval_sfs)
    curve = glm_fit.predict(predict_this);
    
    return pse, se_pse, curve;
    
def pmf_loss(eval_at, n_testResp, n_trials, pmf_model, pmf_means, pmf_slope, lapse_low, lapse_high):

    # n_testResp and n_trials are nCons x nSF, so we loop over all contrasts separately (each gets its own mean param)
    loss_by_con = np.zeros((n_testResp.shape[0], 1));
    nCons = n_testResp.shape[0];
    for c in range(nCons):
        if n_testResp.shape[0] == 1: # only one pmf_means value
            pmf_args = [pmf_means, pmf_slope, lapse_low, lapse_high];
        else: # it's an array - grab the right one
            pmf_args = [pmf_means[c], pmf_slope[c], lapse_low[c], lapse_high[c]];
        # Fit to log2(spatial_frequency)
        pmf_vals = pmf_args[3] + (1-pmf_args[2]-pmf_args[3])*pmf_model(np.log2(eval_at), *pmf_args[0:2]);
        loss_eval = binom.pmf(n_testResp[c], n_trials[c], pmf_vals);
        loss_by_con[c] = sum(-np.log(np.maximum(1e-6, loss_eval))); # ensure no value is zero...can't take log of zero, yadig?
    
    return sum(loss_by_con);
    
def opt_pmf(eval_at, n_testResp, n_trials, pmf_model, nFits):
        
        # params[0:nCon] will be mean (if CDF)
        # params[nCon:2*ncon]   will be slope (sigma)
        # params[2*nCon:3*ncon] will be lapse (floor)
        # params[3nCon:4*con] will be lapse (ceiling)
    nCon = n_testResp.shape[0]; # number of contrast conditions
    obj = lambda params: pmf_loss(eval_at, n_testResp, n_trials, pmf_model, params[0:nCon], params[nCon:2*nCon], params[2*nCon:3*nCon], params[3*nCon:]) 
    
    distr_init = np.repeat(np.median(eval_at), nCon);
    scale_init = np.repeat(1, nCon);
    lapse_init = np.repeat(0.05, nCon);
    init_params = np.hstack((distr_init, scale_init, lapse_init, lapse_init)); # scale = 1, lapses (both) = 0.05
    # init_params = np.hstack((distr_init, scale_init, [0.05, 0.05])); # scale = 1, lapses (both) = 0.05
    methodStr = 'L-BFGS-B';
    z = [eval_at[0], eval_at[-1]]; # ugly python
    sfConstr = [tuple(x) for x in np.broadcast_to(z, (nCon, 2))]
    slopeConstr = [tuple(x) for x in np.broadcast_to([0.25, 6], (nCon, 2))];
    lapseConstr = [tuple(x) for x in np.broadcast_to([0, 0.1], (nCon, 2))];
    # lapseConstr = (0, 0.2)
    
    boundsAll = np.vstack((sfConstr, slopeConstr, lapseConstr, lapseConstr));
    # boundsAll = np.vstack((sfConstr, slopeConstr, [lapseConstr, lapseConstr]));
    boundsAll = [tuple(x) for x in boundsAll]; # turn the (inner) arrays into tuples...
    
    bestOpt = opt.minimize(obj, init_params.tolist(), method=methodStr, bounds=tuple(boundsAll)); # turn the array into a tuple
    
    # doing more fits?
    for i in range(nFits-1):
        distr_init = np.repeat(np.median(eval_at), nCon);
        init_slope = np.repeat(slopeConstr[0][0] + np.diff(slopeConstr[0])*np.random.random(), nCon);
        init_lapse = np.repeat(lapseConstr[0][0]+ np.diff(lapseConstr[0])*np.random.random(), nCon);

        #pdb.set_trace();
        
        init_params = np.hstack((distr_init, init_slope, init_lapse, init_lapse));
        # init_params = np.hstack((distr_arg1, init_slope, np.hstack((init_lapse, init_lapse))));
                           
        if np.mod(i, 2) == 1:
            methodStr = 'L-BFGS-B';
        else:
            methodStr = 'SLSQP';
        
        currOpt = opt.minimize(obj, init_params.tolist(), method=methodStr, bounds=boundsAll);

        if currOpt['fun'] < bestOpt['fun']:
            bestOpt = currOpt;
    
    return bestOpt