
# coding: utf-8

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pltSave
import matplotlib.ticker as ticker
import seaborn as sns
import helper_fcns as hlp
import autoreload
import os
import sys
from scipy.stats import norm, wilcoxon
import statsmodels.api as sm
matplotlib.style.use('plt_style.mplstyle')

import pdb;
import time;

og_time = time.time();

subj = int(sys.argv[1]);
disp = int(sys.argv[2]);
isGeorgeson = int(sys.argv[3]);
bootIter = int(sys.argv[4]);
pmfType = int(sys.argv[5]); # 1 for weibull, 2 for normcdf
if len(sys.argv) > 6:
  sfFlag = str(sys.argv[6]); # special modifier for spatial frequency?
else:
  sfFlag = '';

dataDir = 'data/';
saveDir = dataDir + 'results/';
savePlts = 1;
savePlts_logit = 0;

if pmfType == 1:
  pmf_str = 'weibull';
  weibull = lambda x, *params: 1 - np.exp(-np.power(x/params[0], params[1]));
  true_pmf = lambda params, x: params[3] + (1-params[2]-params[3])*weibull(x, *params[0:2])
elif pmfType == 2:
  pmf_str = 'norm';
  true_pmf = lambda params, x: params[3] + (1-params[2]-params[3])*norm.cdf(x, *params[0:2])

if isGeorgeson:
    whichFiles = 'toAnalyzeS%dD%d%sg.txt' % (subj, disp, sfFlag);
else:
    whichFiles = 'toAnalyzeS%dD%d%s.txt' % (subj, disp, sfFlag);
dataList = open(dataDir + whichFiles, 'r')
dFilesAll = dataList.readlines(-1)
dFilesAll = [x.strip('\n') for x in dFilesAll]
dataList.close()

if not isGeorgeson and subj == 1:
    if disp == 1:
        dFiles = dFilesAll[4:]; # ignore the first 4 sessions - early training
    elif disp == 3:
        dFiles = dFilesAll[1:];
    elif disp == 5:
        dFiles = dFilesAll[2:]; # ignore first 2
else:
    dFiles = dFilesAll;

sfIdx = 0;
conIdx = 1;
dispIdx = 2;
nInd = 3; # 3 indices per stimulus; add 3 to all of the above indices to get the equivalent value for the second stimulus
subjRespIdx = 2*nInd;
refIdx = subjRespIdx + 1; # which interval contained reference is after the subject's response in the data format

stim1 = 39; stim2 = 41; # subject responses (asked if first/left stimulus is of higher frequency than the second/right stimulus)

data = [];
for df in dFiles:
    currData = np.loadtxt(dataDir + df);
    if df == dFiles[0]:
        data = currData;
    else:
        data = np.concatenate((data, currData), axis = 0);

sfVals = np.union1d(np.unique(data[:, sfIdx]), np.unique(data[:, sfIdx+nInd]))
conVals = np.union1d(np.unique(data[:, conIdx]), np.unique(data[:, conIdx+nInd]))
nCons = len(conVals);
dispVals = np.union1d(np.unique(data[:, dispIdx]), np.unique(data[:, dispIdx+nInd]))

allResp = data[:, subjRespIdx];
if any((allResp != stim1) & (allResp != stim2)):
    print('Problem: response which is neither 1 or 2')

s1g = data[:, sfIdx] > data[:, sfIdx+nInd];
s2g = data[:, sfIdx] < data[:, sfIdx+nInd]

# find which one is the reference (as opposed to test) grating
ref = np.zeros(data.shape[0])
ref = data[:, refIdx]

# allows for possibility of more than one refSf...
refSF = np.union1d(np.unique(data[ref==1, sfIdx]), np.unique(data[ref==2, sfIdx+nInd]))

# did the subject perceive the test as higher SF?
testHF = np.zeros(data.shape[0]);
testHF[ref==1] = data[ref==1, subjRespIdx] == stim2;
testHF[ref==2] = data[ref==2, subjRespIdx] == stim1;

testCons = np.zeros(data.shape[0])
testCons[ref==1] = data[ref==1, nInd+conIdx];
testCons[ref==2] = data[ref==2, conIdx];

testSfs = np.zeros(data.shape[0])
testSfs[ref==1] = data[ref==1, nInd+sfIdx];
testSfs[ref==2] = data[ref==2, sfIdx];

checkz = np.zeros((len(testSfs), 1));
for i in range(len(sfVals)):
    checkz[testSfs == sfVals[i]] = 1;


# data analysis - number of trials/responses "test>ref sf", fraction "test>ref" split by sf/con
glmFits = [];
ptSf, stdSf, nTr, nTestResp = (np.zeros((nCons, len(sfVals))) for _ in range(4));
for con in range(nCons):
    for sf in range(len(sfVals)):
        z = (testSfs == sfVals[sf]) & (testCons == conVals[con]) # get trials with the desired test sf/con
        nTr[con][sf] = sum(z); # how many trials in this configuration
        nTestResp[con][sf] = sum(testHF[z]) # how many trials in this configuration with response "test > ref sf"
        p_curr = nTestResp[con][sf] / nTr[con][sf];
        ptSf[con][sf] = p_curr;
        stdSf[con][sf] = np.sqrt(p_curr*(1-p_curr)/nTr[con][sf]); # close form for std of binomially distributed variable

# optimize and bootstrap
nFits = 5; # don't think you need multistart, but build it in anyway
univ_params = 3; # number of "universal parameters", i.e. for all pmf
nParams = nCons + univ_params; # one mean for each contrast; overall slope, lapses
bootTestResps = np.zeros((nCons, len(nTestResp[0]), bootIter))
loss = np.zeros((nCons, 1)); # only one loss value...

if pmfType == 1:
  opt = hlp.opt_pmf(sfVals, nTestResp, nTr, weibull, nFits);
elif pmfType == 2:
  opt = hlp.opt_pmf(sfVals, nTestResp, nTr, norm.cdf, nFits);

params = np.zeros((nCons, univ_params+1)) # 3+1 params per PMF
for c in range(nCons): # now unpack...
    params[c, 0] = opt['x'][c];
    params[c, 1] = opt['x'][nCons+c];
    params[c, 2] = opt['x'][2*nCons+c];
    params[c, 3] = opt['x'][3*nCons+c];
    # expand dimensions of nTestResp so that pmf_loss knows there is only one contrast value to explore
    if pmfType == 1:
      loss[c] = hlp.pmf_loss(sfVals, np.expand_dims(nTestResp[c], 0), np.expand_dims(nTr[c], 0),
                           weibull, params[c, 0], params[c, 1], params[c, 2], params[c, 3]);
    elif pmfType == 2:
      loss[c] = hlp.pmf_loss(sfVals, np.expand_dims(nTestResp[c], 0), np.expand_dims(nTr[c], 0),
                           norm.cdf, params[c, 0], params[c, 1], params[c, 2], params[c, 3]);
    
paramsBoot = np.zeros((nCons, univ_params+1, bootIter));
for iter in range(bootIter):
    bootTestResps[:, :, iter] = np.random.binomial(nTr.astype('int64'), ptSf); # i.e. parameteric bootstrap - draw nTr samples from binomial distribution with probability ptSf
    # expand dimensions of bootTestResps so that pmf_loss knows there is only one contrast value to explore
    if pmfType == 1:
      optBoot = hlp.opt_pmf(sfVals, bootTestResps[:, :, iter], nTr, weibull, 1)
    elif pmfType == 2:
      optBoot = hlp.opt_pmf(sfVals, bootTestResps[:, :, iter], nTr, norm.cdf, 1)
    # unpack paramters
    for c in range(nCons):
        paramsBoot[c, 0, iter] = optBoot['x'][c];
        paramsBoot[c, 1, iter] = optBoot['x'][nCons+c];
        paramsBoot[c, 2, iter] = optBoot['x'][2*nCons+c];
        paramsBoot[c, 3, iter] = optBoot['x'][3*nCons+c];

# compute from bootstrap values...
bootParamMean = paramsBoot.mean(-1);
bootParamStd = paramsBoot.std(-1);
nParams = paramsBoot.shape[1];
bootParamErrBars = np.zeros((2, nCons, nParams)); # plt.errorbar wants (2, N)
for i in range(nParams):
  bootParamErrBars[0, :, i] = [np.percentile(x, 2.5) for x in paramsBoot[:, i, :]]; 
  bootParamErrBars[1, :, i] = [np.percentile(x, 97.5) for x in paramsBoot[:, i, :]];

# reshape number of trials by con/sf so we can get bootstrap estimates of subject's proportion responses "test>ref"
bootPtSf = bootTestResps / np.repeat(nTr.reshape((len(conVals), len(nTestResp[0]), 1)), bootIter, axis=-1); 
bootPtSfmean = bootPtSf.mean(-1);
bootPtSfstd = bootPtSf.std(-1);

# ### Plot all PMF
isPmfLog = 0;

from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica']
rcParams['font.style'] = 'oblique'
rcParams['font.size'] = 15
rcParams['pdf.fonttype'] = 3 # should be 42, but there are kerning issues
rcParams['ps.fonttype'] = 3 # should be 42, but there are kerning issues

sns.set_style('white')

figPMF, allPMF = plt.subplots(1, nCons, figsize=(5*nCons, 4));
sns.despine(offset=10);

if isPmfLog:
    sfValsPlot = np.log2(sfVals);
    refSfPlot = np.log2(refSF);
    pmfPlot = np.log2(np.linspace(sfVals[0], sfVals[-1], 101));
else:
    sfValsPlot = sfVals;
    refSfPlot = refSF;
    pmfPlot = np.linspace(sfVals[0], sfVals[-1], 101);

for c in range(nCons):
    whichCon = c;

    pts = allPMF[c].errorbar(sfValsPlot, ptSf[whichCon], yerr=stdSf[whichCon], clip_on=False, linestyle='none', marker='o');
    refSFline = allPMF[c].plot([refSfPlot, refSfPlot], [0, 1], '--'); 
    pse = allPMF[c].plot([sfValsPlot[0], sfValsPlot[-1]], [0.5, 0.5], '--')
    pmf = allPMF[c].plot(pmfPlot, true_pmf(params[whichCon, :], pmfPlot), '-');
    
    if not isPmfLog:
        allPMF[c].set_xscale('log')
    
    allPMF[c].set_ylim([0, 1]);
    
    # get the ticks correct
    allPMF[c].tick_params(axis='x', which='minor', bottom=False, labelbottom=False); # remove old/original ticks
    allPMF[c].tick_params(labelsize=15, width=2, length=15, which='major')
    sfAsStr = [str(x) for x in np.round(sfValsPlot, 2)];
    allPMF[c].set_xticks(sfValsPlot);
    allPMF[c].set_xticklabels(sfAsStr, rotation='45');
    
    allPMF[c].set_title('Contrast: ' + str(conVals[whichCon]) + ' -- ' + str(sum(nTr[whichCon])) + ' trials');
    if c == 0:
        allPMF[c].legend((pts[0], refSFline[0], pse[0], pmf[0]), ('subject responses', 'reference SF', 'PSE', 'fit pmf'), loc='upper left', fontsize=12);
        allPMF[c].set_xlabel('test SF (c/deg)');
        allPMF[c].set_ylabel('prop "test" higher SF');

figPMF.tight_layout()

# What are the inferred PSEs for each contrast? How about the SF which is perceived at the reference SF?

fit_pse = np.zeros((nCons, 1));
opts = [];
for c in range(nCons):
    fit_pse[c, 0], opts_curr = hlp.find_pse(true_pmf, params[c, :], [sfVals[0], sfVals[-1]])
    opts.append(opts_curr);

# ### Plot of subject bias and sensitivity

pse = np.zeros((len(conVals), 1));
evalSfs = np.arange(sfVals[0], sfVals[-1], 1e-3);
evalPmf = lambda params: true_pmf(params, evalSfs);
#evalPmf = lambda params: params[2] + (1-2*params[2])*norm.cdf(evalSfs, *params[0:2]);
for con in range(len(conVals)):
    pse[con] = evalSfs[np.argmin(abs(evalPmf(params[con]) - 0.5))]

# Statistics on PSE bootstrapped values...
pOfInt = 0; # PSE
outerInd = len(conVals);
pValPSE = np.nan*np.ones((outerInd, outerInd));
for c in range(outerInd):
    for ci in range(len(conVals) - c-1):
        
        z, pValPSE[c, c+ci+1] = wilcoxon(paramsBoot[c, pOfInt, :], paramsBoot[c+ci+1, pOfInt, :]); # 0 is PSE
        pValPSE[c+ci+1, c] = pValPSE[c, c+ci+1]


nParams = params.shape[1];
figParams, paramPlot = plt.subplots(1, nParams, figsize=(5*nParams, 4));

sns.set_style('white')
sns.despine(offset=10);

y_labels = ['SF (cpd)', 'sensitivity', 'lapse rate', 'lapse rate']
titles = ['Subject bias', 'Subject sensitivity', 'Lapse (ceiling)', 'Lapse (floor)']

for param_i in range(nParams):
    
    if param_i == 0:
        pse = paramPlot[param_i].errorbar(conVals, fit_pse[:, 0], np.abs(bootParamErrBars[:, :, param_i]-fit_pse[:, 0]));
        ref_plt = paramPlot[param_i].axhline(refSF, ls='--', c='k');
    else:
        paramPlot[param_i].errorbar(conVals, params[:, param_i], np.abs(bootParamErrBars[:, :, param_i]-params[:, param_i]));
        
    paramPlot[param_i].tick_params(axis='x', which='minor', bottom=False, labelbottom=False); # remove old/original ticks
    paramPlot[param_i].set_xscale('log');
    paramPlot[param_i].set_xlabel('test contrast');
    paramPlot[param_i].set_ylabel(y_labels[param_i]);
    paramPlot[param_i].set_title(titles[param_i]);
    paramPlot[param_i].tick_params(labelsize=15, width=2, length=15, which='major')
    paramPlot[param_i].set_xticks(conVals);
    conAsStr = [str(x) for x in conVals];
    paramPlot[param_i].set_xticklabels(conAsStr);

figParams.tight_layout()

if savePlts:
    allFigs = [figPMF, figParams];
    if isGeorgeson:
        saveName = "results_%sS%dD%d%sg.pdf" % (pmf_str, subj, disp, sfFlag)
    else:
        saveName = "results_%sS%dD%d%s.pdf" % (pmf_str, subj, disp, sfFlag)
    pdf = pltSave.PdfPages(str(saveDir + saveName))
    for fig in range(len(allFigs)):
        pdf.savefig(allFigs[fig], bbox_inches="tight")
    pdf.close()

print('took %d seconds' % (time.time() - og_time));
