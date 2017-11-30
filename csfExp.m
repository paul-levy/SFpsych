% This script will call csfCore to run the core contrast sensitivity
% function experiment; here, we choose the parameters of that experiment
% (i.e. what spatial frequencies, how many trials, etc)

clear;

% useful functions
fnShuffle = @(x) x(randperm(length(x))); % from LH
myRound = @(x, digit) round((x.*10^digit))./10^digit;
sf_round = 2; % just round to X digits...

%% experimental set up
num_trials = 100; % per CSF call
subj_id = 's0l'; % as in, subject01
data_dir = 'data/CSF/';
view_dist_m = 1.14; % 114 cm

% sf information
sf_cent = 2; % in cpd
surr_oct = 3; % go out +/- X octaves
num_steps = 10;

% NOTE: TEST_CON_PERCENT are the current contrast values used in csfCore
TEST_CON_PERCENT = power(2, [-3 : 0.5 : 6.5]);
start_step = randi(length(TEST_CON_PERCENT), [num_steps, 1]);

%% spatial frequency calculation - 
octSeries  = linspace(-surr_oct, surr_oct, num_steps);
freqSeries = myRound(2.^(octSeries(:) + log2(sf_cent)), sf_round);
freqOrder = fnShuffle(freqSeries);

%% Now run the experiments
for f = 1 : length(freqOrder)
    
    csfCore(data_dir, subj_id, view_dist_m, freqOrder(f), start_step(f), num_trials);
    
    waitforbuttonpress();
end



