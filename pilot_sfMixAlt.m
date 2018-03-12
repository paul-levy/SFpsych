%%%%%%%%%%
% SF Mixture piloting program. PL, 20180309.
% Initial code from LH
%%%%%%%%%%

clear;

% useful functions...
myRound = @(x, digit) round((x.*10^digit))./10^digit;
sf_round = 3; % round SF values to the thousandth
scaleGrat = @(sf, con) (255/2)*(1+sf.*con);

%% Experiment parameters
% saving the file
subj = 99;
save_loc = 'data/'; 
meta_loc = 'data/metaData/';
is_pilot = 1;
run_num = 400;
save_meta = 0; % save metadata?

NUM_TRIALS = 150;
REF_CON = 1; % do not change from 1 - 03.12.18
REF_DISP = 4; % 1-4 are our options, as of now (03.09.18)
SF_REF = 1.73; % in cpd - center of spatial frequency vector
incMidSamp = 1; % sample more near (i.e. -1/0/+1 rel. to) the reference?

stimDist = 6; % i.e. 6 degrees in periphy
xsgn = 1; ysgn = -1; % i.e. pos or neg [x/y]

debug_t = 0; % print how long each frame takes?
feedback = 0; % give feedback to the subject?

com_str = 'Pilot in VNL psychophysics room (1139) in dark conditions.';
% comments for meta data file (e.g. room, condition, etc)
%% Design of experiment - (blank?:)stim1:blank:stim2
inter_1 = 0.2; % how many seconds for stimulus 1?
inter_2 = 0.2; % stimulus 2?
inter_blank = 1; % intervening blank?
min_iti = 1; % the minimum inter-trial-interval is 1 second; wait after response if needed

% stimulus information
grat_step = 0.5059; % spacing between adjacent SFs in mixture stimuli in octaves; 0.5059 octaves is spacing in logspace(log10(0.3), log10(10), 11)
cent_step = 0.25; % i.e. space the center of each distribution X octaves apart
n_cent_steps = 11; % must be odd s.t. we can symmetrically go about SF_REF
if iseven(n_cent_steps)
    n_cent_steps = n_cent_steps+1;
end

tf = 5;
tf_spread = tf/5; % draw TF of dispersed gratings from gaussian with this sigma

if REF_DISP == 1
    testDisps = 1; % [1 5];
else
    testDisps = REF_DISP;
end

ori_all = 90; % vertical grating (drift in horizontal)

% stim. location and size
stim_radius = 1; % radius, in degrees
stim_loc = [xsgn*sqrt((stimDist^2)/2), ysgn*sqrt((stimDist^2)/2)];
fp_radius = 0.1; % in degrees
col_fix = [1 1 1]; % for fixation point
slack = 2; % make the grating X times the size of the stencil/aperture...
%% display set up
mon = 2; % in VNL psych room, 0 is "work" monitor, 2 is actual monitor

% Open, init MGL.
if mon == 2
    xPix = 1280; yPix = 960;
    mglOpen(mon, xPix, yPix, 100, 32); %100Hz, bit depth = 32
    scrXcm = 40; scrYcm = 30; % in cm
    scrDist = 114; % in cm
    mglVisualAngleCoordinates(scrDist, [scrXcm scrYcm]);
else
    xPix = 800; yPix = 800;
    mglOpen(mon, xPix, yPix, 60, 32); % 800x800, 60Hz, bit depth = 32
    scrXcm = 16; scrYcm = 10; % in cm
    scrDist = 57; % in cm
    mglVisualAngleCoordinates(scrDist,[scrXcm scrYcm]); % viewing from 57 cm, display is 5cm x 5cm
end
%% For saving results
if is_pilot
    save_base = sprintf('sfPer_s%02d_p%03d', subj, run_num);
    if save_meta
        sm_base = sprintf('META_sfPer_s%02d_p%03d', subj, run_num);
    end
else
    save_base = sprintf('sfPer_s%02d_%03d', subj, run_num);
    if save_meta
        sm_base = sprintf('META_sfPer_s%02d_p%03d', subj, run_num);
    end
end

if isfile([save_loc, save_base, '.txt']) % if it exists, append the time so we don't overwrite
    exp_info = fopen([save_loc, save_base, datestr(now), '.txt'], 'w+');
else
    exp_info = fopen([save_loc, save_base, '.txt'], 'w+');
end

if save_meta
    if isfile([meta_loc, sm_base, '.txt'])
        meta_inf = fopen([meta_loc, sm_base, datestr(now), '.txt'], 'w+');
    else
        meta_inf = fopen([meta_loc, sm_base, '.txt'], 'w+');
    end
    
    fprintf(meta_inf, 'number of trials: %g\nstimulus eccentricity in degrees: %.3f\n', NUM_TRIALS, stimDist);
    fprintf(meta_inf, 'length of interval 1 in s: %.4f\n', inter_1);
    fprintf(meta_inf, 'length of blank interval in s: %.4f\n', inter_blank);
    fprintf(meta_inf, 'length of interval 2 in s: %.4f\n', inter_2);
    fprintf(meta_inf, 'stimulus location x: %.4f\n', stim_loc(1));
    fprintf(meta_inf, 'stimulus location y: %.4f\n', stim_loc(2));
    fprintf(meta_inf, 'stimulus radius degrees: %.3f\n', stim_radius);
    fprintf(meta_inf, 'fixation point radius in degrees: %.3f\n', fp_radius);
    fprintf(meta_inf, 'screen resolution (x, y) in pixels: %g, %g\n', xPix, yPix);
    fprintf(meta_inf, 'screen size (x, y) in cm: %g, %g\n', scrXcm, scrYcm);
    fprintf(meta_inf, 'viewing distance in cm: %g\n', scrDist);
    fprintf(meta_inf, 'subject feedback? %g\n', feedback);
    fprintf(meta_inf, 'start time: %s\n', datestr(now));
    fprintf(meta_inf, 'comments: %s\n', com_str);
    fclose(meta_inf);
end
%% stimulus creation/calculation
% center spatial frequencies
freqMax = 2^(log2(SF_REF) + floor(n_cent_steps/2)*cent_step);
freqMin = 2^(log2(SF_REF) - floor(n_cent_steps/2)*cent_step);
freqCenters = logspace(log10(freqMin), log10(freqMax), n_cent_steps);

% if you want one more point flanking either side ...
% of the reference (log space between adjacent and mid)
if incMidSamp
    logMid = @(a, b) 2^((log2(a) + log2(b))/2);
    ref_ind = find(myRound(freqCenters, sf_round) == myRound(SF_REF, sf_round), 1, 'first'); % always the same; take only 1 if the value is found > 1 time
    lowMid = logMid(freqCenters(ref_ind-1), freqCenters(ref_ind));
    highMid = logMid(freqCenters(ref_ind), freqCenters(ref_ind+1));
    freqCenters = sort([freqCenters, freqCenters(ref_ind), lowMid, lowMid, highMid, highMid]);
end
freqCenters = myRound(freqCenters, sf_round);
    
% dispersion?
num_gratings = 7; % fixed from sfMixAlt physiology experiments
freqMax = 2^(log2(1) + floor(num_gratings/2)*grat_step); % For expl. on log2(1), see comment below sfVec
freqMin = 2^(log2(1) - floor(num_gratings/2)*grat_step);
sfVec = logspace(log10(freqMin), log10(freqMax), num_gratings);
 % relative to a particular sfCenter, sfVec contains the factors which can
 % be used to multiply the sfCenter to create the dispersed grating

% contrast
logSpCons = logspace(log10(0.05), log10(1), 9); % 9 possible contrasts - from sfMixAlt physiology
logSpConsRev = fliplr(logSpCons);
valCons = [9 4 3 2]; % first N of total_cons can be used for corresponding dispersion (inc. disp L to R)
conStart = [0 1 3 3]; % N steps from end of logSpCons list
% indexing into the logSpCons list is simply replication of what is used in
% sfMixAlt (as of 03.09.18 - PL)
conProfile = cell(length(valCons), length(valCons)); % sets all to empty
for dispInd = 1:length(valCons)
    for conInd = 1 : valCons(dispInd)
        currConStart = length(logSpCons) - conStart(dispInd) - (conInd-1);
        if dispInd == 1
            conProfile{dispInd, conInd} = [0 0 0 logSpCons(currConStart) 0 0 0];
        elseif dispInd == 2
            conProfile{dispInd, conInd} = [0 0 logSpCons(currConStart-4) logSpCons(currConStart) logSpCons(currConStart-4) 0 0];            
        elseif dispInd == 3
            conProfile{dispInd, conInd} = [0 logSpCons(currConStart-3)  logSpCons(currConStart-1) logSpCons(currConStart) logSpCons(currConStart-1) logSpCons(currConStart-3) 0];
        elseif dispInd == 4
            conProfile{dispInd, conInd} = [logSpCons(currConStart-4) logSpCons(currConStart-3)  logSpCons(currConStart-2) logSpCons(currConStart) logSpCons(currConStart-2) logSpCons(currConStart-3) logSpCons(currConStart-4)];    
        end
    end
end

warning('Contrast randomization only works under the assumption of one dispersion throughout the experiment');
TEST_CONS = 1:valCons(testDisps); % i.e. how many valid contrasts?indices into conProfile cell above

% Stencils
mglStencilCreateBegin(1);
mglFillOval(stim_loc(1), stim_loc(2), [2*stim_radius 2*stim_radius]);
mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], [1 1 1]); % fixation point
mglStencilCreateEnd;

% Create mean-lum texture and clear screen.
ml = 255*0.5*ones(yPix,xPix);
texml = mglCreateTexture(ml);
mglBltTexture(texml,[0 0]); mglFlush;

% determine order:
rng('shuffle');
which_ref = randi(2, NUM_TRIALS, 1);
test_ind = randi(length(freqCenters), NUM_TRIALS, 1);
test_disp = testDisps(randi(length(testDisps), NUM_TRIALS, 1));
test_con = TEST_CONS(randi(length(TEST_CONS), NUM_TRIALS, 1));
ref_ind = find(freqCenters == SF_REF, 1, 'first'); % always the same; take only one index

if REF_DISP == 1 % Create single sinusoids
    % just fix
    clear grat;
    ori = ori_all;
    for sf_c = 1 : length(freqCenters)
        ph = 360*rand(); % phase in degrees
        % slack*2*stim_radius because the stencil (i.e. stimulus area) is
        % 2*stim_radius wide/tall (diameter = 2*radius, ya!)
        grat{sf_c} = mglMakeGrating(slack*2*stim_radius, slack*2*stim_radius, freqCenters(sf_c), ori, ph);
    end
    % create only the reference texture; others will be computed as needed
    tex_ref = mglCreateTexture(scaleGrat(grat{ref_ind}, REF_CON));
end 

%%%
%% Main loop is here.
%%%
mglBltTexture(texml,[0 0]);
mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], col_fix);
mglFlush;
mglWaitSecs(5);
mglPlaySound('submarine');
mglWaitSecs(2);

n_corr = 0;

for tr_i = 1:NUM_TRIALS

  mglBltTexture(texml,[0 0]); mglFlush;
  mglBltTexture(texml,[0 0]); mglFlush;
  mglStencilSelect(1);

  if which_ref(tr_i) == 1
      sf1 = ref_ind;
      con1 = REF_CON;
      disp1 = REF_DISP;
      if disp1 == 1
          tex1 = tex_ref;
      else % need to create random phases, tf; grab con profile
          ph_grat1 = 360*rand(num_gratings, 1);
          tf_grat1 = random('norm', tf, tf_spread, [num_gratings, 1]);
          profile1 = conProfile{disp1, con1};
      end
          
      sf2 = test_ind(tr_i);
      con2 = test_con(tr_i);
      disp2 = test_disp(tr_i);
      if disp2 == 1
          tex2 = mglCreateTexture(scaleGrat(grat{sf2}, logSpConsRev(con2)));
      else % need to create random phases, tf; calculate con profile
          ph_grat2 = 360*rand(num_gratings, 1);
          tf_grat2 = random('norm', tf, tf_spread, [num_gratings, 1]);
          profile2 = conProfile{disp2, con2};
      end
  else
      sf1 = test_ind(tr_i);
      con1 = test_con(tr_i);
      disp1 = test_disp(tr_i);
      if disp1 == 1
          tex1 = mglCreateTexture(scaleGrat(grat{sf1}, logSpConsRev(con1)));
      else % need to create random phases, tf; calculate con profile
          ph_grat1 = 360*rand(num_gratings, 1);
          tf_grat1 = random('norm', tf, tf_spread, [num_gratings, 1]);
          profile1 = conProfile{disp1, con1};
      end
      
      sf2 = ref_ind;
      con2 = REF_CON;
      disp2 = REF_DISP;
      if disp2 == 1
        tex2 = tex_ref;
      else % need to create random phases, tf; calculate con profile
          ph_grat2 = 360*rand(num_gratings, 1);
          tf_grat2 = random('norm', tf, tf_spread, [num_gratings, 1]);
          profile2 = conProfile{disp2, con2};
      end
  end
  
  % first stim interval
  inter1 = clock;
  while (etime(clock,inter1) < inter_1)
 
    elapsed_time_s = etime(clock,inter1);

    if debug_t, tic; end; % for timing...
    elpsed = clock;
    if disp1 > 1 % i.e. dispersed grating
        % need to calculate
        ori = ori_all;
        sfVecCurr = sfVec * freqCenters(sf1);
        conVec = profile1; % total contrast levels built in
        for grat = 1 : num_gratings            
            curr_grat{grat} = scaleGrat(mglMakeGrating(2*stim_radius, ...
                2*stim_radius, sfVecCurr(grat), ori, ...
                mod(ph_grat1(grat) - 360*elapsed_time_s*tf_grat1(grat), 360)), conVec(grat));
            if grat == 1
                curr_stim = curr_grat{grat} - 255*0.5; % subtract off mean luminance
            else
                curr_stim = curr_stim + (curr_grat{grat} - 255*0.5);
            end
        end
        curr_stim = curr_stim + 255*0.5;
        
        tex1 = mglCreateTexture(curr_stim);
        mglBltTexture(tex1, [stim_loc(1) stim_loc(2)], 0, 0); % center
    else
        curr_x = (stim_loc(1)+stim_radius) + mod(elapsed_time_s*tf/freqCenters(sf1), (slack-1)*2*stim_radius);
        mglBltTexture(tex1, [curr_x stim_loc(2)], 1, 0); % right-align the texture
    end
    
    mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], col_fix);
    mglFlush;
    
    if debug_t, toc; end;
    
  end
  
  % blank interval
  interBlank = clock;
  while (etime(clock,interBlank) < inter_blank)
 
    elapsed_time_s = etime(clock,interBlank);
    col_fix = [1 1 1]; % for fixation point

    if debug_t, tic; end % for timing...
    elpsed = clock;
    mglBltTexture(texml,[0 0]);
    mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], col_fix);
    mglFlush;
    if debug_t, toc; end

  end
  
  % second stimulus interval
  inter2 = clock;
  while (etime(clock,inter2) < inter_2)
 
    elapsed_time_s = etime(clock,inter2);
    col_fix = [1 1 1]; % for fixation point

    if debug_t, tic; end; % for timing...
    elpsed = clock;
    if disp2 > 1 % i.e. dispersed grating
        % need to calculate
        ori = ori_all;
        sfVecCurr = sfVec * freqCenters(sf2);
        conVec = profile2; % total contrast built in
        for grat = 1 : num_gratings            
            curr_grat{grat} = scaleGrat(mglMakeGrating(2*stim_radius, ...
                2*stim_radius, sfVecCurr(grat), ori, ...
                mod(ph_grat2(grat) - 360*elapsed_time_s*tf_grat2(grat), 360)), conVec(grat));
            if grat == 1
                curr_stim = curr_grat{grat} - 255*0.5; % subtract off mean luminance
            else
                curr_stim = curr_stim + (curr_grat{grat} - 255*0.5);
            end
        end
        curr_stim = curr_stim + 255*0.5;
        
        tex2 = mglCreateTexture(curr_stim);
        mglBltTexture(tex2, [stim_loc(1) stim_loc(2)], 0, 0); % center
    else
        curr_x = (stim_loc(1)+stim_radius) + mod(elapsed_time_s*tf/freqCenters(sf2), (slack-1)*2*stim_radius);
        mglBltTexture(tex2, [curr_x stim_loc(2)], 1, 0); % right-align the texture
    end
    mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], col_fix);
    mglFlush;
        
    if debug_t, toc; end;

  end
  
  begin_iti = clock; % see min_iti parameter above
  
  mglBltTexture(texml,[0 0]); mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], col_fix); 
    mglFlush;
  mglBltTexture(texml,[0 0]); mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], col_fix); 
    mglFlush;
  mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], col_fix);

  respKeys = 0 * mglGetKeys();
  KEYS_RESPONSE = [39 41]; % 'j' and 'k'
  % wait until you get the right input from the subject...
  while (sum(respKeys([KEYS_RESPONSE])) < 1), respKeys = mglGetKeys(); end

  IX_FEEDBACK_POS = 2;
  IX_FEEDBACK_NEG = 1;
  neg = 1;
  if (freqCenters(sf1) == freqCenters(sf2)) && randi(2) > 1.5 % if it's a "tie", play the sounds 50% of the time
      neg = 0;
  end
  if (respKeys(KEYS_RESPONSE(1)) == 1) && (freqCenters(sf1) > freqCenters(sf2))
      neg = 0;
  end
  if (respKeys(KEYS_RESPONSE(2)) == 1) && (freqCenters(sf2) > freqCenters(sf1))
      neg = 0;
  end
  
  if feedback && ~neg
      mglPlaySound(IX_FEEDBACK_POS)
      mglWaitSecs(0.500); % wait a bit so we actually get the sound...
  end
  
  fprintf(exp_info, '%g %g %d %g %g %d %d %g\n', freqCenters(sf1), con1, disp1, ...
      freqCenters(sf2), con2, disp2, find(respKeys==1), which_ref(tr_i));
  fprintf('tr %d: response: %g ...sf1 %g and sf2 %g...right? %d\n', tr_i, find(respKeys == 1), sf1, sf2, ~neg);
  
  if ~neg
    n_corr = n_corr + 1;
  end
  
  % before going on, make sure it's been at least "min_iti" seconds long
  while (etime(clock,begin_iti) < min_iti)
      continue;
  end
  
end
%%%
% End main loop.
%%%

% give the user feedback
mglStencilSelect(0); % no stencil...
mglTextSet('Helvectica', 40, [0 0 0]);
perCorr = mglText(sprintf('%.2f%% correct', 100*n_corr/NUM_TRIALS));
mglBltTexture(texml,[0 0]);
mglBltTexture(perCorr, [0 0], 'center', 'center');
mglFlush;

mglWaitSecs(10);

mglClose;
