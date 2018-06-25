%%%%%%%%%%
% SF Mixture piloting program. PL, 20180623.
% In the style of Georgeson, 1985 (i.e. simulatenous 2-AFC)

% This is an alternate formulation of the stimulus set where we choose the
% spacing between SF centers based on multiples of a Weber fraction, not
% logarithmic spacing over a defined range.
% This version will be more difficult, so use once the subject is trained.
%%%%%%%%%%

clear;
% useful functions...
myRound = @(x, digit) round((x.*10^digit))./10^digit;
scaleGrat = @(sf, con) (255/2)*(1+sf.*con);

%% Experiment parameters
% saving the file
subj = 1;
save_loc = 'data/'; meta_loc = 'data/metaData/';
is_pilot = 1;
run_num = 67;
save_meta = 1; % save metadata?

NUM_TRIALS = 150;
REF_DISP = 1;

stimDist = 4; % i.e. 6 degrees in periphy
xsgn = [-1 1]; ysgn = [-1, -1]; % i.e. pos or neg [x/y]
% [1st, 2nd] stimulus

debug_t = 0; % print how long each frame takes?
feedback = 0; % give feedback to the subject?

com_str = 'Pilot in VNL psychophysics room (1139) in dark conditions. Georgeson style (one interval). SF spacing based on Weber fraction.';
% comments for meta data file (e.g. room, condition, etc)
%% Design of experiment - (blank?:)stim1:blank:stim2
stim_dur = 0.32; % how many seconds for stimulus?
min_iti = 1; % the minimum inter-trial-interval is 1 second; wait after response if needed

% stimulus information
% sf
SF_REF = 6;
% SF_REF = 1.5;
weber_frac = 0.05; % 5 percent change
num_steps = 5; % N steps above, N steps below
incRef = 0; % include refStimulus in the list of test SF?
% tf 
tf = 3;
tf_spread = tf/5; % draw TF of dispersed gratings from gaussian with this sigma

if REF_DISP <= 3
    testDisps = REF_DISP; % [1 5];
    TEST_CONS = [0.04 0.08 0.16 0.32];
else
    testDisps = REF_DISP;
    TEST_CONS = [0.16 0.32 0.64 0.96];
end

REF_CON = TEST_CONS(end);

ori_all = 90; % vertical grating (drift in horizontal)

% stim. location and size
stim_radius = 1; % radius, in degrees
stim_loc1 = [xsgn(1)*sqrt((stimDist^2)/2), ysgn(1)*sqrt((stimDist^2)/2)];
stim_loc2 = [xsgn(2)*sqrt((stimDist^2)/2), ysgn(2)*sqrt((stimDist^2)/2)];
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
    save_base = sprintf('sfPer_s%02d_p%03dg', subj, run_num);
    if save_meta
        sm_base = sprintf('META_sfPer_s%02d_p%03dg', subj, run_num);
    end
else
    save_base = sprintf('sfPer_s%02d_%03dg', subj, run_num);
    if save_meta
        sm_base = sprintf('META_sfPer_s%02d_p%03dg', subj, run_num);
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
    fprintf(meta_inf, 'length of stimulus presentation in s: %.4f\n', stim_dur);
    fprintf(meta_inf, 'stimulus location x1: %.4f\n', stim_loc1(1));
    fprintf(meta_inf, 'stimulus location y1: %.4f\n', stim_loc1(2));
    fprintf(meta_inf, 'stimulus location x2: %.4f\n', stim_loc2(1));
    fprintf(meta_inf, 'stimulus location y2: %.4f\n', stim_loc2(2));
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
sf_round = 2; % just round to X digits...                                  
freqSeries = SF_REF*(1+(weber_frac*(-num_steps:num_steps)));
freqSeries = myRound(freqSeries, sf_round); % before removing 
ref_ind = find(freqSeries == SF_REF, 1, 'first'); % always the same; take only one index
% which indices/SFs are valid for test stimulus?
SFtestInds = 1:length(freqSeries);
if incRef == 0
    SFtestInds = setxor(SFtestInds, ref_ind);
end
    
% dispersion?
num_families = 5;
surr_oct = 1.5; % +- 1 octave relative to sfRef
num_gratings = 9; % how many samples around (and including) sfRef
octSeries  = linspace(surr_oct, -surr_oct, num_gratings);
sfVec = 2.^(octSeries(:) + log2(SF_REF));

% contrast
spreadVec = logspace(log10(.125), log10(1.25), num_families);
profTemp   = normpdf(octSeries, 0, spreadVec(REF_DISP));
profile    = profTemp/sum(profTemp);
if REF_DISP == 1 % if it's first family, we make it just 1 grating...
    profile = round(profile);
end
conVec = profile .* REF_CON;
%% Stencils
mglStencilCreateBegin(1);
mglFillOval(stim_loc1(1), stim_loc1(2), [2*stim_radius 2*stim_radius]);
mglFillOval(stim_loc2(1), stim_loc2(2), [2*stim_radius 2*stim_radius]);
mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], [1 1 1]); % fixation point
mglStencilCreateEnd;

% mglStencilCreateBegin(2);
% mglFillOval(stim_loc2(1), stim_loc2(2), [2*stim_radius 2*stim_radius]);
% mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], [1 1 1]); % fixation point
% mglStencilCreateEnd;
%% Create mean-lum texture and clear screen.
ml = 255*0.5*ones(yPix,xPix);
texml = mglCreateTexture(ml);
mglBltTexture(texml,[0 0]); mglFlush;

% determine order:
rng('shuffle');
which_ref = randi(2, NUM_TRIALS, 1);
test_ind = randi(length(SFtestInds), NUM_TRIALS, 1);
test_disp = testDisps(randi(length(testDisps), NUM_TRIALS, 1));
test_con = TEST_CONS(randi(length(TEST_CONS), NUM_TRIALS, 1));

if REF_DISP == 1 % Create single sinusoids
    % just fix
    ori = ori_all;
    for sf_c = 1 : length(freqSeries)
        ph = 360*rand(); % phase in degrees
        % slack*2*stim_radius because the stencil (i.e. stimulus area) is
        % 2*stim_radius wide/tall (diameter = 2*radius, ya!)
        grat{sf_c} = mglMakeGrating(slack*2*stim_radius, slack*2*stim_radius, freqSeries(sf_c), ori, ph);
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

  if which_ref(tr_i) == 1
      sf1 = ref_ind;
      con1 = REF_CON;
      disp1 = REF_DISP;
      if disp1 == 1
          tex1 = tex_ref;
      else % need to create random phases, tf; calculate con profile
          ph_grat1 = 360*rand(num_gratings, 1);
          tf_grat1 = random('norm', tf, tf_spread, [num_gratings, 1]);
          profTemp   = normpdf(octSeries, 0, spreadVec(disp1));
          profile1    = profTemp/sum(profTemp);
      end
          
      sf2 = SFtestInds(test_ind(tr_i));
      con2 = test_con(tr_i);
      disp2 = test_disp(tr_i);
      if disp2 == 1
          tex2 = mglCreateTexture(scaleGrat(grat{sf2}, con2));
      else % need to create random phases, tf; calculate con profile
          ph_grat2 = 360*rand(num_gratings, 1);
          tf_grat2 = random('norm', tf, tf_spread, [num_gratings, 1]);
          profTemp   = normpdf(octSeries, 0, spreadVec(disp2));
          profile2    = profTemp/sum(profTemp);
      end
  else
      sf1 = SFtestInds(test_ind(tr_i));
      con1 = test_con(tr_i);
      disp1 = test_disp(tr_i);
      if disp1 == 1
          tex1 = mglCreateTexture(scaleGrat(grat{sf1}, con1));
      else % need to create random phases, tf; calculate con profile
          ph_grat1 = 360*rand(num_gratings, 1);
          tf_grat1 = random('norm', tf, tf_spread, [num_gratings, 1]);
          profTemp   = normpdf(octSeries, 0, spreadVec(disp1));
          profile1    = profTemp/sum(profTemp);
      end
      
      sf2 = ref_ind;
      con2 = REF_CON;
      disp2 = REF_DISP;
      if disp2 == 1
        tex2 = tex_ref;
      else % need to create random phases, tf; calculate con profile
          ph_grat2 = 360*rand(num_gratings, 1);
          tf_grat2 = random('norm', tf, tf_spread, [num_gratings, 1]);
          profTemp   = normpdf(octSeries, 0, spreadVec(disp2));
          profile2    = profTemp/sum(profTemp);
      end
  end
  
  stim_on = clock;
  while (etime(clock,stim_on) < stim_dur)
 
    % first stim location
    mglStencilSelect(1);
    elapsed_time_s = etime(clock,stim_on);
    
    if disp1 > 1 % i.e. dispersed grating
        % need to calculate
        ori = ori_all;
        sfVec = 2.^(octSeries(:) + log2(freqSeries(sf1)));
        conVec = profile1 .* con1;
        for grat = 1 : num_gratings            
            curr_grat{grat} = scaleGrat(mglMakeGrating(2*stim_radius, ...
                2*stim_radius, sfVec(grat), ori, ...
                mod(ph_grat1(grat) - 360*elapsed_time_s*tf_grat1(grat), 360)), conVec(grat));
            if grat == 1
                curr_stim = curr_grat{grat} - 255*0.5; % subtract off mean luminance
            else
                curr_stim = curr_stim + (curr_grat{grat} - 255*0.5);
            end
        end
        curr_stim = curr_stim + 255*0.5;
        
        tex1 = mglCreateTexture(curr_stim);
        mglBltTexture(tex1, [stim_loc1(1) stim_loc1(2)], 0, 0); % center
    else
        curr_x = (stim_loc1(1)+stim_radius) + mod(elapsed_time_s*tf/freqSeries(sf1), (slack-1)*2*stim_radius);
        mglBltTexture(tex1, [curr_x stim_loc1(2)], 1, 0); % right-align the texture
    end
    
    mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], col_fix);
    mglFlush;
    
    if debug_t, toc; end;
        
    % second stim location
%     mglStencilSelect(2);
    elapsed_time_s = etime(clock,stim_on);
 
    if disp2 > 1 % i.e. dispersed grating
        % need to calculate
        ori = ori_all;
        sfVec = 2.^(octSeries(:) + log2(freqSeries(sf2)));
        conVec = profile2 .* con2;
        for grat = 1 : num_gratings            
            curr_grat{grat} = scaleGrat(mglMakeGrating(2*stim_radius, ...
                2*stim_radius, sfVec(grat), ori, ...
                mod(ph_grat2(grat) - 360*elapsed_time_s*tf_grat2(grat), 360)), conVec(grat));
            if grat == 1
                curr_stim = curr_grat{grat} - 255*0.5; % subtract off mean luminance
            else
                curr_stim = curr_stim + (curr_grat{grat} - 255*0.5);
            end
        end
        curr_stim = curr_stim + 255*0.5;
        
        tex2 = mglCreateTexture(curr_stim);
        mglBltTexture(tex2, [stim_loc2(1) stim_loc2(2)], 0, 0); % center
    else
        curr_x = (stim_loc2(1)+stim_radius) + mod(elapsed_time_s*tf/freqSeries(sf2), (slack-1)*2*stim_radius);
        mglBltTexture(tex2, [curr_x stim_loc2(2)], 1, 0); % right-align the texture
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
  % wait until you get the right input from the subject...0
  while (sum(respKeys([KEYS_RESPONSE])) < 1), respKeys = mglGetKeys(); end

  IX_FEEDBACK_POS = 2;
  IX_FEEDBACK_NEG = 1;
  neg = 1;
  if (freqSeries(sf1) == freqSeries(sf2)) && randi(2) > 1.5 % if it's a "tie", play the sounds 50% of the time
      neg = 0;
  end
  if (respKeys(KEYS_RESPONSE(1)) == 1) && (freqSeries(sf1) > freqSeries(sf2))
      neg = 0;
  end
  if (respKeys(KEYS_RESPONSE(2)) == 1) && (freqSeries(sf2) > freqSeries(sf1))
      neg = 0;
  end
  
  if feedback && ~neg
      mglPlaySound(IX_FEEDBACK_POS)
      mglWaitSecs(0.500); % wait a bit so we actually get the sound...
  end
  
  fprintf(exp_info, '%g %g %d %g %g %d %d %g\n', freqSeries(sf1), con1, disp1, ...
      freqSeries(sf2), con2, disp2, find(respKeys==1), which_ref(tr_i));
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
