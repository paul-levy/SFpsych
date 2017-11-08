%%%%%%%%%%
% SF Mixture piloting program. PL, 20171024.
% Initial code from LH
%%%%%%%%%%

% useful functions...
myRound = @(x, digit) round((x.*10^digit))./10^digit;
scaleGrat = @(sf, con) (255/2)*(1+sf.*con);

%% Experiment parameters
% saving the file
subj = 1;
save_loc = 'data/'; meta_loc = 'data/metaData/';
is_pilot = 1;
run_num = 17;
save_meta = 1; % save metadata?

NUM_TRIALS = 150;
REF_CON = 1;
REF_DISP = 1; % for now, just do 1 or 5...

stimDist = 6; % i.e. 6 degrees in periphy
xsgn = 1; ysgn = -1; % i.e. pos or neg [x/y]

debug_t = 0; % print how long each frame takes?
feedback = 0; % give feedback to the subject?

com_str = 'Pilot in VNL psychophysics room (1139) in dark conditions';
% comments for meta data file (e.g. room, condition, etc)

%% Design of experiment - (blank?:)stim1:blank:stim2
inter_1 = 0.2; % how many seconds for stimulus 1?
inter_2 = 0.2; % stimulus 2?
inter_blank = 1; % intervening blank?

% stimulus information
SF_REF = 3;
stim_oct = 1.25; % what separation (in octaves) +- for end points of center sf rel. to SF_REF?
num_steps = 9;
tf = 5;

TEST_CONS = [0.05 0.1 0.33 1];

if REF_DISP == 1
    testDisps = 1; % [1 5];
else
    testDisps = 5;
end

% stim. location and size
stim_radius = 2; % radius, in degrees
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

if isfile([save_loc, save_base]) % if it exists, append the time so we don't overwrite
    exp_info = fopen([save_loc, save_base, datestr(now), '.txt'], 'w+');
else
    exp_info = fopen([save_loc, save_base, '.txt'], 'w+');
end

if save_meta
    if isfile([meta_loc, sm_base])
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
sf_round = 2; % just round to X digits...
lower_cent = 2^(log2(SF_REF) - stim_oct);
higher_cent = 2^(log2(SF_REF) + stim_oct);
freqSeries = myRound(logspace(log10(lower_cent), log10(higher_cent), num_steps), sf_round);
% if you want one more point flanking either side ...
% of the reference (log space between adjacent and mid)
logMid = @(a, b) 2^((log2(a) + log2(b))/2);
ref_ind = find(freqSeries == SF_REF); % always the same...
lowMid = logMid(freqSeries(ref_ind-1), freqSeries(ref_ind));
highMid = logMid(freqSeries(ref_ind), freqSeries(ref_ind+1));
freqSeries = sort([freqSeries, lowMid, lowMid, lowMid, highMid, highMid, highMid]);

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
conVec = profile.*REF_CON;
if REF_DISP == 1 % if it's first family, we make it just 1 grating...
    conVec = round(conVec);
end

% Stencils
mglStencilCreateBegin(1);
mglFillOval(stim_loc(1), stim_loc(2), [stim_radius stim_radius]);
mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], [1 1 1]); % fixation point
mglStencilCreateEnd;

% Create mean-lum texture and clear screen.
ml = 255*0.5*ones(yPix,xPix);
texml = mglCreateTexture(ml);
mglBltTexture(texml,[0 0]); mglFlush;

% determine order:
rng('shuffle');
which_ref = randi(2, NUM_TRIALS, 1);
test_ind = randi(length(freqSeries), NUM_TRIALS, 1);
test_disp = testDisps(randi(length(testDisps), NUM_TRIALS, 1));
test_con = TEST_CONS(randi(length(TEST_CONS), NUM_TRIALS, 1));
ref_ind = find(freqSeries == SF_REF); % always the same...

if REF_DISP == 1 % Create single sinusoids
    % just fix
    ori = 90; ph = 0;
    for sf_c = 1 : length(freqSeries)
        grat{sf_c} = mglMakeGrating(slack*stim_radius, slack*stim_radius, freqSeries(sf_c), ori, ph);
        % tex{sf_c} = mglCreateTexture(scaleGrat(grat{sf_c}, TOTAL_CON);
    end
end
% create only the reference texture; others will be computed as needed
tex_ref = mglCreateTexture(scaleGrat(grat{ref_ind}, REF_CON));

%%%
% Main loop is here.
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
      tex1 = tex_ref;
      con1 = REF_CON;
      disp1 = REF_DISP;

      sf2 = test_ind(tr_i);
      con2 = test_con(tr_i);
      tex2 = mglCreateTexture(scaleGrat(grat{sf2}, con2));
      disp2 = test_disp(tr_i);
  else
      sf1 = test_ind(tr_i);
      con1 = test_con(tr_i);
      tex1 = mglCreateTexture(scaleGrat(grat{sf1}, con1));
      disp1 = test_disp(tr_i);
      
      sf2 = ref_ind;
      tex2 = tex_ref;
      con2 = REF_CON;
      disp2 = REF_DISP;
  end
  
  % first stim interval
  inter1 = clock;
  while (etime(clock,inter1) < inter_1)
 
    elapsed_time_s = etime(clock,inter1);

    if debug_t, tic; end; % for timing...
    elpsed = clock;
%     curr_x = stim_loc(1);
    curr_x = mod(stim_loc(1) + elapsed_time_s*tf/freqSeries(sf1), slack*stim_loc(1));
    mglBltTexture(tex1, [curr_x stim_loc(2)]);
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
    curr_x = stim_loc(1);
    curr_x = mod(stim_loc(1) + elapsed_time_s*tf/freqSeries(sf2), slack*stim_loc(1));
    mglBltTexture(tex2,[curr_x stim_loc(2)]);
    mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], col_fix);
    mglFlush;
    if debug_t, toc; end;

  end
%   mglDeleteTexture(tex2);
  
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
  if (sf1 == sf2) && randi(2) > 1.5 % if it's a "tie", play the sounds 50% of the time
      neg = 0;
  end
  if (respKeys(KEYS_RESPONSE(1)) == 1) && (sf1 > sf2)
      neg = 0;
  end
  if (respKeys(KEYS_RESPONSE(2)) == 1) && (sf2 > sf1)
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
