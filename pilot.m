%%%%%%%%%%
% SF Mixture piloting program. PL, 20171024.
% Initial code from LH
%%%%%%%%%%

% useful function...
myRound = @(x, digit) round((x.*10^digit))./10^digit;

% Constants.
debug_t = 0;
feedback = 1;

% Design of experiment - (blank?:)stim1:blank:stim2
NUM_TRIALS = 15;
inter_1 = 0.2; % how many seconds for stimulus 1?
inter_2 = 0.2; % stimulus 2?
inter_blank = 0.1; % intervening blank?

% stim. location and size
stim_radius = 0.5; % radius, in degrees
stim_loc = [1, -1];
fp_radius = 0.02; % in degrees
col_fix = [1 1 1]; % for fixation point
tf = 5; % set tf = cps
slack = 2; % make the grating X times the size of the stencil/aperture...

% center spatial frequencies
SF_REF = 5; % in cpd
sf_round = 2; % just round to X digits...
num_steps = 11;
stim_oct = 0.75; % +- centers will be spaced +/-X octave(s) apart
lower_cent = 2^(log2(SF_REF) - 1);
higher_cent = 2^(log2(SF_REF) + 1);
freqSeries = myRound(logspace(log10(lower_cent), log10(higher_cent), num_steps), sf_round);

% dispersion?
DISP_LEVEL = 1; % for now, just do 1 or 5...
num_families = 5;
surr_oct = 1.5; % +- 1 octave relative to sfRef
num_gratings = 9; % how many samples around (and including) sfRef
octSeries  = linspace(surr_oct, -surr_oct, num_gratings);
sfVec = 2.^(octSeries(:) + log2(SF_REF));

% contrast
TOTAL_CON = 1;
spreadVec = logspace(log10(.125), log10(1.25), num_families);
profTemp   = normpdf(octSeries, 0, spreadVec(DISP_LEVEL));
profile    = profTemp/sum(profTemp);
conVec = profile.*TOTAL_CON;
if DISP_LEVEL == 1 % if it's first family, we make it just 1 grating...
    conVec = round(conVec);
end

mon = 2; % in VNL psych room, 0 is "work" monitor, 2 is actual monitor

% Open, init MGL.
if mon == 2
    mglOpen(mon, 800, 800, 60, 32); %60Hz, bit depth = 32
else
    mglOpen(mon, 800, 800, 60, 32); % 800x800, 60Hz, bit depth = 32
end
mglVisualAngleCoordinates(57,[5 5]); % viewing from 57 cm, display will be in 5 x 5 cm box

% Stencils
mglStencilCreateBegin(1);
mglFillOval(stim_loc(1), stim_loc(2), [stim_radius stim_radius]);
mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], [1 1 1]); % fixation point
mglStencilCreateEnd;

% Create mean-lum texture and clear screen.
ml = 255*0.5*ones(800,800);
texml = mglCreateTexture(ml);
mglBltTexture(texml,[0 0]); mglFlush;

if DISP_LEVEL == 1 % Create single sinusoids
    % just fix
    ori = 90; ph = 0;
    for sf_c = 1 : length(freqSeries)
        grat = mglMakeGrating(slack*stim_radius, slack*stim_radius, freqSeries(sf_c), ori, ph);
        tex{sf_c} = mglCreateTexture(TOTAL_CON.*(255.*0.5.*(1+grat)));
    end
end

% determine order:
which_ref = randi(2, NUM_TRIALS, 1);
test_ind = randi(length(freqSeries), NUM_TRIALS, 1);
ref_ind = find(freqSeries == SF_REF); % always the same...


%%%
% Main loop is here.
%%%
mglWaitSecs(5);
mglPlaySound('submarine');

for tr_i = 1:NUM_TRIALS

  mglBltTexture(texml,[0 0]); mglFlush;
  mglBltTexture(texml,[0 0]); mglFlush;
  mglStencilSelect(1);

  if which_ref(tr_i) == 1
      sf1 = ref_ind;
      sf2 = test_ind(tr_i);
  else
      sf1 = test_ind(tr_i);
      sf2 = ref_ind;
  end
  
  % first stim interval
  inter1 = clock;
  while (etime(clock,inter1) < inter_1)
 
    elapsed_time_s = etime(clock,inter1);

    if debug_t, tic; end; % for timing...
    elpsed = clock;
%     curr_x = stim_loc(1);
    curr_x = mod(stim_loc(1) + elapsed_time_s*tf/freqSeries(sf1), slack*stim_loc(1));
    mglBltTexture(tex{sf1}, [curr_x stim_loc(2)]);
    mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], col_fix);
%     mglPolygon([-2.5 -2.5 -2.0 -2.0], [2.0 2.5 2.5 2.0], [1 1 1]);
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
    mglBltTexture(tex{sf2},[curr_x stim_loc(2)]);
    mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], col_fix);
    mglFlush;
    if debug_t, toc; end;

  end
  
  mglBltTexture(texml,[0 0]); mglFlush;
  mglBltTexture(texml,[0 0]); mglFlush;
  mglPolygon(0.5*[-1 -1 1 1], 0.5*[-1 1 1 -1], col_fix);

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
      mglWaitSecs(0.250); % wait a bit so we actually get the sound...
  end
        
  fprintf('tr %d: response: %g ...sf1 %g and sf2 %g...right? %d\n', tr_i, find(respKeys == 1), sf1, sf2, ~neg);
  
end
%%%
% End main loop.
%%%

mglClose;
