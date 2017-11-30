function behavior = csfCore(dataLoc,subjectID,distanceViewing_m,sfGrating_cpd,stepStart,numTrialsTot)
%%%
% v1.0 20160617 Luke Hallum
% v1.1 2017???? Gerick Lee
% v1.2 20171130 Paul Levy
%%%%%%%%%%
global MGL;
mglListener('init');

%%%
%% Experiment parameters and constants.
%%%%%%%%%%
fp_pix_radius = 4;
feedback = 1;
save_meta = 1; 

NUM_EASY_TRIALS = 0;%5;
ECCENTRICITY_STIMULI_DEG = 4;
WIDTH_STIMULUS_DEG = 2;
POSITIONS_STIMULI = [7 3 11]; % 0 is east, 1 has angle of elevation = 30deg, etc.
                              % 0-indexed from a cartesian heading of 0 deg
                              % order pertains to subjects responses 'j', 'k', 'l'
RES_DISPLAY_CM = [30 40];
xPix = 1280; yPix = 960;
RES_DISPLAY_PIX = [yPix xPix];
RES_BACKGROUND_PIX = RES_DISPLAY_PIX;
DURATION_BLANK1_S = 0.4;
DURATION_BLANK2_S = 0.4;
DURATION_PROBE_S = 0.15;
% DURATION_PROBE_S = 1; % for pilot foveal testing (see filename, too)

IX_FEEDBACK_POS = 2;
IX_FEEDBACK_NEG = 1;
TEST_CON_PERCENT = power(2, [-3 : 0.5 : 6.5]);
% [0.125 0.1875 0.25 0.375 0.5 0.75 1 1.5 2 3 4 6 8 12 16 24 32 48 64 96];
%%%%%%%%%%

com_str = '1139, dark room';

%%%
% Seed rand -- see 'help rand'.
%%%%%%%%%%
rand('twister',sum(100*clock));
%%%%%%%%%%
%%%
%% Convert stim positions to coordinates in degrees. Also, get coordinates for fiducial markers.
%%%%%%%%%%
x_stim_deg = zeros(1,length(POSITIONS_STIMULI));
y_stim_deg = x_stim_deg;
for iipos = 1:length(POSITIONS_STIMULI)
  x_stim_deg(iipos) = real(ECCENTRICITY_STIMULI_DEG*exp(i*POSITIONS_STIMULI(iipos)*pi/6));
  y_stim_deg(iipos) = imag(ECCENTRICITY_STIMULI_DEG*exp(i*POSITIONS_STIMULI(iipos)*pi/6));
end
x_fiducial_deg = zeros(1,length(POSITIONS_STIMULI));
y_fiducial_deg = x_fiducial_deg;
for iipos = 1:length(POSITIONS_STIMULI)
  x_fiducial_deg(iipos) = real((ECCENTRICITY_STIMULI_DEG-0.66*WIDTH_STIMULUS_DEG)*exp(i*POSITIONS_STIMULI(iipos)*pi/6));
  y_fiducial_deg(iipos) = imag((ECCENTRICITY_STIMULI_DEG-0.66*WIDTH_STIMULUS_DEG)*exp(i*POSITIONS_STIMULI(iipos)*pi/6));
end
%%%%%%%%%%

%%%
% Make a coordinate system for gratings. The _resolution_ of this coordinate
% system is unchanged from block to block. The unit of measure here is
% 'cycles', which depends upon the spatial frequency being tested in this block
% of trials -- an input arg to the function. We'll use the DISK for windowing
% gratings -- its edges are smooth.
%%%%%%%%%%
resGrating_pix = repmat(distanceViewing_m/0.57*WIDTH_STIMULUS_DEG/RES_DISPLAY_CM(1)*RES_DISPLAY_PIX(1),[1 2]);
[x,y] = meshgrid(linspace(-sfGrating_cpd*WIDTH_STIMULUS_DEG/2, sfGrating_cpd*WIDTH_STIMULUS_DEG/2, resGrating_pix(1)));
y = flipud(y);
%
[xx,yy] = meshgrid(linspace(-0.5,0.5,resGrating_pix(1)));
rr = sqrt(xx.^2 + yy.^2);
DISK = double(rr < 0.5) .* (0.5 + 0.5*cos(2*pi*1/0.2*(max(0.4, rr) - 0.4)));
%REC_APERTURE = double((xx > -0.3125) & (xx < 0.3125) & (yy > -0.1875/2) & (yy < 0.1875/2));
%REC_APERTURE = conv2(REC_APERTURE, ones(10,10)/100, 'same');
%%%%%%%%%%

%%%
%% Information for saving files; save metadata if requested
saveDir = dataLoc;
metaDir = [saveDir, 'metaData/'];

filename = sprintf('csf_oddball_%s_sf%s', subjectID, sprintf('%04.0f', round(1e2*sfGrating_cpd)));
if isfile([saveDir, filename, '.txt'])
    filename = [filename, datestr(now)];
end

if save_meta
    meta_fn = ['META_' filename];

    if isfile([metaDir, meta_fn, '.txt'])
        meta_inf = fopen([metaDir, meta_fn, datestr(now), '.txt'], 'w+');
    else
        meta_inf = fopen([metaDir, meta_fn, '.txt'], 'w+');
    end
    
    fprintf(meta_inf, 'number of trials: %g\nstimulus eccentricity in degrees: %.3f\n', numTrialsTot, ECCENTRICITY_STIMULI_DEG);
    fprintf(meta_inf, 'length of blank pre-interval in s: %.4f\n', DURATION_BLANK1_S);
    fprintf(meta_inf, 'length of blank post-interval in s: %.4f\n', DURATION_BLANK2_S);
    fprintf(meta_inf, 'length of stimulus interval in s: %.4f\n', DURATION_PROBE_S);
    fprintf(meta_inf, 'stimulus location x1: %.4f\n', x_stim_deg(1));
    fprintf(meta_inf, 'stimulus location y1: %.4f\n', y_stim_deg(1));
    fprintf(meta_inf, 'stimulus location x2: %.4f\n', x_stim_deg(2));
    fprintf(meta_inf, 'stimulus location y2: %.4f\n', y_stim_deg(2));
    fprintf(meta_inf, 'stimulus location x3: %.4f\n', x_stim_deg(3));
    fprintf(meta_inf, 'stimulus location y3: %.4f\n', y_stim_deg(3));
    fprintf(meta_inf, 'stimulus radius degrees: %.3f\n', WIDTH_STIMULUS_DEG/2);
    fprintf(meta_inf, 'fixation cross radius in pixels: %.3f\n', fp_pix_radius);
    fprintf(meta_inf, 'screen resolution (x, y) in pixels: %g, %g\n', xPix, yPix);
    fprintf(meta_inf, 'screen size (x, y) in cm: %g, %g\n', RES_DISPLAY_CM(2), RES_DISPLAY_CM(1));
    fprintf(meta_inf, 'viewing distance in cm: %g\n', distanceViewing_m*1e2);
    fprintf(meta_inf, 'subject feedback? %g\n', feedback);
    fprintf(meta_inf, 'start time: %s\n', datestr(now));
    fprintf(meta_inf, 'comments: %s\n', com_str);
    fclose(meta_inf);
end
%% Design matrix.
%%%%%%%%%%
fnShuffle = @(x) x(randperm(length(x)));
% pseudorandom
vectorDesignIXPositionOddball = fnShuffle(repmat(1:length(POSITIONS_STIMULI), [1 ceil(numTrialsTot/length(POSITIONS_STIMULI))]));
vectorDesignIXPositionOddball = vectorDesignIXPositionOddball(1:numTrialsTot);
vectorDesignOddballIsGrating = fnShuffle(repmat([0 1], [1 ceil(numTrialsTot/2)]));
vectorDesignOddballIsGrating = vectorDesignOddballIsGrating(1:numTrialsTot);
%%%%%%%%%%

%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%
%% Set up display
mglOpen();
mglVisualAngleCoordinates(100*distanceViewing_m, fliplr(RES_DISPLAY_CM));
%%%
% Gamma correction - completed 9/16 by GL, checked Fall, 2017
%%%%
EXP_R = 1/2.3307;
EXP_G = 1/2.3724;
EXP_B = 1/2.3787;
redGammaTable = (0:1/255:1).^EXP_R;
greenGammaTable = (0:1/255:1).^EXP_G;
blueGammaTable = (0:1/255:1).^EXP_B;
mglSetGammaTable(redGammaTable,greenGammaTable,blueGammaTable);
%%%%%%%%%%

%%%
%% Create mean luminance, etc.
%%%%%%%%%%
fieldMeanLumAndFixation = 0.5*ones(RES_BACKGROUND_PIX);
% vertical portion of cross
fieldMeanLumAndFixation(floor(yPix/2-fp_pix_radius):ceil(yPix/2+fp_pix_radius),floor(xPix/2)) = 0;
% horizontal portion of cross
fieldMeanLumAndFixation(floor(yPix/2),floor(xPix/2-fp_pix_radius):ceil(xPix/2+fp_pix_radius)) = 0;
fieldMeanLumAndFixation2 = -1*(fieldMeanLumAndFixation - 0.5) + 0.5;

mglgray = mglCreateTexture(round(255.0*fieldMeanLumAndFixation));
mglgray2 = mglCreateTexture(round(255.0*fieldMeanLumAndFixation2));
mglBltTexture(mglgray2,[0 0]); mglFlush; % Set both display buffers to gray.
mglBltTexture(mglgray2,[0 0]); mglFlush; %
%%%%%%%%%%

%%%%%%%%%%
%%%%%%%%%%
%%%
%% The format of matrix 'behavior'.
%%%%%%%%%%
IX_BEH_LEVEL = 1;
IX_BEH_JKL = 2;
IX_BEH_CORRECT = 3;
IX_BEH_ODDBALL_IS_GRATING = 4;
IX_BEH_POSITION_ODDBALL = 5;
behavior = [-1 -1 -1 -1 -1];
%%%%%%%%%%

%%%
%% Start the experiment: wait for key press....
%%%%%%%%%%
mglBltTexture(mglgray2,[0 0]); mglFlush();
KEYS_RESPONSE = [39 41 38]; % j, k, l
%respKeys = 0 * mglGetKeys();
%disp('Awaiting key press...')
%while (sum(respKeys([KEYS_RESPONSE])) < 1), respKeys = mglGetKeys(); end
%%%%%%%%%%
for iitrial = 1:length(vectorDesignIXPositionOddball)

  %%%
  % Get contrast of the grating for this trial, then make the grating.
  %%%%%%%%%%
  if (iitrial == 1), this_step = stepStart; end
  if (iitrial > NUM_EASY_TRIALS)
    if (behavior(iitrial,IX_BEH_CORRECT) == 1), this_step = max(1,this_step-1); end
    if (behavior(iitrial,IX_BEH_CORRECT) == 0), this_step = min(length(TEST_CON_PERCENT),this_step+1); end
  end
  this_con_percent = TEST_CON_PERCENT(this_step);
  if (iitrial <= NUM_EASY_TRIALS), this_con_percent = 80; end
  this_grating = this_con_percent/100*cos(2*pi*(y-rand));
  this_grating = DISK .* this_grating;
  %%%%%%%%%%

  %%%
  % Create textures. Organize.
  %%%%%%%%%%
  % Target goes into 'this_stim1_'.
  if (vectorDesignOddballIsGrating(iitrial) == 1)
    this_stim1_ = this_grating;
    this_stim2_ = zeros(size(x));
    this_stim3_ = zeros(size(x));
  end
  if (vectorDesignOddballIsGrating(iitrial) == 0)
    this_stim1_ = zeros(size(x));
    this_stim2_ = this_grating;
    this_stim3_ = this_grating;
  end
  % Down-sample in preparation for horizontal dithering.
  this_stim1_ = this_stim1_(:,1:3:size(x,2));
  this_stim2_ = this_stim2_(:,1:3:size(x,2));
  this_stim3_ = this_stim3_(:,1:3:size(x,2));
  % Target into stimulus site 1, 2 or 3. Others get distractors, in no particular order.
  if (vectorDesignIXPositionOddball(iitrial) == 1)
    this_stim1 = mglCreateTexture(double(csfDitherH(255.0*(0.5 + 0.5*this_stim1_))));
    this_stim2 = mglCreateTexture(double(csfDitherH(255.0*(0.5 + 0.5*this_stim2_))));
    this_stim3 = mglCreateTexture(double(csfDitherH(255.0*(0.5 + 0.5*this_stim3_))));
  end
  if (vectorDesignIXPositionOddball(iitrial) == 2)
    this_stim1 = mglCreateTexture(double(csfDitherH(255.0*(0.5 + 0.5*this_stim3_))));
    this_stim2 = mglCreateTexture(double(csfDitherH(255.0*(0.5 + 0.5*this_stim1_))));
    this_stim3 = mglCreateTexture(double(csfDitherH(255.0*(0.5 + 0.5*this_stim2_))));
  end
  if (vectorDesignIXPositionOddball(iitrial) == 3)
    this_stim1 = mglCreateTexture(double(csfDitherH(255.0*(0.5 + 0.5*this_stim2_))));
    this_stim2 = mglCreateTexture(double(csfDitherH(255.0*(0.5 + 0.5*this_stim3_))));
    this_stim3 = mglCreateTexture(double(csfDitherH(255.0*(0.5 + 0.5*this_stim1_))));
  end
  %%%%%%%%%%
  
  %% dithering test section:
%   t3 = double(csfDitherH(255.0*(0.5 + 0.5*this_stim1_)));
%   t3(t3 ~= 127) = 127;
% %   t3(t3 ~= 128) = 128;
% %   t3 = zeros(size(t3));
%   t3(:, 1 : 3 : 258) = 128;
%   t3(:, 2 : 3 : 258) = 128;
% %   t3([1 : 10, end - 9 : end], [1 end]) = 0;
%   mglblack = mglCreateTexture(zeros(900));
% 
%   this_stim1 = mglCreateTexture(t3);
%   this_stim2 = mglCreateTexture(zeros(size(t3)));
%   this_stim3 = mglCreateTexture(zeros(size(t3)));
% %   mglBltTexture([mglblack; this_stim1; this_stim2; this_stim3],[0 0; 0 0; 90 90; 90 90]); mglFlush();
%   mglBltTexture([mglblack; this_stim1],[0 0; 0 0]); mglFlush();
  % end dithering test section
  %%

  %%%
  % Wait for key press then render.
  %%%%%%%%%%
  respKeys = 0 * mglGetKeys();
  while (sum(respKeys([KEYS_RESPONSE])) < 1), respKeys = mglGetKeys(); end
  mglBltTexture(mglgray,[0 0]); mglFlush;
  if (DURATION_BLANK1_S > 0)
    mglBltTexture(mglgray,[0 0]); mglFlush(); mglWaitSecs(DURATION_BLANK1_S);
  end
  mglBltTexture([mglgray; this_stim1; this_stim2; this_stim3],[0 0; x_stim_deg(1) y_stim_deg(1); x_stim_deg(2) y_stim_deg(2); x_stim_deg(3) y_stim_deg(3)]); mglFlush(); mglWaitSecs(DURATION_PROBE_S+double(iitrial<NUM_EASY_TRIALS));
  mglBltTexture(mglgray,[0 0]); mglFlush(); mglWaitSecs(DURATION_BLANK2_S);
  %%%%%%%%%%
  %%%
  % Get behavior.
  %%%%%%%%%%
  respKeys = 0 * mglGetKeys();
  while (sum(respKeys([KEYS_RESPONSE])) < 1), respKeys = mglGetKeys(); end
  behavior(end+1,IX_BEH_LEVEL) = this_con_percent;
  if (respKeys(KEYS_RESPONSE(1)) == 1), behavior(end,IX_BEH_JKL) = 1; end
  if (respKeys(KEYS_RESPONSE(2)) == 1), behavior(end,IX_BEH_JKL) = 2; end
  if (respKeys(KEYS_RESPONSE(3)) == 1), behavior(end,IX_BEH_JKL) = 3; end
  behavior(end,IX_BEH_CORRECT) = 0;
  if ((vectorDesignIXPositionOddball(iitrial) == 1) & (behavior(end,IX_BEH_JKL) == 1)), behavior(end,IX_BEH_CORRECT) = 1; end
  if ((vectorDesignIXPositionOddball(iitrial) == 2) & (behavior(end,IX_BEH_JKL) == 2)), behavior(end,IX_BEH_CORRECT) = 1; end
  if ((vectorDesignIXPositionOddball(iitrial) == 3) & (behavior(end,IX_BEH_JKL) == 3)), behavior(end,IX_BEH_CORRECT) = 1; end
  behavior(end,IX_BEH_ODDBALL_IS_GRATING) = vectorDesignOddballIsGrating(iitrial);
  behavior(end,IX_BEH_POSITION_ODDBALL) = POSITIONS_STIMULI(vectorDesignIXPositionOddball(iitrial));
  %%%%%%%%%%

  %%%
  % Feedback.
  %%%%%%%%%%
  if feedback
      if(behavior(end,IX_BEH_CORRECT) > 0), mglPlaySound(IX_FEEDBACK_POS);
      else mglPlaySound(IX_FEEDBACK_NEG);
      end
  end
  pause(0.25);
  mglBltTexture(mglgray2,[0 0]); mglFlush();
  %%%%%%%%%%

  mglDeleteTexture(this_stim1);
  mglDeleteTexture(this_stim2);
  mglDeleteTexture(this_stim3);
end

%% Save and close out of MGL
behavior = behavior(2:end,:);
mglClose;
dlmwrite(sprintf('%s.txt',[saveDir, filename]), behavior, 'delimiter',' ')

return;

%%%
% Functions...
%%%%%%%%%%

function x_ = csfDitherH(x)
% This function dithers horizontally. The input image is floating point and
% ranges on [0.0, 255.0]; its resolution, originally [nrow ncol], grows to
% [nrow 3*ncol]; the output is quantized to integers on [0, 255].

  x_ = uint8(floor(kron(x, [1 1 1])));
  x_(:,1:3:size(x_,2)) = x_(:,1:3:size(x_,2)) + uint8(round(3*rem(x,1)) > 0);
  x_(:,2:3:size(x_,2)) = x_(:,2:3:size(x_,2)) + uint8(round(3*rem(x,1)) > 1);
  x_(:,3:3:size(x_,2)) = x_(:,3:3:size(x_,2)) + uint8(round(3*rem(x,1)) > 2);

return;

