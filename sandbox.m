%%%%%%%%%%
% Code for exploring contrast and sf perception
%%%%%%%%%%

clear;

% useful functions...
myRound = @(x, digit) round((x.*10^digit))./10^digit;
scaleGrat = @(sf, con) (255/2)*(1+sf.*con);

%% Experiment parameters
init_con = 0.16;
TEST_CONS = [0.01*2.^(0:6), 1];
con_ind = find(TEST_CONS == init_con);

init_sf = 3;
stim_oct = 2;
num_steps = 9;

REF_DISP = 5;

% center spatial frequencies
sf_round = 2; % just round to X digits...
lower_cent = 2^(log2(init_sf) - stim_oct);
higher_cent = 2^(log2(init_sf) + stim_oct);
freqSeries = myRound(logspace(log10(lower_cent), log10(higher_cent), num_steps), sf_round);

sf_ind = find(freqSeries == init_sf);

ori_all = 90; % use 90 degrees (vertical)

stimDist = 6; % i.e. 6 degrees in periphy
xsgn = 1; ysgn = -1; % i.e. pos or neg [x/y]

% stimulus information
tf = 5;
tf_spread = tf/5; % sigma will be this

% stim. location and size
stim_radius = 1; % radius, in degrees
stim_loc = [xsgn*sqrt((stimDist^2)/2), ysgn*sqrt((stimDist^2)/2)];
fp_radius = 0.1; % in degrees
col_fix = [1 1 1]; % for fixation point
slack = 2; % make the grating X times the size of the stencil/aperture...must be >1

refrac = 0.5; % X second(s) wait required for next adjustment following an initial one

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

%% stimulus creation/calculation
% dispersion?
num_families = 5;
surr_oct = 1.5; % +- 1 octave relative to sfRef
num_gratings = 9; % how many samples around (and including) sfRef
octSeries  = linspace(surr_oct, -surr_oct, num_gratings);
sfVec = 2.^(octSeries(:) + log2(init_sf));

% contrast
spreadVec = logspace(log10(.125), log10(1.25), num_families);
profTemp   = normpdf(octSeries, 0, spreadVec(REF_DISP));
profile    = profTemp/sum(profTemp);
if REF_DISP == 1 % if it's first family, we make it just 1 grating...
    profile = round(profile);
end
conVec = profile .* init_con;

% Stencils
mglStencilCreateBegin(1);
mglFillOval(stim_loc(1), stim_loc(2), [2*stim_radius 2*stim_radius]);
mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], [1 1 1]); % fixation point
mglStencilCreateEnd;

% Create mean-lum texture and clear screen.
ml = 255*0.5*ones(yPix,xPix);
texml = mglCreateTexture(ml);
mglBltTexture(texml,[0 0]); mglFlush;
%%%
%% Main stuff is here.
%%%

mglStencilSelect(0); % no stencil...
mglTextSet('Helvectica', 40, [0 0 0]);
instr = mglText(sprintf('a/d to lower/raise spatial frequency\ns/w to lower/raise contrast'));
mglBltTexture(texml,[0 0]);
mglBltTexture(instr, [0 0], 'center', 'center');
mglFlush;

mglBltTexture(texml,[0 0]);
mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], col_fix);
mglFlush;
% mglWaitSecs(5);
mglPlaySound('submarine');

mglBltTexture(texml,[0 0]); mglFlush;
mglBltTexture(texml,[0 0]); mglFlush;
mglStencilSelect(1);

% put the stimulus on the screen
% create the initial stimulus here
% slack*2*stim_radius because the stencil (i.e. stimulus area) is
% 2*stim_radius wide/tall (diameter = 2*radius, ya!)
if REF_DISP == 1
    ph = 360*rand(); % phase in degrees
    ori = ori_all;
    curr_grat = mglMakeGrating(slack*2*stim_radius, slack*2*stim_radius, init_sf, ori, ph);
    curr_tex = mglCreateTexture(scaleGrat(curr_grat, init_con));
else
    ori = ori_all;
    for grat = 1 : num_gratings
        ph_grat(grat) = 360*rand(); % phase in degrees
        tf_grat(grat) = random('norm', tf, tf_spread);
        curr_grat{grat} = scaleGrat(mglMakeGrating(2*stim_radius, ...
            2*stim_radius, sfVec(grat), ori, ph_grat(grat)), conVec(grat));
        if grat == 1
            curr_stim = curr_grat{grat} - 255*0.5; % subtract off mean luminance
        else
            curr_stim = curr_stim + (curr_grat{grat} - 255*0.5);
        end
    end
    curr_stim = curr_stim + 255*0.5;
    curr_tex = [];
end
    
KEYS_RESPONSE = mglCharToKeycode({'a', 'd', 's', 'w', '.'}); % see "instr" above

goOn = 1;
inter_on = clock;
last_change = -Inf.*ones(size(inter_on));
while (goOn)

%     tic;
    elapsed_time_s = etime(clock,inter_on);  

    if isempty(curr_tex) % i.e. dispersed grating...
        curr_tex = mglCreateTexture(curr_stim);
        mglBltTexture(curr_tex, [stim_loc(1) stim_loc(2)], 0, 0); % center
    else
        curr_x = (stim_loc(1)+stim_radius) + mod(elapsed_time_s*tf/freqSeries(sf_ind), (slack-1)*2*stim_radius);
        mglBltTexture(curr_tex, [curr_x stim_loc(2)], 1, 0); % right-align
    end
    mglPolygon(fp_radius*[-1 -1 1 1], fp_radius*[-1 1 1 -1], col_fix);
    mglFlush;

    respKeys = mglGetKeys();
    
    refracWait = etime(clock, last_change);
    if sum(respKeys(KEYS_RESPONSE)) == 1 && ((refracWait > refrac) || all(isinf(last_change)))
        % then user has changed the stim and we're ready to do it       
        fprintf('existing sf %.02f, con %.02f\n', freqSeries(sf_ind), TEST_CONS(con_ind));
        if respKeys(KEYS_RESPONSE(1)) == 1 % lower SF
            sf_ind = max(1, sf_ind - 1); % ensure we don't go below lowest sf_ind
        elseif respKeys(KEYS_RESPONSE(2)) == 1 % raise SF
            sf_ind = min(sf_ind + 1, length(freqSeries)); % ensure we don't go above highest sf_ind
        elseif respKeys(KEYS_RESPONSE(3)) == 1 % lower contrast
            con_ind = max(1, con_ind - 1); 
        elseif respKeys(KEYS_RESPONSE(4)) == 1 % raise contrast
            con_ind = min(con_ind + 1, length(TEST_CONS));
        elseif respKeys(KEYS_RESPONSE(5)) == 1 % exit
            goOn = 0;
        end
        
        fprintf('updated sf %.02f, con %.02f\n', freqSeries(sf_ind), TEST_CONS(con_ind));
        
        if REF_DISP == 1
            ph = 360*rand();
            curr_grat = mglMakeGrating(slack*2*stim_radius, slack*2*stim_radius, freqSeries(sf_ind), ori, ph);
            curr_tex = mglCreateTexture(scaleGrat(curr_grat, TEST_CONS(con_ind)));
        else % recalculate SF or CON vec
            ph = 360*rand(); % phase in degrees
            % RECALCULATE distributions
            sfVec = 2.^(octSeries(:) + log2(freqSeries(sf_ind)));
            conVec = profile .* TEST_CONS(con_ind);
            
            % Fresh gratings
            ori = ori_all;
            for grat = 1 : num_gratings
                curr_grat{grat} = scaleGrat(mglMakeGrating(2*stim_radius, ...
                    2*stim_radius, sfVec(grat), ori, ...
                    mod(ph_grat(grat) + 360*elapsed_time_s*tf_grat(grat), 360)), conVec(grat));
                if grat == 1
                    curr_stim = curr_grat{grat} - 255*0.5; % subtract off mean luminance
                else
                    curr_stim = curr_stim + (curr_grat{grat} - 255*0.5);
                end
            end
            curr_stim = curr_stim + 255*0.5;
            mglDeleteTexture(curr_tex); % clear mem/cache
            curr_tex = [];
        end
            
        inter_on = clock;
        last_change = clock;
    elseif REF_DISP > 1 % we update each grating. slow, but what we do for now
        ori = ori_all;
        for grat = 1 : num_gratings
            curr_grat{grat} = scaleGrat(mglMakeGrating(2*stim_radius, ...
                2*stim_radius, sfVec(grat), ori, ...
                mod(ph_grat(grat) + 360*elapsed_time_s*tf_grat(grat), 360)), conVec(grat));
            if grat == 1
                curr_stim = curr_grat{grat} - 255*0.5; % subtract off mean luminance
            else
                curr_stim = curr_stim + (curr_grat{grat} - 255*0.5);
            end
        end
        curr_stim = curr_stim + 255*0.5;
        mglDeleteTexture(curr_tex); % clear mem/cache
        curr_tex = [];
    end
    
%     toc
end

mglClose;
