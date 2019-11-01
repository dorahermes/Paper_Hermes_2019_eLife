%
% Filter patterns with gabor patches for:
%
% Hermes D, Petridou N, Kay K, Winawer J. 2019 An image-computable model
% for the stimulus selectivity of gamma oscillations. eLife 2019;8:e47035.
% DOI: https://doi.org/10.7554/eLife.47035
%
%
% The was adapted from the code used in: 
%
% Kay KN, Winawer J, Rokem A, Mezer A, Wandell BA (2013) A Two-Stage
% Cascade Model of BOLD Responses in Human Visual Cortex. PLoS Comput Biol
% 9(5): e1003079. https://doi.org/10.1371/journal.pcbi.1003079
%
% Dora Hermes, 2017

clear all

% set paths:
rootPath = gammaModelPath();
dataDir = fullfile(rootPath,'data');

%%
%% %%%%%% START preprocess images %%%%%%%%
%%

% load images
load(fullfile(dataDir,'stimuli','task-soc_stimuli.mat'),'stimuli')


%% %%%%%% preprocess images - part 1 is fast %%%%%%%%

%%%% DOWNSAMPLE TO INCREASE SPEED
% resize the stimuli to 240 x 240 to reduce computational time.
% use single-format to save memory.
temp = zeros(240,240,size(stimuli,3),'single');
for p=1:size(stimuli,3)
  temp(:,:,p) = imresize(single(stimuli(:,:,p)),[240 240],'cubic');
end
stimulus = temp;
clear temp;

%%%% RESCALE
% ensure that all values are between 0 and 254.
% rescale values to the range [0,1].
% subtract off the background luminance (0.5).
% after these steps, all pixel values will be in the
% range [-.5,.5] with the background corresponding to 0.
stimulus(stimulus < 0) = 0;
stimulus(stimulus > 254) = 254;
stimulus = stimulus/254 - 0.5;

%%%% ZERO PAD
% pad the stimulus with zeros (to reduce edge effects).
% the new resolution is 270 x 250 (15-pixel padding on each side).
stimulus = placematrix(zeros(270,270,size(stimulus,3),'single'),stimulus);

% inspect one of the stimuli
figure;
imagesc(stimulus(:,:,15));
axis image tight;
caxis([-.5 .5]);
colormap(gray);
colorbar;
title('Stimulus');

%% %%%%%% preprocess images - part 2 is timeconsuming %%%%%%%%

% Apply Gabor filters to the stimuli.  filters occur at different positions,
% orientations, and phases.  there are several parameters that govern the
% design of the filters:
filt_prop.cycles = 60*(270/240);    %   the number of cycles per image is 60*(270/240)
filt_prop.bandwidth = -1;           %   the spatial frequency bandwidth of the filters is 1 octave
filt_prop.spacings=1;               %   the separation of adjacent filters is 1 std dev of the Gaussian envelopes
                                    %     (this results in a 135 x 135 grid of positions)
filt_prop.orientations=8;           %   filters occur at 8 orientations
filt_prop.phases=2;                 %   filters occur at 2 phases (between 0 and pi)
filt_prop.thres=0.01;               %   the Gaussian envelopes are thresholded at .01
filt_prop.scaling=2;                %   filters are scaled to have an equivalent Michelson contrast of 1
filt_prop.mode=0;                   %   the dot-product between each filter and each stimulus is computed

% after this step, stimulus is images x phases*orientations*positions.
stimulus = applymultiscalegaborfilters(reshape(stimulus,270*270,[])', ...
  filt_prop.cycles,filt_prop.bandwidth,filt_prop.spacings,filt_prop.orientations,...
  filt_prop.phases,filt_prop.thres,filt_prop.scaling,filt_prop.mode);

save(fullfile(dataDir,'derivatives','gaborFilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

%% %%%%%% preprocess images - part 3 is fast %%%%%%%%

load(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

% compute the square root of the sum of the squares of the outputs of
% quadrature-phase filter pairs (this is the standard complex-cell energy model).
% after this step, stimulus is images x orientations*positions.
stimulus = sqrt(blob(stimulus.^2,2,2));

% compute the population term in the divisive-normalization equation.
% this term is simply the average across the complex-cell outputs
% at each position (averaging across orientation).
stimulusPOP = blob(stimulus,2,8)/8;

% repeat the population term for each of the orientations
stimulusPOP = upsamplematrix(stimulusPOP,8,2,[],'nearest');

% apply divisive normalization to the complex-cell outputs.  there are two parameters
% that influence this operation: an exponent term (r) and a semi-saturation term (s).
% the parameter values specified here were determined through a separate fitting
% procedure (see paper for details).  for the purposes of this script, we will
% simply hard-code the parameter values here and not worry about attempting to fit
% the parameters.
r = 1;
s = 0.5;
stimulus = stimulus.^r ./ (s.^r + stimulusPOP.^r);
clear stimulusPOP;

% sum across orientation.  after this step, stimulus is images x positions.
stimulus = blob(stimulus,2,8);
save(fullfile(dataDir,'derivatives','gaborFilt','task-soc_stimuli_gaborFilt02.mat'),'stimulus')

%%
% inspect one of the stimuli
figure;
mx = max(abs(stimulus(:)));
imagesc(reshape(stimulus(12,:),[135 135]));
axis image tight;
caxis([0 mx]);
colormap(gray);
colorbar;
title('Stimulus');
    
%%
%% %%%%%% END preprocess images %%%%%%%%
%%





