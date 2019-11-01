
% This script will generate the panels of Figure 8 from: 
%
% Hermes D, Petridou N, Kay K, Winawer J. 2019 An image-computable model
% for the stimulus selectivity of gamma oscillations. eLife 2019;8:e47035.
% DOI: https://doi.org/10.7554/eLife.47035
%
% Dora Hermes, 2019

clear all

% set paths:
rootPath = gammaModelPath();
dataDir = fullfile(rootPath,'data');

%% Load original images
origIm = load(fullfile(dataDir,'stimuli','task-soc_stimuli.mat'),'stimuli');

%% Load preprocessed images

load(fullfile(dataDir,'derivatives','gaborFilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

% compute the square root of the sum of the squares of the outputs of
% quadrature-phase filter pairs (this is the standard complex-cell energy model).
% after this step, stimulus is images x orientations*positions.
stimulus = sqrt(blob(stimulus.^2,2,2));
%%%-----> this is this input to the OV model

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
stimulusSOC = stimulus.^r ./ (s.^r + stimulusPOP.^r);
clear stimulusPOP;

% sum across orientation.  after this step, stimulus is images x positions.
imEnergyMean = blob(stimulusSOC,2,8);
%%%-----> this is this input to the SOC model

res = sqrt(size(stimulus,2)/8);

%% Display NBF model results for one electrode
% Figure 8C-D

%%%%% Pick a subject:
subjects = {'19','24','1001'};
s = 2; 
subj = subjects{s};
elec = 45;% used sub-24 electrode 45

analysisType = 'spectra200';
modelType = 'OVmodel';

% Load NBF model fit:
load(fullfile(dataDir,'derivatives','gaborFilt','OV_gamma',...
    ['sub' subj '_el' int2str(elec) '_' analysisType '_' modelType]),...
    'seed_params','cross_OVparams','cross_OVestimate','ov_exponents','train_OVperformance')
ov_exp_used = find(ov_exponents==.5);

% Load ECoG data:
dataFitName = fullfile(dataDir,'derivatives','preprocessing',...
    ['sub-' subj '_ses-01_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
load(dataFitName)

% Calculate ECoG gamma power percent signal change:
bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
ecog_bb = mean(100*(10.^(resamp_parms(:,:,2)-bb_base)-1),2);
ecog_bb_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,2),.16,2)) ...
    squeeze(quantile(resamp_parms(:,:,2),.84,2))]'-bb_base)-1);
ecog_g = mean(100*(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1),2);
ecog_g_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.16,2)) ...
    squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.84,2))]')-1);
ylims = [min(ecog_g_err(:)) max([ecog_g_err(:)])];

res = sqrt(size(imEnergyMean,2));  % resolution of the pre-processed stimuli
[~,xx,yy] = makegaussian2d(res,2,2,2,2);

% Define a Gaussian centered at the prf:
gaufun1 = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),...
    pp(3),xx,yy,0,0)/(2*pi*pp(3)^2));

im_nrs = [10 74:78];

% Plot filtered (contrast) images with prf
figure('Position',[0 0 800 600])
max_prf_images = 0.05;
for kk = 1:length(im_nrs)
    % images + prf
    subplot(3,length(im_nrs),kk)
    imagesc(double(origIm.stimuli(:,:,im_nrs(kk))),[0 max(origIm.stimuli(:))]),hold on
    % Move pRF in 135x135 to original size of 800x800 
    %   half of the screen is 60 in filtered and 400 in original
    orig_x = 400 + (seed_params(1)-res/2) * 400./60;
    orig_y = 400 + (seed_params(2)-res/2) * 400./60;
    orig_sigma = seed_params(3) * 400./60;
    axis image

    numPoints = 50;
    c.th = linspace(0,2*pi, numPoints);
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*orig_sigma);
    plot(c.x + orig_y, c.y + orig_x, 'y','LineWidth',2) % this is just reversed because plot and imagesc are opposite, checked this with contour
    [c.x, c.y] = pol2cart(c.th, 2*ones(1,numPoints)*orig_sigma);
    plot(c.x + orig_y, c.y + orig_x, 'y:','LineWidth',2) % this is just reversed because plot and imagesc are opposite, checked this with contour
    axis off
    xlim([orig_y-3*orig_sigma orig_y+3*orig_sigma])
    ylim([orig_x-3*orig_sigma orig_x+3*orig_sigma])
    title(['im ' int2str(im_nrs(kk))])
    
%     % contrast images + prf
%     subplot(3,length(im_nrs),length(im_nrs)+kk)
%     imagesc(reshape(imEnergyMean(im_nrs(kk),:),res,res),[0 2]);
%     axis image
%     hold on
%     numPoints = 50;
%     c.th = linspace(0,2*pi, numPoints);
%     [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*seed_params(3));
%     plot(c.x + seed_params(2), c.y + seed_params(1), 'r') % this is just reversed because plot and imagesc are opposite, checked this with contour
%     xlim([seed_params(2)-3*seed_params(3) seed_params(2)+3*seed_params(3)])
%     ylim([seed_params(1)-3*seed_params(3) seed_params(1)+3*seed_params(3)])
%     
    colormap gray
end

%%% Bar plot gamma
subplot(3,length(im_nrs),length(im_nrs)*2+[1 2]),hold on
bar([1:length(im_nrs)],ecog_g(im_nrs),1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0]);
plot([1:length(im_nrs); 1:length(im_nrs)],ecog_g_err(:,im_nrs),'k');
% plot([1:length(im_nrs)],cross_OVestimate(im_nrs,ov_exp_used),'r','LineWidth',1)
plot([1:length(im_nrs)],cross_OVestimate(im_nrs,ov_exp_used),'r.','MarkerSize',20)
% ylim([0 max(max(ecog_g_err(:,im_nrs)))])
xlim([0 length(im_nrs)+1])
set(gca,'XTick',1:length(im_nrs),'XTickLabel',im_nrs)


set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',fullfile(dataDir,'derivatives','figures',...
        ['Figure8d_sub-' subj '_el' int2str(elec)]))
print('-dpng','-r300',fullfile(dataDir,'derivatives','figures',...
        ['Figure8d_sub-' subj '_el' int2str(elec)]))

%% plot large images so we can zoom into the pRF

max_prf_images = 0.05;
im_nrs = [10 74:78];

for kk = 1:length(im_nrs)
    % images + prf
    figure('Position',[0 0 500 500])
    imagesc(double(origIm.stimuli(:,:,im_nrs(kk))),[0 max(origIm.stimuli(:))]),hold on
    % Move pRF in 135x135 to original size of 800x800 
    %   half of the screen is 60 in filtered and 400 in original
    orig_x = 400 + (seed_params(1)-res/2) * 400./60;
    orig_y = 400 + (seed_params(2)-res/2) * 400./60;
    orig_sigma = seed_params(3) * 400./60;
    axis image
    colormap gray
    
    numPoints = 50;
    c.th = linspace(0,2*pi, numPoints);
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*orig_sigma);
    plot(c.x + orig_y, c.y + orig_x, 'y','LineWidth',2) % this is just reversed because plot and imagesc are opposite, checked this with contour
    [c.x, c.y] = pol2cart(c.th, 2*ones(1,numPoints)*orig_sigma);
    plot(c.x + orig_y, c.y + orig_x, 'y:','LineWidth',2) % this is just reversed because plot and imagesc are opposite, checked this with contour
    
    axis off
    
    set(gcf,'PaperPositionMode','auto')
    print('-depsc','-r300',fullfile(dataDir,'derivatives','figures',...
            ['Figure8c_sub-' int2str(subj) 'el' int2str(elec) '_im' int2str(im_nrs(kk))]))
    print('-dpng','-r300',fullfile(dataDir,'derivatives','figures',...
            ['Figure8c_sub-' int2str(subj) 'el' int2str(elec) '_im' int2str(im_nrs(kk))]))
end

%% Display NBF model results in 2 stimuli for several electrodes/subjects
% Figure 8A-B

%%%%% Pick a subject:
subject_ind = [19 24];
electrodes = [115 46];

elec_colors = {'c','g'};

figure('Position',[0 0 600 200])
    
for ll = 1:length(electrodes)
    
    subj = subject_ind(ll);
    elec = electrodes(ll);
    
    analysisType = 'spectra200';
    modelType = 'OVmodel';

    % Load NBF model fit:
    load(fullfile(dataDir,'derivatives','gaborFilt','OV_gamma',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
        'seed_params','cross_OVparams','cross_OVestimate','ov_exponents','train_OVperformance')
    ov_exp_used = find(ov_exponents==.5);

    % Load ECoG data:
    dataFitName = fullfile(dataDir,'derivatives','preprocessing',...
        ['sub-' int2str(subj) '_ses-01_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)

    % Calculate ECoG gamma power percent signal change:
    bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
    ecog_bb = mean(100*(10.^(resamp_parms(:,:,2)-bb_base)-1),2);
    ecog_bb_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,2),.16,2)) ...
        squeeze(quantile(resamp_parms(:,:,2),.84,2))]'-bb_base)-1);
    ecog_g = mean(100*(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1),2);
    ecog_g_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.16,2)) ...
        squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.84,2))]')-1);
    ylims = [min(ecog_g_err(:)) max([ecog_g_err(:)])];

    res = sqrt(size(imEnergyMean,2));  % resolution of the pre-processed stimuli
    [~,xx,yy] = makegaussian2d(res,2,2,2,2);

    % Define a Gaussian centered at the prf:
    gaufun1 = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),...
        pp(3),xx,yy,0,0)/(2*pi*pp(3)^2));

    im_nrs = [10 50];

    % Plot gamma responses for these two images
    subplot(length(electrodes),6,6*ll-5),hold on
    bar([1:length(im_nrs)],ecog_g(im_nrs),1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0]);
    plot([1:length(im_nrs); 1:length(im_nrs)],ecog_g_err(:,im_nrs),'k');
    plot(1:length(im_nrs),cross_OVestimate(im_nrs,ov_exp_used)' ,'r.','MarkerSize',20)
    xlim([0 length(im_nrs)+1]),ylim([0 max(max(ecog_g_err(:,im_nrs)))])
    set(gca,'XTick',1:length(im_nrs),'XTickLabel',im_nrs)
    
    %%%%% Plot filtered (contrast) images with prf
    for kk = 1:length(im_nrs)
        % images + prf
        subplot(length(electrodes),6,6*ll-4+kk)
        imagesc(double(origIm.stimuli(:,:,im_nrs(kk))),[0 max(origIm.stimuli(:))]),hold on
        % Move pRF in 135x135 to original size of 800x800 
        %   half of the screen is 60 in filtered and 400 in original
        orig_x = 400 + (seed_params(1)-res/2) * 400./60;
        orig_y = 400 + (seed_params(2)-res/2) * 400./60;
        orig_sigma = seed_params(3) * 400./60;
        axis image
        numPoints = 50;
        c.th = linspace(0,2*pi, numPoints);
        % 1 sd
        [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*orig_sigma);
        plot(c.x + orig_y, c.y + orig_x, 'Color',elec_colors{ll},'LineWidth',2) % this is just reversed because plot and imagesc are opposite, checked this with contour
        % 2 sd
        [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*2*orig_sigma);
        plot(c.x + orig_y, c.y + orig_x, ':','Color',elec_colors{ll},'LineWidth',2) % this is just reversed because plot and imagesc are opposite, checked this with contour
        axis off
        ylabel(['im ' int2str(im_nrs(kk))])

        % zoom
        subplot(length(electrodes),6,6*ll-2+kk)
        imagesc(double(origIm.stimuli(:,:,im_nrs(kk))),[0 max(origIm.stimuli(:))]),hold on
        % Move pRF in 135x135 to original size of 800x800 
        %   half of the screen is 60 in filtered and 400 in original
        orig_x = 400 + (seed_params(1)-res/2) * 400./60;
        orig_y = 400 + (seed_params(2)-res/2) * 400./60;
        orig_sigma = seed_params(3) * 400./60;
        axis image
        numPoints = 50;
        c.th = linspace(0,2*pi, numPoints);
        % 1 sd
        [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*orig_sigma);
        plot(c.x + orig_y, c.y + orig_x,'Color',elec_colors{ll},'LineWidth',2) % this is just reversed because plot and imagesc are opposite, checked this with contour
        % 2 sd
        [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*2*orig_sigma);
        plot(c.x + orig_y, c.y + orig_x,':','Color',elec_colors{ll},'LineWidth',2) % this is just reversed because plot and imagesc are opposite, checked this with contour
        xlim([orig_y-3*orig_sigma orig_y+3*orig_sigma])
        ylim([orig_x-3*orig_sigma orig_x+3*orig_sigma])
        axis off
    end

end
colormap gray
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',fullfile(dataDir,'derivatives','figures',...
        ['Figure8a_sub1-2']))
print('-dpng','-r300',fullfile(dataDir,'derivatives','figures',...
        ['Figure8a_sub1-2']))
