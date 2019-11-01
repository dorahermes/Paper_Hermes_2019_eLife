% This script will make pannels for Figures 2 and 6 and supplements to
% Figure 6
%
% Hermes D, Petridou N, Kay K, Winawer J. 2019 An image-computable model
% for the stimulus selectivity of gamma oscillations. eLife 2019;8:e47035.
% DOI: https://doi.org/10.7554/eLife.47035
%
%
% Dora Hermes, 2017

clear all

% set paths:
rootPath = gammaModelPath();
dataDir = fullfile(rootPath,'data');

%% Load preprocessed images:

load(fullfile(dataDir,'derivatives','gaborFilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

% compute the square root of the sum of the squares of the outputs of
% quadrature-phase filter pairs (this is the standard complex-cell energy model).
% after this step, stimulus is images x orientations*positions.
stimulus = sqrt(blob(stimulus.^2,2,2));
% the square root may not be typical for complex cell energy model

% we filtered the image with 8 orientations:
nrOrientations = 8;
res = sqrt(size(stimulus,2)/nrOrientations);  % resolution of the pre-processed stimuli

%%
%% Get data for figure 2-6: average across subjects (Fig2) + example electrodes (Fig6).
% Figure 2: gamma average across subjects 
% Figure 6: gamma/OV model in 6 example electrodes 

%%%%% Subjects/electrodes:
subject_ind = [19 19 19 19 19 19  ... % S1
    24 24 ... % S2
    1001 1001 1001 1001 1001 1001 1001 1001]; % S3
electrodes = [107 108 109 115 120 121 ... % S1
    45 46 ... % S2
    49 50 52 57 58 59 60]; % S3

% initialize parameters across subjects
OV_COD_all = zeros(length(electrodes),1);
OVparams_all = zeros(length(electrodes),5);
OVestimate_all = zeros(length(electrodes),86);
ecog_g_all = zeros(length(electrodes),86);
ecog_g_err_all = zeros(length(electrodes),2,86);
ecog_g_base = zeros(length(electrodes),1);
ecog_g_err_base = zeros(length(electrodes),2,1);
pboot = zeros(length(electrodes),86);
ci_base = zeros(length(electrodes),1);

for ll = 1:length(electrodes)
    
    subj = subject_ind(ll);
    elec = electrodes(ll);

    analysisType = 'spectra200';
    modelType = 'OVmodel';

    res = sqrt(size(stimulus,2)/8);

    % Load model fit:
    load(fullfile(dataDir,'derivatives','gaborFilt','OV_gamma',...
            ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
            'seed_params','cross_OVparams','cross_OVestimate','ov_exponents','train_OVperformance')
    % Median model parameters across leav-one-out (actually, only the gain)
    OVparams_all(ll,:) = median(cross_OVparams(:,:,ov_exponents==.5),1);
    OVestimate_all(ll,:) = cross_OVestimate(:,:,ov_exponents==.5);
    
    % Load ecog stimulus data:
    dataFitName = fullfile(dataDir,'derivatives','preprocessing',...
        ['sub-' int2str(subj) '_ses-01_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)
    % Load ecog baseline data:
    dataBaseFitName = fullfile(dataDir,'derivatives','preprocessing',...
        ['sub-' int2str(subj) '_ses-01_task-soc_allruns_' analysisType '_fitBaseEl' int2str(elec) '.mat']);
    load(dataBaseFitName)

    % Gamma power percent signal change:
    ecog_g = 100*(mean(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1,2));
    ecog_g_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.16,2)) ...
        squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.84,2))]')-1);
    ecog_g_all(ll,:) = ecog_g;
    ecog_g_err_all(ll,:,:) = ecog_g_err;
    
    ecog_g_base(ll) = 100*(mean(10.^(resamp_parms_baseline(:,:,3)./resamp_parms_baseline(:,:,5))-1,2));
    ecog_g_err_base(ll,:) = 100*(10.^([squeeze(quantile(resamp_parms_baseline(:,:,3)./resamp_parms_baseline(:,:,5),.16,2)) ...
        squeeze(quantile(resamp_parms_baseline(:,:,3)./resamp_parms_baseline(:,:,5),.84,2))]')-1);

    % run bootstrap test:
    ecog_g_100 = 100*(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1);
    ecog_g_base_100 = 100*(10.^(resamp_parms_baseline(:,:,3)./resamp_parms_baseline(:,:,5))-1);
    [zz,upbase_ci] = ecog_bootstrapStat(ecog_g_100,ecog_g_base_100,'bootdiff');

    pboot(ll,:) = zz; clear zz
    ci_base(ll) = upbase_ci; clear upbase_ci
    
    % model performance (no mean subtracted)
    r2 = calccod(squeeze(cross_OVestimate(:,:,ov_exponents==.5)),ecog_g,[],0,0);
    OV_COD_all(ll) = r2;
end 

disp(['mean OV model performance: COD = ' int2str(mean(OV_COD_all))])

%% Plot the mean gamma power across all electrodes (Figure 2E)

figure('Position',[0 0 600 600])
subplot(8,5,2:5),hold on

% plot stimulus data and error:
bar(mean(ecog_g_all,1),1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0],'LineWidth',1);
g_group_err = std(ecog_g_all,1)./sqrt(size(ecog_g_all,1));
g_group_up = mean(ecog_g_all,1)+g_group_err;
g_group_low = mean(ecog_g_all,1)-g_group_err;
plot([1:86; 1:86],[g_group_up; g_group_low],'k');

% plot baseline gamma and sterr
g_base_mean = mean(ecog_g_base,1);
g_base_err = std(ecog_g_base,1)./sqrt(size(ecog_g_base,1));
g_base_up = mean(ecog_g_base,1)+g_base_err;
g_base_low = mean(ecog_g_base,1)-g_base_err;

ylims = [min(g_group_low(:)) max(g_group_up(:))];

% plot stimulus cutoffs
stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
for k = 1:length(stim_change)
    plot([stim_change(k) stim_change(k)],ylims(1,:),'Color',[.5 .5 .5],'LineWidth',2)
end

set(gca,'XTick',[])
xlim([0 87])
ylabel(['average'])
ylim([min(g_group_low)-10 max(g_group_up)+10])
% save Figure 2E
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['Figure2E_' modelType]))
print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['Figure2E_' modelType]))


%%
%%%%%% plot example electrodes and average across electrodes
% example electrodes
example_els = [3 5 7 8 9 14];

figure('Position',[0 0 600 600])
for ll = 1:length(example_els)
    elec = electrodes(example_els(ll));
    
    %%% PLOT GAUSSIAN
    subplot(8,5,5*ll-4),% NO HOLD ON HERE, only after imagesc, otherwise Y axis inverted
    res = sqrt(size(stimulus,2)/8);
    [~,xx,yy] = makegaussian2d(res,2,2,2,2);
    imagesc(ones(size(xx)),[0 1]);
    axis image, hold on, colormap gray
    plot([res/2 res/2],[1 res],'k'),plot([1 res],[res/2 res/2],'k')
    %%% Where is the seed prf:
    numPoints = 50;
    c.th = linspace(0,2*pi, numPoints);
    % plot 1 sd
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*OVparams_all(example_els(ll),3));
    plot(c.x + OVparams_all(example_els(ll),2), c.y + OVparams_all(example_els(ll),1), 'r') % this is just reversed because plot and imagesc are opposite, checked this with contour
    % plot 2 sd
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*2*OVparams_all(example_els(ll),3));
    plot(c.x + OVparams_all(example_els(ll),2), c.y + OVparams_all(example_els(ll),1), 'r:') % this is just reversed because plot and imagesc are opposite, checked this with contour
    axis off
    
    %%% PLOT GAMMA POWER AND OV FIT
    subplot(8,5,5*ll-3:5*ll),hold on
    bar(ecog_g_all(example_els(ll),:),1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0]);
    plot([1:86; 1:86],squeeze(ecog_g_err_all(example_els(ll),:,:)),'k');
    plot(OVestimate_all(example_els(ll),:)','r','LineWidth',2)
        
    % plot bootstrap diff test
    if ~isempty(find(pboot(example_els(ll),:)<.05,1))
        plot(find(pboot(example_els(ll),:)<.05),0,'b.','MarkerSize',10)
    end

    % set ylim
    ylims = [min(min(ecog_g_err_all(example_els(ll),:,:))) max(max(ecog_g_err_all(example_els(ll),:,:)))];
    ylim(ylims(1,:))
    % plot stimulus cutoffs
    stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
    for k = 1:length(stim_change)
        plot([stim_change(k) stim_change(k)],ylims(1,:),'Color',[.5 .5 .5],'LineWidth',2)
    end
    set(gca,'XTick',[])
    xlim([0 87])
    ylabel(['el' int2str(elec) ' R^2=' int2str(OV_COD_all(example_els(ll),1))])
end

% save Figure 6
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['Figure6_' modelType]))
print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['Figure6_' modelType]))

%%
%% Figure 6 supplement: OV model results for all electrodes, all subjects

%%%%% Subjects/electrodes:
subject_ind = [19 19 19 19 19 19  ... % S1
    24 24 ... % S2
    1001 1001 1001 1001 1001 1001 1001 1001]; % S3
electrodes = [107 108 109 115 120 121 ... % S1
    45 46 ... % S2
    49 50 52 57 58 59 60]; % S3

plot_nr = 0; 
figure_nr = 1;
figure('Position',[0 0 470 600])    
for ll = 1:length(electrodes)
    plot_nr = plot_nr + 1;
    
    subj = subject_ind(ll);
    elec = electrodes(ll);

    analysisType = 'spectra200';
    modelType = 'OVmodel';

    res = sqrt(size(stimulus,2)/8); % filtered stimulus resolution, assuming 8 orientations

    %%%%% Load model fit:
    load(fullfile(dataDir,'derivatives','gaborFilt','OV_gamma',...
            ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
            'seed_params','cross_OVparams','cross_OVestimate','ov_exponents','train_OVperformance')

    %%%%% Load ecog data:
    dataFitName = fullfile(dataDir,'derivatives','preprocessing',...
        ['sub-' int2str(subj) '_ses-01_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)

    %%%%% Gamma power percent signal change:
    ecog_g = 100*(mean(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1,2));
    ecog_g_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.16,2)) ...
        squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.84,2))]')-1);
    ecog_g_base = 100*(median(10.^(resamp_parms_baseline(:,:,3)./resamp_parms_baseline(:,:,5))-1,2));
    ecog_g_err_base = 100*(10.^([squeeze(quantile(resamp_parms_baseline(:,:,3)./resamp_parms_baseline(:,:,5),.16,2)) ...
        squeeze(quantile(resamp_parms_baseline(:,:,3)./resamp_parms_baseline(:,:,5),.84,2))]')-1);

    % Calculate ylims for gamma range:
    ylims = [min(ecog_g_err(:)) max([ecog_g_err(:)])];
    
    %%%%% Make plot for this electrode:
    subplot(8,1,plot_nr),hold on   
    bar(ecog_g,1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0]);
    plot([1:86; 1:86],ecog_g_err,'k');
    % Add gamma estimate (cross validated):
    plot(squeeze(cross_OVestimate(:,:,ov_exponents==.5)),'r','LineWidth',2)
    % Plot stimulus cutoffs:
    stim_change = [38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
    for k = 1:length(stim_change)
        plot([stim_change(k) stim_change(k)],ylims(1,:),'Color',[.5 .5 .5],'LineWidth',2)
    end
    xlim([0 87]), ylim(ylims(1,:))
    set(gca,'XTick',[])
    ylabel('gamma')
    
    % plot bootstrap diff test
    if ~isempty(find(pboot(ll,:)<.05,1))
        plot(find(pboot(ll,:)<.05),0,'b.','MarkerSize',10)
    end

    %%%%% Add model performance in label
    r2 = calccod(squeeze(cross_OVestimate(:,:,ov_exponents==.5)),ecog_g,[],0,0);
    ylabel(['el' int2str(elec) ' R^2=' int2str(r2)])
    
    if mod(ll,8)==0 && ll<length(electrodes)% save figure and make a new one every 8 electrodes
        % Save Figure S8
        set(gcf,'PaperPositionMode','auto')
        print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
                ['Figure6_Supplement_elset' int2str(figure_nr) '_' modelType]))
        print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
                ['Figure6_Supplement_elset' int2str(figure_nr) '_' modelType]))

        % Make a new figure
        figure_nr = figure_nr +1;
        figure('Position',[0 0 470 600])
        % reset the subplot number
        plot_nr = 0;
%         
    elseif ll==length(electrodes)% save figure after last electrode
        % Save Figure S9
        set(gcf,'PaperPositionMode','auto')
        print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
                ['Figure6_Supplement_elset' int2str(figure_nr) '_' modelType]))
        print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
                ['Figure6_Supplement_elset' int2str(figure_nr) '_' modelType]))
    end
    
end

%%

%% Look at where the Gaussian is for all electrodes

%%%%% Subjects/electrodes:
subject_ind = [19 19 19 19 19 19  ... % S1
    24 24 ... % S2
    1001 1001 1001 1001 1001 1001 1001 1001]; % S3
electrodes = [107 108 109 115 120 121 ... % S1
    45 46 ... % S2
    49 50 52 57 58 59 60]; % S3

plot_nr = 0; 
figure_nr = 1;
figure('Position',[0 0 200 600])    
for ll = 1:length(electrodes)
    plot_nr = plot_nr + 1;
    
    subj = subject_ind(ll);
    elec = electrodes(ll);

    analysisType = 'spectra200';
    modelType = 'OVmodel';
    res = sqrt(size(stimulus,2)/8);

     %%%%% Load model fit:
    load(fullfile(dataDir,'derivatives','gaborFilt','OV_gamma',...
            ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
            'seed_params','cross_OVparams','cross_OVestimate','ov_exponents','train_OVperformance')

    %%% LOOK AT WHERE THE GAUSSIAN IS
    subplot(8,1,plot_nr) % NO HOLD ON HERE, only after imagesc, otherwise Y axis inverted
    [~,xx,yy] = makegaussian2d(res,2,2,2,2);
    imagesc(ones(size(xx)),[0 1]);
    axis image, hold on, colormap gray
    plot([res/2 res/2],[1 res],'k'),plot([1 res],[res/2 res/2],'k')

    %%% Where is the seed prf:
    numPoints = 50;
    c.th = linspace(0,2*pi, numPoints);
    % plot 1 sd
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*seed_params(3));
    plot(c.x + seed_params(2), c.y + seed_params(1), 'r') % this is just reversed because plot and imagesc are opposite, checked this with contour
    % plot 2 sd
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*2*seed_params(3));
    plot(c.x + seed_params(2), c.y + seed_params(1), 'r:') % this is just reversed because plot and imagesc are opposite, checked this with contour
    axis off
    
    if mod(ll,8)==0 && ll<length(electrodes)% save figure and make a new one every 8 electrodes
        % Save Figure S8a
        set(gcf,'PaperPositionMode','auto')
        print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
                ['FigureS8a_elset' int2str(figure_nr) '_' modelType]))
        print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
                ['FigureS8a_elset' int2str(figure_nr) '_' modelType]))

        % and make a new figure
        figure_nr = figure_nr +1;
        figure('Position',[0 0 200 600])
        % reset the subplot number
        plot_nr = 0;
%         
    elseif ll==length(electrodes)% save figure after last electrode
        % Save Figure S9a
        set(gcf,'PaperPositionMode','auto')
        print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
                ['FigureS9a_elset' int2str(figure_nr) '_' modelType]))
        print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
                ['FigureS9a_elset' int2str(figure_nr) '_' modelType]))
    end
end

%%
%% reviewer comment about visualizing errors:

% all electrodes in 1
figure('Position',[0 0 300 300]),hold on
for ll = 1:length(electrodes)
    plot([0:max(ecog_g_all(:))],[0:max(ecog_g_all(:))],'k','LineWidth',1)
    plot(OVestimate_all',ecog_g_all','o')
end
xlabel('OV prediction')
ylabel('gamma power')
axis square
axis tight 
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['FigureReview_' modelType]))
print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['FigureReview_' modelType]))
        
% all electrodes in 1
figure('Position',[0 0 300 300]),hold on
for ll = 1:length(electrodes)
    plot([0:max(ecog_g_all(:))],[0:max(ecog_g_all(:))],'k','LineWidth',1)
    plot(OVestimate_all',ecog_g_all','o')
end
xlabel('OV prediction')
ylabel('gamma power')
axis square
xlim([0 500]),ylim([0 500])

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['FigureReview_' modelType '_zoom']))
print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['FigureReview_' modelType '_zoom']))


