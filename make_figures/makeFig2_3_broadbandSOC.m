% This script will make pannels for Figures 2 and 3 and supplements to
% Figure 3
%
% Hermes D, Petridou N, Kay K, Winawer J. 2019 An image-computable model
% for the stimulus selectivity of gamma oscillations. eLife 2019;8:e47035.
% DOI: https://doi.org/10.7554/eLife.47035
%
%
% The SOC model was adapted from:
% Kay KN, Winawer J, Rokem A, Mezer A, Wandell BA (2013) A Two-Stage
% Cascade Model of BOLD Responses in Human Visual Cortex. PLoS Comput Biol
% 9(5): e1003079. https://doi.org/10.1371/journal.pcbi.1003079
%
% Dora Hermes, 2017

clear all

% set paths:
rootPath = gammaModelPath();
dataDir = fullfile(rootPath,'data');

%% Load preprocessed images and divisive normalization:

load(fullfile(dataDir,'derivatives','gaborFilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

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
imEnergyMean = blob(stimulus,2,8);

%%
%% Get data for figure 2-3: average across subjects (Fig2) + example electrodes (Fig3).

% Get average broadband across all electrodes/subjects:
%%%%% Subjects/electrodes:
subject_ind = [19 19 19 19 19 19  ... % S1
    24 24 ... % S2
    1001 1001 1001 1001 1001 1001 1001 1001]; % S3
electrodes = [107 108 109 115 120 121 ... % S1
    45 46 ... % S2
    49 50 52 57 58 59 60]; % S3

% initialize parameters across subjects
socParams_all = zeros(length(electrodes),6);
socCOD_all = zeros(length(electrodes),2);
SOCestimate_all = zeros(length(electrodes),86);
ecog_bb_all = zeros(length(electrodes),86);
ecog_bb_err_all = zeros(length(electrodes),2,86);
ecog_bb_base = zeros(length(electrodes),1);
ecog_bb_err_base = zeros(length(electrodes),2,1);
pboot = zeros(length(electrodes),86);
ci_base = zeros(length(electrodes),1);

% Get average broadband across all electrodes/subjects:
for ll = 1:length(electrodes)

    subj = subject_ind(ll);
    elec = electrodes(ll);
    
    analysisType = 'spectra200';
    modelType = 'fitSOCbbpower';

    res = sqrt(size(imEnergyMean,2));

    % load model fit
    load(fullfile(dataDir,'derivatives','gaborFilt','SOC_broadband',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
        'seeds','cross_SOCparams','cross_SOCestimate')

    % load ecog stimulus data (1000 bootstraps):
    dataFitName = fullfile(dataDir,'derivatives','preprocessing',...
        ['sub-' int2str(subj) '_ses-01_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)
    % Load ecog baseline data (1000 bootstraps):
    dataBaseFitName = fullfile(dataDir,'derivatives','preprocessing',...
        ['sub-' int2str(subj) '_ses-01_task-soc_allruns_' analysisType '_fitBaseEl' int2str(elec) '.mat']);
    load(dataBaseFitName)

    % get ecog power percent signal change
    bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
    ecog_bb = mean(100*(10.^(resamp_parms(:,:,2)-bb_base)-1),2);
    ecog_bb_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,2),.16,2)) ...
        squeeze(quantile(resamp_parms(:,:,2),.84,2))]'-bb_base)-1);   
    ecog_bb_all(ll,:) = ecog_bb;
    ecog_bb_err_all(ll,:,:) = ecog_bb_err;
        
    ecog_bb_base(ll) = mean(100*(10.^(resamp_parms_baseline(:,:,2)-bb_base)-1),2);
    ecog_bb_err_base(ll,:) = 100*(10.^([squeeze(quantile(resamp_parms_baseline(:,:,2),.16,2)) ...
        squeeze(quantile(resamp_parms_baseline(:,:,2),.84,2))]'-bb_base)-1);   
    
    % run bootstrap test:
    ecog_bb_100 = 100*(10.^(resamp_parms(:,:,2)-bb_base)-1);
    ecog_bb_base_100 = 100*(10.^(resamp_parms_baseline(:,:,2)-bb_base)-1);
    [zz,upbase_ci] = ecog_bootstrapStat(ecog_bb_100,ecog_bb_base_100,'bootdiff');
    pboot(ll,:) = zz; clear zz
    ci_base(ll) = upbase_ci; clear upbase_ci
    
    % get mean model parameters and plot prediction
    cross_SOCparams(cross_SOCparams(:,6)<0,6) = 0; % restrictrange at 0
    cross_SOCparams(cross_SOCparams(:,6)>1,6) = 1; % restrictrange at 1
    cross_SOCparams(:,3) = abs(cross_SOCparams(:,3)); % size>0
    % write out median model parameters across 86 leave 1 outs
    socParams_all(ll,:) = median(cross_SOCparams);

    % write out coefficient of determination
    socCOD_all(ll,1) = calccod(cross_SOCestimate,ecog_bb,[],0,1); % subtract mean 
    socCOD_all(ll,2) = calccod(cross_SOCestimate,ecog_bb,[],0,0); % predict mean + variance

    % write out leave 1 out estimate across electrodes 
    SOCestimate_all(ll,:) = cross_SOCestimate;
end

disp(['mean OV model performance: COD = ' int2str(mean(socCOD_all(:,2)))])

%%
%% Figure 2E mean broadband across electrodes

ylims = [min(ecog_bb_err(:)) max(ecog_bb_err(:))];

% Plot the mean broadband power across all electrodes (Figure 3a)
figure('Position',[0 0 470 600])
subplot(8,1,1),hold on
bar(mean(ecog_bb_all,1),1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0],'LineWidth',1);
bb_group_err = std(ecog_bb_all,1)./sqrt(size(ecog_bb_all,1));
bb_group_up = mean(ecog_bb_all,1)+bb_group_err;
bb_group_low = mean(ecog_bb_all,1)-bb_group_err;
plot([1:86; 1:86],[bb_group_up; bb_group_low],'k');
% plot(mean(SOCestimate_all,1)','r','LineWidth',2)

% plot stimulus cutoffs
stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
for k = 1:length(stim_change)
    plot([stim_change(k) stim_change(k)],ylims(1,:),'Color',[.5 .5 .5],'LineWidth',2)
end
set(gca,'XTick',[])
xlim([0 87])
ylabel(['average'])
ylim([min(bb_group_low)-10 max(bb_group_up)+10])

% save Figure 2E
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['Figure2E_' modelType]))
print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['Figure2E_' modelType]))     

    
%% Plot example electrodes (Figure 3)
example_els = [3 5 7 8 9 14];

figure('Position',[0 0 470 600])
for ll = 1:length(example_els)
    elec = electrodes(example_els(ll));
    
    %%% PLOT BROADBAND POWER AND SOC FIT
    subplot(8,1,ll),hold on
    bar(ecog_bb_all(example_els(ll),:),1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0]);
    plot([1:86; 1:86],squeeze(ecog_bb_err_all(example_els(ll),:,:)),'k');
    plot(SOCestimate_all(example_els(ll),:)','r','LineWidth',2)
        
    % plot bootstrap p<0.05
    if ~isempty(find(pboot(example_els(ll),:)<0.05,1))
        plot(find(pboot(example_els(ll),:)<0.05),0,'b.','MarkerSize',10)
    end
    
    % set ylim
    ylims = [min(min(ecog_bb_err_all(example_els(ll),:,:))) max(max(ecog_bb_err_all(example_els(ll),:,:)))];
    ylim(ylims(1,:))
    % plot stimulus cutoffs
    stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
    for k = 1:length(stim_change)
        plot([stim_change(k) stim_change(k)],ylims(1,:),'Color',[.5 .5 .5],'LineWidth',2)
    end
    set(gca,'XTick',[])
    xlim([0 87])
    ylabel(['el' int2str(elec) ' R^2=' int2str(socCOD_all(example_els(ll),2))])

end

% save Figure 3
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['Figure3_' modelType]))
print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['Figure3_' modelType]))


%% Figure 3 supplement: all electrodes all subjects

%%%%% Subjects/electrodes:
subject_ind = [19 19 19 19 19 19  ... % S1
    24 24 ... % S2
    1001 1001 1001 1001 1001 1001 1001 1001]; % S3
electrodes = [107 108 109 115 120 121 ... % S1
    45 46 ... % S2
    49 50 52 57 58 59 60]; % S3

socParams_all = zeros(length(electrodes),6);
socCOD_all = zeros(length(electrodes),2);

plot_nr = 0; 
figure_nr = 1;
figure('Position',[0 0 470 600])
for ll = 1:length(electrodes)
    clear resamp_parms_baseline resamp_parms
    plot_nr = plot_nr + 1;

    subj = subject_ind(ll);
    elec = electrodes(ll);
    
    analysisType = 'spectra200';
    modelType = 'fitSOCbbpower';

    res = sqrt(size(imEnergyMean,2));
    
    % load model fit
    load(fullfile(dataDir,'derivatives','gaborFilt','SOC_broadband',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
        'seeds','cross_SOCparams','cross_SOCestimate')

    % load ecog data:
    dataFitName = fullfile(dataDir,'derivatives','preprocessing',...
        ['sub-' int2str(subj) '_ses-01_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)

    % Load ecog baseline data:
    dataBaseFitName = fullfile(dataDir,'derivatives','preprocessing',...
        ['sub-' int2str(subj) '_ses-01_task-soc_allruns_' analysisType '_fitBaseEl' int2str(elec) '.mat']);
    load(dataBaseFitName)
    
    % get ecog power percent signal change
    bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
    ecog_bb = mean(100*(10.^(resamp_parms(:,:,2)-bb_base)-1),2);
    ecog_bb_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,2),.16,2)) ...
        squeeze(quantile(resamp_parms(:,:,2),.84,2))]'-bb_base)-1);
    
    ecog_bb_base = mean(100*(10.^(resamp_parms_baseline(:,:,2)-bb_base)-1),2);
    ecog_bb_err_base = 100*(10.^([squeeze(quantile(resamp_parms_baseline(:,:,2),.16,2)) ...
        squeeze(quantile(resamp_parms_baseline(:,:,2),.84,2))]'-bb_base)-1);   

    
    %%% PLOT BROADBAND POWER AND SOC FIT
    subplot(8,1,plot_nr),hold on
    bar(ecog_bb,1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0]);
    plot([1:86; 1:86],ecog_bb_err,'k');
    plot(cross_SOCestimate' ,'r','LineWidth',2)
    
    % plot bootstrap p<0.05
    if ~isempty(find(pboot(ll,:)<0.05,1))
        plot(find(pboot(ll,:)<0.05),0,'b.','MarkerSize',10)
    end
    
    % set ylim
    ylims = [min(ecog_bb_err(:)) max(ecog_bb_err(:))];
    ylim(ylims(1,:))
    % plot stimulus cutoffs
    stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
    for k = 1:length(stim_change)
        plot([stim_change(k) stim_change(k)],ylims(1,:),'Color',[.5 .5 .5],'LineWidth',2)
    end
    set(gca,'XTick',[])
    xlim([0 87])

    % get mean model parameters and plot prediction
    cross_SOCparams(cross_SOCparams(:,6)<0,6) = 0; % restrictrange at 0
    cross_SOCparams(cross_SOCparams(:,6)>1,6) = 1; % restrictrange at 1
    cross_SOCparams(:,3) = abs(cross_SOCparams(:,3)); % size>0
    socParams_all(ll,:) = median(cross_SOCparams);

    socCOD_all(ll,1) = calccod(cross_SOCestimate,ecog_bb,[],0,1);
    socCOD_all(ll,2) = calccod(cross_SOCestimate,ecog_bb,[],0,0);
    
    ylabel(['el' int2str(elec) ' R^2=' int2str(socCOD_all(ll,2))])

    if mod(ll,8)==0 && ll<length(electrodes)% save figure and make a new one every 8 electrodes
        % save the figure
        set(gcf,'PaperPositionMode','auto')
        print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
                ['Figure3_supplement_elset' int2str(figure_nr) '_' modelType]))
        print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
                ['Figure3_supplement_elset' int2str(figure_nr) '_' modelType]))

        % and make a new figure
        figure_nr = figure_nr +1;
        figure('Position',[0 0 470 600])
        % reset the subplot number
        plot_nr = 0;
        
    elseif ll==length(electrodes)% save figure after last electrode
        % save the last figure
        set(gcf,'PaperPositionMode','auto')
        print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
                ['Figure3_supplement_elset' int2str(figure_nr) '_' modelType]))
        print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
                ['Figure3_supplement_elset' int2str(figure_nr) '_' modelType]))
    end
end


