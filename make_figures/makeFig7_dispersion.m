% This script will make a figure of the R^2 of the OV and SOC model
% predictions of the gamma and brodband signals. This is provided in the
%
% Hermes D, Petridou N, Kay K, Winawer J. 2019 An image-computable model
% for the stimulus selectivity of gamma oscillations. eLife 2019;8:e47035.
% DOI: https://doi.org/10.7554/eLife.47035
%
% The SOC model was adapted from:
% Kay KN, Winawer J, Rokem A, Mezer A, Wandell BA (2013) A Two-Stage
% Cascade Model of BOLD Responses in Human Visual Cortex. PLoS Comput Biol
% 9(5): e1003079. https://doi.org/10.1371/journal.pcbi.1003079
%
% Dora Hermes, 2019

clear all

% set paths:
rootPath = gammaModelPath();
dataDir = fullfile(rootPath,'data');

%%
%% Get ECoG broadband and gamma data

%%%%% Subjects/electrodes:
subject_ind = [19 19 19 19 19 19  ... % S1
    24 24 ... % S2
    1001 1001 1001 1001 1001 1001 1001 1001]; % S3
electrodes = [107 108 109 115 120 121 ... % S1
    45 46 ... % S2
    49 50 52 57 58 59 60]; % S3

% initialize parameters across subjects
SOCbb_COD_all = zeros(length(electrodes),1);
SOCbbestimate_all = zeros(length(electrodes),86);

SOCg_COD_all = zeros(length(electrodes),1);
SOCgestimate_all = zeros(length(electrodes),86);

OVbb_COD_all = zeros(length(electrodes),1);
OVbbestimate_all = zeros(length(electrodes),86);

OVg_COD_all = zeros(length(electrodes),1);
OVgestimate_all = zeros(length(electrodes),86);

MEANbb_COD_all = zeros(length(electrodes),1);
MEANbbestimate_all = zeros(length(electrodes),86);

MEANg_COD_all = zeros(length(electrodes),1);
MEANgestimate_all = zeros(length(electrodes),86);

ecog_bb_all = zeros(length(electrodes),86);
ecog_bb_err_all = zeros(length(electrodes),2,86);
ecog_g_all = zeros(length(electrodes),86);
ecog_g_err_all = zeros(length(electrodes),2,86);

for ll = 1:length(electrodes)
    
    subj = subject_ind(ll);
    elec = electrodes(ll);
    
    % 1) get SOC broadband model fits
    analysisType = 'spectra200';
    modelType = 'fitSOCbbpower';
    SOCbb = load(fullfile(dataDir,'derivatives','gaborFilt','SOC_broadband',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
        'seeds','cross_SOCparams','cross_SOCestimate');
    
    % 2) get SOC gamma model fits
    analysisType = 'spectra200';
    modelType = 'fitSOCgammapower';
    SOCgamma = load(fullfile(dataDir,'derivatives','gaborFilt','SOC_gamma',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
        'seeds','cross_SOCparams','cross_SOCestimate');
    
    % 3) get OV broadband model fits
    analysisType = 'spectra200';
    modelType = 'OVmodel';
    OVbb = load(fullfile(dataDir,'derivatives','gaborFilt','OV_broadband',...
            ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
            'seed_params','cross_OVparams','cross_OVestimate','ov_exponents','train_OVperformance');

    % 4) get OV gamma model fits
    analysisType = 'spectra200';
    modelType = 'OVmodel';
    OVgamma = load(fullfile(dataDir,'derivatives','gaborFilt','OV_gamma',...
            ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
            'seed_params','cross_OVparams','cross_OVestimate','ov_exponents','train_OVperformance');

    % 5) get MEAN broadband model fits
    analysisType = 'spectra200';
    modelType = 'MeanBB';
    MEANbb = load(fullfile(dataDir,'derivatives','gaborFilt','mean_broadband',...
            ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
            'cross_MeanPred');

    % 6) get MEAN gamma model fits
    analysisType = 'spectra200';
    modelType = 'MeanGamma';
    MEANgamma = load(fullfile(dataDir,'derivatives','gaborFilt','mean_gamma',...
            ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
            'cross_MeanPred');
                
    % save model estimate outputs:
    SOCbbestimate_all(ll,:) = SOCbb.cross_SOCestimate;
    SOCgestimate_all(ll,:) = SOCgamma.cross_SOCestimate;
    OVbbestimate_all(ll,:) = OVbb.cross_OVestimate(:,:,OVbb.ov_exponents==.5);
    OVgestimate_all(ll,:) = OVgamma.cross_OVestimate(:,:,OVgamma.ov_exponents==.5);
    MEANbbestimate_all(ll,:) = MEANbb.cross_MeanPred;
    MEANgestimate_all(ll,:) = MEANgamma.cross_MeanPred;
    
    % load ecog stimulus data (1000 bootstraps):
    dataFitName = fullfile(dataDir,'derivatives','preprocessing',...
        ['sub-' int2str(subj) '_ses-01_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)

    % Broadband power percent signal change:
    bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
    ecog_bb = mean(100*(10.^(resamp_parms(:,:,2)-bb_base)-1),2);
    ecog_bb_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,2),.16,2)) ...
        squeeze(quantile(resamp_parms(:,:,2),.84,2))]'-bb_base)-1);   
    ecog_bb_all(ll,:) = ecog_bb;
    ecog_bb_err_all(ll,:,:) = ecog_bb_err;
    
    % Gamma power percent signal change:
    ecog_g = 100*(mean(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1,2));
    ecog_g_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.16,2)) ...
        squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.84,2))]')-1);
    ecog_g_all(ll,:) = ecog_g;
    ecog_g_err_all(ll,:,:) = ecog_g_err;

    % save model fit outputs, no means subtracted:
    SOCbb_COD_all(ll) = calccod(SOCbb.cross_SOCestimate,ecog_bb,[],0,0); % predict mean + variance
    SOCg_COD_all(ll) = calccod(SOCgamma.cross_SOCestimate,ecog_g,[],0,0); % predict mean + variance

    OVbb_COD_all(ll) = calccod(OVbb.cross_OVestimate(:,:,OVbb.ov_exponents==.5),ecog_bb,[],0,0);
    OVg_COD_all(ll) = calccod(OVgamma.cross_OVestimate(:,:,OVgamma.ov_exponents==.5),ecog_g,[],0,0);

    MEANbb_COD_all(ll) = calccod(MEANbb.cross_MeanPred,ecog_bb,[],0,0);
    MEANg_COD_all(ll) = calccod(MEANgamma.cross_MeanPred,ecog_g,[],0,0);

end 

mean([SOCbb_COD_all SOCg_COD_all OVbb_COD_all OVg_COD_all MEANbb_COD_all MEANg_COD_all])

% test whether SOC is better for broadband than OV:
[~,p,hi,stats] = ttest(fisherz(SOCbb_COD_all./100),fisherz(OVbb_COD_all./100));

% test whether SOC is better for broadband than mean:
[~,p,hi,stats] = ttest(fisherz(SOCbb_COD_all./100),fisherz(MEANbb_COD_all./100));


% test whether OV is better for gamma than SOC:
[~,p,hi,stats] = ttest(fisherz(SOCg_COD_all./100),fisherz(OVg_COD_all./100));

% test whether OV is better for gamma than mean:
[~,p,hi,stats] = ttest(fisherz(SOCg_COD_all./100),fisherz(MEANg_COD_all./100));


%% plot OV vs SOC models

figure('Position',[0 0 300 150]),hold on

plot(OVbb_COD_all,SOCbb_COD_all,'.','MarkerSize',5,'Color',[0 0 1])
plot(OVg_COD_all,SOCg_COD_all,'.','MarkerSize',5,'Color',[1 0 0])

% plot error bars across the datapoints
nn = length(SOCbb_COD_all);
group_means = mean([SOCbb_COD_all,SOCg_COD_all,OVbb_COD_all,OVg_COD_all]);
group_low = mean([SOCbb_COD_all,SOCg_COD_all,OVbb_COD_all,OVg_COD_all]) - ...
    2*std([SOCbb_COD_all,SOCg_COD_all,OVbb_COD_all,OVg_COD_all])./sqrt(nn);
group_up = mean([SOCbb_COD_all,SOCg_COD_all,OVbb_COD_all,OVg_COD_all]) + ...
    2*std([SOCbb_COD_all,SOCg_COD_all,OVbb_COD_all,OVg_COD_all])./sqrt(nn);

plot([group_means(3) group_means(3)],[group_means(1) group_means(1)],'.','MarkerSize',20,'Color',[0 0 .8])
plot([group_means(3) group_means(3)],[group_low(1) group_up(1)],'Color',[0 0 .8],'Linewidth',2)
plot([group_low(3) group_up(3)],[group_means(1) group_means(1)],'Color',[0 0 .8],'Linewidth',2)

plot([group_means(4) group_means(4)],[group_means(2) group_means(2)],'.','MarkerSize',20,'Color',[.8 0 0])
plot([group_means(4) group_means(4)],[group_low(2) group_up(2)],'Color',[.8 0 0],'Linewidth',2)
plot([group_low(4) group_up(4)],[group_means(2) group_means(2)],'Color',[.8 0 0],'Linewidth',2)

xlim([0 100]),ylim([0 100])
plot([0 100],[0 100],'k')
legend('broadband','gamma','Location','eastoutside')
xlabel('R^2 OV model')
ylabel('R^2 SOC model')
set(gca,'XTick',[0:20:100],'YTick',[0:20:100])
axis square

% Save Dispersion Figure
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['Figure7']))
print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['Figure7']))

    
%% plot Spearman

figure('Position',[0 0 300 150]),hold on

plot(OVbb_Spearman_all.^2,SOCbb_Spearman_all.^2,'b.','MarkerSize',10)
plot(OVg_Spearman_all.^2,SOCg_Spearman_all.^2,'r.','MarkerSize',10)

% plot error bars across the datapoints
nn = length(SOCbb_Spearman_all);
group_means = mean([SOCbb_Spearman_all,SOCg_Spearman_all,OVbb_Spearman_all,OVg_Spearman_all]);
group_low = mean([SOCbb_Spearman_all,SOCg_Spearman_all,OVbb_Spearman_all,OVg_Spearman_all]) - ...
    2*std([SOCbb_Spearman_all,SOCg_Spearman_all,OVbb_Spearman_all,OVg_Spearman_all])./sqrt(nn);
group_up = mean([SOCbb_Spearman_all,SOCg_Spearman_all,OVbb_Spearman_all,OVg_Spearman_all]) + ...
    2*std([SOCbb_Spearman_all,SOCg_Spearman_all,OVbb_Spearman_all,OVg_Spearman_all])./sqrt(nn);

plot([group_means(3) group_means(3)],[group_means(1) group_means(1)],'.','MarkerSize',20,'Color',[0 0 .8])
plot([group_means(3) group_means(3)],[group_low(1) group_up(1)],'Color',[0 0 .8],'Linewidth',2)
plot([group_low(3) group_up(3)],[group_means(1) group_means(1)],'Color',[0 0 .8],'Linewidth',2)

plot([group_means(4) group_means(4)],[group_means(2) group_means(2)],'.','MarkerSize',20,'Color',[.8 0 0])
plot([group_means(4) group_means(4)],[group_low(2) group_up(2)],'Color',[.8 0 0],'Linewidth',2)
plot([group_low(4) group_up(4)],[group_means(2) group_means(2)],'Color',[.8 0 0],'Linewidth',2)

xlim([0 1]),ylim([0 1])
plot([0 1],[0 1],'k')
legend('broadband','gamma','Location','eastoutside')
xlabel('Spearman OV model')
ylabel('Spearman SOC model')
set(gca,'XTick',[0:.2:1],'YTick',[0:.2:1])
axis square
 
% Save Dispersion Figure
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['Dispersion_Spearman']))
print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['Dispersion_Spearman']))

%% plot Scatter prediction & data

figure('Position',[0 0 750 250]),hold on

subplot(1,3,1)
plot(SOCbbestimate_all',ecog_bb_all','.','MarkerSize',10)

subplot(1,3,2)
plot(OVgestimate_all',ecog_g_all','.','MarkerSize',10)

subplot(1,3,3)
plot(repmat([1:15],15,1),repmat([1:15],15,1),'.','MarkerSize',10)
title('electrode numbers/colors')

legend({'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'},'Location','eastoutside')

% % Save Dispersion Figure
% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
%         ['Dispersion_Spearman']))
% print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
%         ['Dispersion_Spearman']))


%% plot histogram of errors, prediction - data

figure('Position',[0 0 750 400])

subplot(2,2,1)
for kk = 1:15
    histogram(SOCbbestimate_all(kk,:)-ecog_bb_all(kk,:),[-1000:20:1000])
    hold on
end
subplot(2,2,3)
histogram(OVbbestimate_all(:)-ecog_bb_all(:),[-1000:20:1000])


subplot(2,2,2)
histogram(SOCgestimate_all(:)-ecog_g_all(:),[-1000:20:1000])
subplot(2,2,4)
histogram(OVgestimate_all(:)-ecog_g_all(:),[-1000:20:1000])

%%

figure('Position',[0 0 750 400])

subplot(2,2,1)
[n,x] = hist(SOCbbestimate_all'-ecog_bb_all',[-500:10:500]);
bar(x,n,10)
title('SOC model error (estimate - broadband)')
xlabel('error')
ylabel('# stimuli')

subplot(2,2,3)
[n,x] = hist(OVbbestimate_all'-ecog_bb_all',[-500:10:500]);
bar(x,n,10)
title('OV model error (estimate - broadband)')
xlabel('error')
ylabel('# stimuli')

subplot(2,2,2)
[n,x] = hist(SOCgestimate_all'-ecog_g_all',[-500:10:500]);
bar(x,n,10)
title('SOC model error (estimate - gamma)')
xlabel('error')
ylabel('# stimuli')

subplot(2,2,4)
[n,x] = hist(OVgestimate_all'-ecog_g_all',[-500:10:500]);
bar(x,n,10)
title('OV model error (estimate - gamma)')
xlabel('error')
ylabel('# stimuli')


% Save Dispersion Figure
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['error_histograms']))
print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['error_histograms']))
