%
%
% This script generate the power spectra for example electrodes as shown in
% Figure 1B from:
%
% Hermes D, Petridou N, Kay K, Winawer J. 2019 An image-computable model
% for the stimulus selectivity of gamma oscillations. eLife 2019;8:e47035.
% DOI: https://doi.org/10.7554/eLife.47035
%
% Dora Hermes 2019

% clear all
clearvars

% set paths:
rootPath = gammaModelPath();
dataDir = fullfile(rootPath,'data');


%%
%% Plot example channel noise and grating
%%

subjects = {'19','24','1001'};
electrodes = {109,45,50}; % example electrodes per subject
analysisType = 'spectra200';

figure('Position',[0 0 500 150]),hold on

for s = 1:length(subjects)
    subj = subjects{s};
    elec = electrodes{s};

    % load all data:
    dataName = fullfile(dataDir,'derivatives','preprocessing',['sub-' subj],'ses-01','ieeg',...
        ['sub-' subj '_ses-01_task-soc_allruns_' analysisType '.mat']);
    load(dataName,'f','spectra','spectra_off','stims','runs','electrodes_incl')   
    elec_nr = find(electrodes_incl == elec);
    
    % get data for this electrode
    data_fft = squeeze(spectra(elec_nr,:,:));
    data_fft_off = squeeze(spectra_off(elec_nr,:,:));
    
    % ranges used for fitting
    f_use4fit = [30:57 65:115 126:175 186:200];
    f_sel = ismember(f,f_use4fit);
    nr_boots = 10;
    
    stims_plot = [45 83];
    stims_color = {[1 .1 .1],[0 .3 .9]};

    data_base_boot = zeros(nr_boots,size(data_fft_off,2));
    % bootstrap data and plot:
    for ii = 1:nr_boots
        % from baseline data, random sample with replacement
        trial_set = randsample(size(data_fft_off,1),size(data_fft_off,1),true);

        % average across resampled trials
        data_base_boot(ii,:) = mean(data_fft_off(trial_set,:),1);
    end
    
    subplot(1,3,s),hold on
    fill([f; f(end:-1:1)],...
        [quantile(data_base_boot,.16) quantile(data_base_boot(:,end:-1:1),.84)],...
        [0 1 0]);%[.5 .5 .5])

    % mean baseline for fit
    data_base = mean(data_fft_off,1); % baseline
    [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width,fit_f2] = ...
        ecog_fitgamma(f,f_use4fit,data_base,data_base);
    plot(f,10.^(out_exp(2)-out_exp(1)*log10(f)),'k:','LineWidth',1)
    
    resamp_parms = NaN(length(stims_plot),6);
    for k = 1:length(stims_plot)
        % get stimulus data
        data_fit = data_fft(stims==stims_plot(k),:); % stimuli

        % fit stimulus data
        [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width,fit_f2] = ...
            ecog_fitgamma(f,f_use4fit,data_base,mean(data_fit,1));
        resamp_parms(k,1) = out_exp(1); % this is the baseline slope
        resamp_parms(k,2) = bb_amp;
        resamp_parms(k,3) = gamma_amp; % actual gamma amplitude is gamma_amp./gamma_width
        resamp_parms(k,4) = gamma_freq;
        resamp_parms(k,5) = gamma_width;
        resamp_parms(k,6) = out_exp(2); % this is the baseline intercept

        % plot fit to the mean
        plot(f,10.^fit_f2,':','Color',stims_color{k},'LineWidth',1)

        % bootstrap data and plot CI
        data_boot = zeros(nr_boots,size(data_fit,2));
        for ii = 1:nr_boots
            % from baseline data, random sample with replacement
            trial_set = randsample(size(data_fit,1),size(data_fit,1),true);

            % average across resampled trials
            data_boot(ii,:) = mean(data_fit(trial_set,:),1);

        end
        fill([f; f(end:-1:1)],...
            [quantile(data_boot,.84) quantile(data_boot(:,end:-1:1),.16)],...
            [.5 .5 .5],'EdgeColor',stims_color{k})
    end

    if isequal(subj,'19')
        xlim([30 200]),ylim([10.^-2 10.^1.6])
    elseif isequal(subj,'24')
        xlim([30 200]),ylim([10.^-1.5 10.^2.5])
    elseif isequal(subj,'1001')
        xlim([30 200]),ylim([10.^-1.2 10.^2.5])
    end

    set(gca,'XTick',[50 100 200],...
        'YTick',[10.^-1 10.^0 10.^1 10.^2])
    set(gca,'xscale','log','yscale','log','TickLength',[0.05 .1])
    xlabel('Frequency (Hz)'),ylabel('Power')
end
% 
% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',fullfile(dataDir,'derivatives','figures',...
%     ['Figure1B']))
% print('-depsc','-r300',fullfile(dataDir,'derivatives','figures',...
%     ['Figure1B']))
% 

