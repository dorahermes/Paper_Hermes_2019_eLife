%
% This script will generate the time-frequency plots for example electrodes
% in Figure 1c from: 
%
% Hermes D, Petridou N, Kay K, Winawer J. 2019 An image-computable model
% for the stimulus selectivity of gamma oscillations. eLife 2019;8:e47035.
% DOI: https://doi.org/10.7554/eLife.47035
%
% dhermes 2019 

clear all

% set paths:
rootPath = gammaModelPath();
dataDir = fullfile(rootPath,'data');

%% load data and epoch

subjects = {'19','24','1001'};
electrodes = {109,45,50}; % example electrodes per subject

% Conditions to plot: grating, noise and blank (blank to get colorbar)
conds_plot = {[45],[83],87};

% Make a nice colormap
cm1 = [repmat([0 0 0],100,1)];
cm1(1:40,1) = [0.7]';
cm1(1:40,2) = [0.7:-0.6/39:0.1]';
cm1(1:40,3) = [0.7:-0.7/39:0]';
cm1(40:100,1) = [0.7:(1-0.7)/60:1]';
cm1(40:100,2) = [0.1:.9/60:1]';
cm2 = [repmat([0 0 0],100,1)];
cm2(1:30,3) = [0.7]';
cm2(1:30,1) = [0.7:-0.7/29:0]';
cm2(1:30,2) = [0.7:-0.7/29:0]';
cm2(30:100,3) = [0.7:(1-0.7)/70:1]';
cm2(30:100,2) = [0:1/70:1]';
cm = [cm2(end:-1:1,:); cm1];

for s = 1:length(subjects)

    % subject name
    subj = subjects{s};
    dataName = dir(fullfile(dataDir,'derivatives','preprocessing',['sub-' subj],'ses-01','ieeg',...
        ['sub-' subj '_ses-01_task-soc_run-*_ieeg_preproc.mat']));
    eventsName = dir(fullfile(dataDir,['sub-' subj],'ses-01','ieeg',...
        ['sub-' subj '_ses-01_task-soc_run-*_events.tsv']));
    nr_runs = length(dataName);

    data_epoch_all = [];
    stim_all = [];
    
    % loop runs and add all to one matrix
    for data_nr = 1:nr_runs
        %%%% load preprocessed data
        load(fullfile(dataName(data_nr).folder,dataName(data_nr).name));

        %%%% load stimulus information
        stim = readtable(fullfile(eventsName(data_nr).folder,eventsName(data_nr).name),...
            'FileType','text','Delimiter','\t');

        %%%% notch filter data at 60, 120 and 180 Hz
        if ismember(subj,{'19','24'})
            data = ecog_notch(data',srate,60)';
        else
            disp('no notch filter, data are pretty noise-free')
        end

        %%%% make epochs
        onset_trial = round(stim.onset*srate);%from seconds to samples
        epoch_l = 1.5; % epoch length: -0.5:1 sec
        pre_stim = .5;

        data_epoch = zeros(size(data,1),length(onset_trial),epoch_l*srate);
        for elec = 1:size(data,1)
            for l = 1:length(onset_trial)
                data_epoch(elec,l,:) = ...
                    data(elec,onset_trial(l)-pre_stim*srate+1:onset_trial(l)+(epoch_l-pre_stim)*srate);
            end
        end
        % define t - time vector for each epoch
        t = [1:epoch_l*srate]/srate - pre_stim;  

        data_epoch_all = cat(2,data_epoch_all,data_epoch);
        stim_all = cat(1,stim_all,stim.trial_type);

    end
    data_epoch = data_epoch_all;
    clear data_epoch_all data stim

    %%%%% Plot ersp for an electrode for a grating and noise pattern
    elec = electrodes{s}; 
    elec_nr = find(electrodes_incl == elec);
    
    movingwin = [.200 .05];
    params.pad = -1;
    params.tapers = [3 5];
    params.fpass = [0 200];
    params.Fs = srate;
    params.trialave = 0;

    data_temp = squeeze(data_epoch(elec_nr,:,:));
    data_baseline = squeeze(data_epoch(elec_nr,stim_all==87,:));

    % calculate baseline
    params.trialave = 1;
    [S1b,t_tf_b,f_b] = mtspecgramc(data_baseline(:,t>.25 & t<.5)',movingwin,params);
    S1b = mean(S1b,1);

    params.trialave = 0;
    figure('Color',[1 1 1],'Position',[0 0 260 170])
    %%%%% all responses:
    for kk = 1:length(conds_plot)
        data2use = squeeze(data_temp(ismember(stim_all,conds_plot{kk}),:))';

        % calculate spectgram
        [S1,t_tf,f] = mtspecgramc(data2use,movingwin,params);
        t_tf = t_tf+t(1);
        % normalize wrt baseline
        S1 = nanmean(S1,3);
        S1 = S1./repmat(S1b,[size(S1,1),1]);

        subplot(1,length(conds_plot),kk)
        imagesc(t_tf,f,log10(S1)',[-1.5 1.5])
        axis xy
        colormap(cm)
        set(gca,'XTick',[0 .5])
        xlim([-.1 .9])
        title(['stim ' int2str(conds_plot{kk})])
    end

    colorbar

    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r300',fullfile(dataDir,'derivatives','figures',...
        ['Figure1C_sub-' subj '_el' int2str(elec)]))
    print('-depsc','-r300',fullfile(dataDir,'derivatives','figures',...
        ['Figure1C_sub-' subj '_el' int2str(elec)]))
end


%% 
%% phase locked power plot
%%

disp('Making Figure Supplement 1')

% https://groups.google.com/forum/#!topic/analyzingneuraltimeseriesdata/yqrV_rPu_OA
% From Mike X Cohen
% "Separating the phase-locked from the non-phase-locked signal is another
% area where the literature hasn't settled on a consistent solution.
% Nonetheless, the approach I used in the book is the most common one. The
% total signal is simply the time-frequency decomposition of each trial and
% then the average over trials. To compute the non-phase-locked signal, you
% first compute the ERP and then subtract the ERP from each trial, and then
% perform time-frequency decomposition of the residual (the assumption is
% that the ERP is the phase-locked signal, and thus by subtracting the
% phase-locked from the total, you are left with the non-phase-locked
% part). After computing the total and the non-phase-locked time-frequency
% power, subtract the non-phase-locked from the total. This gives you the
% phase-locked power. You might think that you should just take the
% time-frequency decomposition of the ERP, but this generally produces
% results that are very difficult to interpret."


%% load data and epoch

subjects = {'19','24','1001'};
electrodes = {109,45,50}; % example electrodes per subject

% Conditions to plot: grating, noise and blank (blank to get colorbar)
conds_plot = {[45],[83],87};

% Make a nice colormap
cm1 = [repmat([0 0 0],100,1)];
cm1(1:40,1) = [0.7]';
cm1(1:40,2) = [0.7:-0.6/39:0.1]';
cm1(1:40,3) = [0.7:-0.7/39:0]';
cm1(40:100,1) = [0.7:(1-0.7)/60:1]';
cm1(40:100,2) = [0.1:.9/60:1]';
cm2 = [repmat([0 0 0],100,1)];
cm2(1:30,3) = [0.7]';
cm2(1:30,1) = [0.7:-0.7/29:0]';
cm2(1:30,2) = [0.7:-0.7/29:0]';
cm2(30:100,3) = [0.7:(1-0.7)/70:1]';
cm2(30:100,2) = [0:1/70:1]';
cm = [cm2(end:-1:1,:); cm1];

for s = 1:length(subjects)

    % subject name
    subj = subjects{s};
    dataName = dir(fullfile(dataDir,'derivatives','preprocessing',['sub-' subj],'ses-01','ieeg',...
        ['sub-' subj '_ses-01_task-soc_run-*_ieeg_preproc.mat']));
    eventsName = dir(fullfile(dataDir,['sub-' subj],'ses-01','ieeg',...
        ['sub-' subj '_ses-01_task-soc_run-*_events.tsv']));
    nr_runs = length(dataName);

    data_epoch_all = [];
    stim_all = [];
    
    % loop runs and add all to one matrix
    for data_nr = 1:nr_runs
        %%%% load preprocessed data
        load(fullfile(dataName(data_nr).folder,dataName(data_nr).name));

        %%%% load stimulus information
        stim = readtable(fullfile(eventsName(data_nr).folder,eventsName(data_nr).name),...
            'FileType','text','Delimiter','\t');

        %%%% notch filter data at 60, 120 and 180 Hz
        if ismember(subj,{'19','24'})
            data = ecog_notch(data',srate,60)';
        else
            disp('no notch filter, data are pretty noise-free')
        end

        %%%% make epochs
        onset_trial = round(stim.onset*srate);%from seconds to samples
        epoch_l = 1.5; % epoch length: -0.5:1 sec
        pre_stim = .5;

        data_epoch = zeros(size(data,1),length(onset_trial),epoch_l*srate);
        for elec = 1:size(data,1)
            for l = 1:length(onset_trial)
                data_epoch(elec,l,:) = ...
                    data(elec,onset_trial(l)-pre_stim*srate+1:onset_trial(l)+(epoch_l-pre_stim)*srate);
            end
        end
                
        % define t - time vector for each epoch
        t = [1:epoch_l*srate]/srate - pre_stim;  

        data_epoch_all = cat(2,data_epoch_all,data_epoch);
        stim_all = cat(1,stim_all,stim.trial_type);

    end
    data_epoch = data_epoch_all;
    clear data_epoch_all data stim

    %%%%% Plot ersp for an electrode for a grating and noise pattern
    elec = electrodes{s}; 
    elec_nr = find(electrodes_incl == elec);

    % get data for this channel
    data_epoch_elec = data_epoch(elec_nr,:,:);
    
    % regress ERP out (includes baseline subtraction)
    t_base = find(t>-.3 & t<0);
    data_epoch_elecregerp = ecog_regresserp(data_epoch_elec,t_base,stim_all);
    
    % squeeze the data to be trials x time
    data_epoch_elec = squeeze(data_epoch_elec);
    data_epoch_elecregerp = squeeze(data_epoch_elecregerp);
    
    % specify multitaper parameters
    movingwin = [.200 .05];
    params.pad = -1;
    params.tapers = [3 5];
    params.fpass = [0 200];
    params.Fs = srate;
    params.trialave = 0;

    % calculate baseline spectrum for plotting
    data_baseline = data_epoch_elec(stim_all==87,:);
    data_baseline_regerp = data_epoch_elecregerp(stim_all==87,:);  
    params.trialave = 1;
    [S1b,~,~] = mtspecgramc(data_baseline(:,t>.25 & t<.5)',movingwin,params);
    [S1b_regerp,t_tf_b,f_b] = mtspecgramc(data_baseline_regerp(:,t>.25 & t<.5)',movingwin,params);


    params.trialave = 0;
    figure('Color',[1 1 1],'Position',[0 0 260 510])
    %%%%% all responses:
    for kk = 1:length(conds_plot)
        
        % 1) spectrum total
        data2use = data_epoch_elec(ismember(stim_all,conds_plot{kk}),:)';
        % calculate total spectgram
        [S1total,t_tf,f] = mtspecgramc(data2use,movingwin,params);
        t_tf = t_tf+t(1);
        S1total = nanmean(S1total,3);
        
        % 2) spectrum after regressing ERP
        data2use_regerp = data_epoch_elecregerp(ismember(stim_all,conds_plot{kk}),:)';       
        % calculate spectgram minus phaselocked
        [S1regerp,t_tf,f] = mtspecgramc(data2use_regerp,movingwin,params);
        t_tf = t_tf+t(1);
        S1regerp = nanmean(S1regerp,3);
        
        % normalize wrt baseline
        S1total = S1total./repmat(S1b,[size(S1total,1),1]);
        S1regerp = S1regerp./repmat(S1b_regerp,[size(S1regerp,1),1]);
                
        % plot total power
        subplot(3,length(conds_plot),kk)
        imagesc(t_tf,f,log10(S1total'),[-1.5 1.5])
        axis xy
        colormap(cm)
        set(gca,'XTick',[-5 5])
        xlim([-.1 .9])
        title(['total s' int2str(conds_plot{kk})])
        
        % plot non-phaselocked power
        subplot(3,length(conds_plot),3+kk)
        imagesc(t_tf,f,log10(S1regerp'),[-1.5 1.5])
        axis xy
        colormap(cm)
        set(gca,'XTick',[-5 5])
        xlim([-.1 .9])
        title(['nonphl s' int2str(conds_plot{kk})])

        % plot phase locked power
        % subtract non-phaselocked from total power
        S1phase = log10(S1total)-log10(S1regerp);
        
        subplot(3,length(conds_plot),6+kk)
        imagesc(t_tf,f,S1phase',[-1.5 1.5])
        axis xy
        colormap(cm)
        set(gca,'XTick',[0 .5])
        xlim([-.1 .9])
        title(['phl s' int2str(conds_plot{kk})])
    end

    colorbar

    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r300',fullfile(dataDir,'derivatives','figures',...
        ['Figure1_Supp_sub-' subj '_el' int2str(elec) '_phaselocked']))
    print('-depsc','-r300',fullfile(dataDir,'derivatives','figures',...
        ['Figure1_Supp_sub-' subj '_el' int2str(elec) '_phaselocked']))
end
