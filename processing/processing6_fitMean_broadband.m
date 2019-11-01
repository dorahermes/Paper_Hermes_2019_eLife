%
% Mean model as control
%
% Dora Hermes, 2019

clear all

% set paths:
rootPath = gammaModelPath();
dataDir = fullfile(rootPath,'data');

%% Load ECoG data and fit

%%%%% Pick a subject:
subjects = {'19','24','1001'};

for s = 1:3
    subj = subjects{s};

    if isequal(subj,'19') % S1
        im_deg = rad2deg(atan(17.9./50));
        electrodes = [107 108 109 115 120 121]; 
    elseif isequal(subj,'24') % S2
        im_deg = rad2deg(atan(17.9./45));
        electrodes = [45 46]; 
    elseif isequal(subj,'1001') % S3
        im_deg = rad2deg(atan(17.9./50));
        electrodes = [49 50 52 57 58 59 60]; 
    end
    
    for el = 1:length(electrodes)
        elec = electrodes(el);

        % Choose an analysis type:
        analysisType = 'spectra200';

        % load ecog stimulus data (1000 bootstraps):
        dataFitName = fullfile(dataDir,'derivatives','preprocessing',...
            ['sub-' subj '_ses-01_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
        load(dataFitName)

        % Broadband power percent signal change per stimulus:
        bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
        ecog_bb = 100*(10.^(resamp_parms(:,:,2)-bb_base)-1);

        %% Fit mean model

        % Initalize outputs:
        cross_MeanPred = zeros(size(ecog_bb,1),1);

        % Estimate gain, leave one out every time for cross validation
        for kk = 1:size(ecog_bb,1) % number of stimuli, leave out kk   
            % training stimuli (kk is left out)
            trainSet = setdiff([1:size(ecog_bb,1)],kk);

            % get the mean of training stimuli
            kkEstimate = mean(mean(ecog_bb(trainSet,:),2),1);

            % this is the estimate for the leftout stimulus
            cross_MeanPred(kk,1) = kkEstimate;
        end

        save(fullfile(dataDir,'derivatives','gaborFilt','mean_broadband',...
            ['sub' subj '_el' int2str(elec) '_' analysisType '_MeanBB']),...
            'cross_MeanPred')    
    end
end

