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

for s = 1:length(subjects) 
    subj = subjects{s};

    if isequal(subj,'19') % S1
        im_deg = rad2deg(atan(17.9./50));
        electrodes = [107 108 109 115 120 121]; 
    elseif isequal(subj,'24') % S2
        im_deg = rad2deg(atan(17.9./45));
        electrodes = [45 46];
    elseif isequal(subj,'1001') % S3
        im_deg = rad2deg(atan(20./60));
        electrodes = [49 50 52 57 58 59 60]; 
    end

    for el = 1:length(electrodes)
        elec = electrodes(el);

        % Choose an analysis type:
        analysisType = 'spectra200';

        % Load ecog stimulus data:
        dataFitName = fullfile(dataDir,'derivatives','preprocessing',...
            ['sub-' subj '_ses-01_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
        load(dataFitName)

        % Gamma power percent signal change:
        ecog_g = 100*(mean(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1,2));
    
        %% Fit mean model

        % Initalize outputs:
        cross_MeanPred = zeros(size(ecog_g,1),1);
        
            % Estimate gain, leave one out every time for cross validation
            for kk = 1:size(ecog_g,1) % number of stimuli, leave out kk   
                % training stimuli (kk is left out)
                trainSet = setdiff([1:size(ecog_g,1)],kk);

                % get the mean of training stimuli
                kkEstimate = mean(mean(ecog_g(trainSet,:),2),1);

                % this is the estimate for the leftout stimulus
                cross_MeanPred(kk,1) = kkEstimate;
            end

        save(fullfile(dataDir,'derivatives','gaborFilt','mean_gamma',...
            ['sub' subj '_el' int2str(elec) '_' analysisType '_MeanGamma']),...
            'cross_MeanPred')    
        
        disp(['Done fitting mean model for subj ' subj ' el ' int2str(elec)])
    end
end
