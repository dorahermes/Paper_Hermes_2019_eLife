
% This is the code that was used to preprocess the ECoG data. This script
% is saved to show what we did before the ECoG data were deposited in OSF
% for the following publication:
%
%
% Hermes D, Petridou N, Kay K, Winawer J. 2019 An image-computable model
% for the stimulus selectivity of gamma oscillations. eLife 2019;8:e47035.
% DOI: https://doi.org/10.7554/eLife.47035
%
% dhermes 2018/2019 UMC Utrecht, Mayo Clinic

clear all

% set paths:
rootPath = gammaModelPath();
dataDir = fullfile(rootPath,'data');

sourceDir = '../../data/OVmodel/';

%% task-soc downsample and CAR

subjects = {'19','24','1001'};

% electrodes used in the paper 
electrodes = {[107 108 109 115 120 121],[45 46],[49 50 52 57 58 59 60]};

for s = 1:length(subjects)
    % subject name
    subj = subjects{s};
    dataName = dir(fullfile(sourceDir,['sub-' subj],...
        ['sub-' subj '_ses-01_task-soc_run-*_ieeg.eeg']));
    nr_runs = length(dataName);
    dataChannelsName = dir(fullfile(dataDir,['sub-' subj],['ses-01'],'ieeg',...
        ['sub-' subj '_ses-01_task-soc_run-*_channels.tsv']));
    
    %%%% THIS LOOP IS FOR PREPROCESSING DATA ACROSS RUNS
    for data_nr = 1:nr_runs 
        %%%% load raw data
        disp(['loading data set ' int2str(data_nr)])
        data = ft_read_data(fullfile(sourceDir,['sub-' subj],...
            dataName(data_nr).name));
        data_hdr = ft_read_header(fullfile(sourceDir,['sub-' subj],...
            dataName(data_nr).name));
        srate = round(data_hdr.Fs); % get sampling frequency

        %%%% get channel info 
        disp(['loading channels.tsv set ' int2str(data_nr)])
        channel_info = readtable(fullfile(dataDir,['sub-' subj],'ses-01','ieeg/',...
            dataChannelsName(data_nr).name),...
            'FileType','text','Delimiter','\t');
        % make a vector of channels to exclude from common average
        % referencing
        exclude_channels = zeros(length(channel_info.status),1);
        for k = 1:length(channel_info.status)
            if isequal(channel_info.status{k},'bad')  
                exclude_channels(k) = 1;
            end
            if ~isequal(channel_info.type{k},'ECOG') && ~isequal(channel_info.type{k},'SEEG')
                exclude_channels(k) = 1;
            end
        end
        include_channels = find(exclude_channels==0);
        exclude_channels = find(exclude_channels>0);
        
        % rereference
        disp(['CAR'])
        [data] = ecog_CarRegress(data,include_channels);
       
        % resample / cut awkward sampling rates down to 1000 Hz
        if ismember(subj,{'19','24'})
            sr_temp = srate*(1000/1031)*(1000/1480); 
            disp(['downsampling to ' int2str(sr_temp) ' Hz'])
            data = resample(data',1000,1031)';
            data = resample(data',1000,1480)';
            srate = floor(sr_temp); clear sr_temp
            disp('data downsampled')
        end
        
        % only save the included channels on V1/V2/V3 with pRFs within the
        % stimulus
        data = data(electrodes{s},:);
        electrodes_incl = electrodes{s};
        
        [~,a,~] = fileparts(dataName(data_nr).name);
        dataNamePrep = fullfile(dataDir,'derivatives','preprocessing',['sub-' subj],'ses-01','ieeg',...
            [a '_preproc.mat']);

        %%%% save data
        save(dataNamePrep,'-v7.3','data','srate','electrodes_incl')
        clear data
        
    end

end

