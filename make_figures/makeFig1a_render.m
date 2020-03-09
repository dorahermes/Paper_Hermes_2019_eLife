%
% This script will generate Figure 1A from: 
%
% Hermes D, Petridou N, Kay K, Winawer J. 2019 An image-computable model
% for the stimulus selectivity of gamma oscillations. eLife 2019;8:e47035.
% DOI: https://doi.org/10.7554/eLife.47035
% 
%
%
% The cortical surface is rendered with visual maps extracted by the
% updated Benson atlas. We add ECoG electrodes that were used for analyses
% (located on V1, V2, V3 with receptive fields within the presented images)
% 
% Citation for updated Benson maps:
% Benson and Winawer. eLife 2018;7:e40224. DOI: https://doi.org/10.7554/eLife.40224
% 
% Original citation for Benson maps:
% Noah C. Benson, Omar H. Butt, Ritobrato Datta, Petya D. Radoeva, David H. Brainard, Geoffrey K. Aguirre
% The Retinotopic Organization of Striate Cortex Is Well Predicted by Surface Topology
% Current Biology, Volume 22, Issue 23, 4 December 2012, Pages 2284
% 
% dhermes 2018 

% clear all
clearvars

% set paths:
rootPath = gammaModelPath();
dataDir = fullfile(rootPath,'data');

%% Render Benson Areas with electrodes

subjects = {'19','24','1001'};
hemi = {'L','R','R'};
hemi_s = {'l','r','r'};

% best viewing angle for left/right hemisphere where we can see electrodes
% for each subject
popout = true; % 'false' = turned off
if popout == true
    warning('off', 'backtrace') % rjc
    warning('electrodes pop out in the viewing direction, do not rotate the brain using the mouse, use the code instead or turn off popout')
end % if
v_dirs = [64 -23;270 0;-60 -15];

% electrodes used in the paper 
electrodes = {[107 108 109 115 120 121],[45 46],[49 50 52 57 58 59 60]};
electrodes_label = {["107" "108" "109" "115" "120" "121"],["45" "46"],...
    ["49" "50" "52" "57" "58" "59" "60"]};

Benson_Area_Names = {'V1','V2','V3','hV4','V01','V02','V3a','V3b','LO1','LO2','TO1','T02'};

for s = 1:length(subjects)
    % subject code
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataDir,'derivatives','surfaces',['sub-' subj],...
        ['sub-' subj '_ses-01_T1w_pial.' hemi{s} '.surf.gii']);
    % load gifti:
    g = gifti(dataGiiName);

    % electrode locations name:
    dataLocName = fullfile(dataDir,['sub-' subj],'ses-01','ieeg',...
        ['sub-' subj '_ses-01_electrodes.tsv']);
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];
    
    % surface labels 
    surface_labels_name = fullfile(dataDir,'derivatives','surfaces',['sub-' subj],...
        [hemi_s{s} 'h.benson14_varea.mgz']);
    surface_labels = MRIread(surface_labels_name);
    vert_label = surface_labels.vol(:);

    % load the color map
    cmap = make_BensonColormap;
    
    % get viewing angle for this subject
    v_d = v_dirs(s,:);

    % make sure electrodes pop out
    % don't rotate the brains if popout is true
    b_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
    if popout == true
        a_offset = b_offset;
        b_offset = 2*a_offset;
    else
        a_offset = 0;
    end % if
    els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1); % electrode position
    els_label = elecmatrix+repmat(b_offset,size(elecmatrix,1),1);

    figure
    ecog_RenderGiftiLabels(g,vert_label,cmap,Benson_Area_Names)

    el_add(els(electrodes{s},:),[.99 .99 .99],40) % add electrode positions
    el_add(els(electrodes{s},:),'k',30) % add electrode positions
    ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
    set(gcf,'PaperPositionMode','auto')
%     print('-dpng','-r300',fullfile(dataDir,'derivatives','render',['sub-' subj],...
%         ['sub-' subj '_BensonAreas__v' int2str(v_d(1)) '_' int2str(v_d(2))]))
end
   
% [EOF]