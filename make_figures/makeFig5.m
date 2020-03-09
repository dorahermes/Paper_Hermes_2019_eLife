%
% This script make figure 5 from:
%
% Hermes D, Petridou N, Kay K, Winawer J. 2019 An image-computable model
% for the stimulus selectivity of gamma oscillations. eLife 2019;8:e47035.
% DOI: https://doi.org/10.7554/eLife.47035
%
% 
%%
%% Visualize energy per orientation
%
% clear all
clearvars

% set paths:
rootPath = gammaModelPath();
dataDir = fullfile(rootPath,'data');

% load images
load(fullfile(dataDir,'stimuli','task-soc_stimuli.mat'),'stimuli')
load(fullfile(dataDir,'derivatives','gaborFilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

stimulus = sqrt(blob(stimulus.^2,2,2)); % quadrature pairs
% sum across orientation.  after this step, stimulus is images x positions.
imEnergyMean = blob(stimulus,2,8);

numOrientations = 8;

res = sqrt(size(stimulus,2)/numOrientations);

figure('Position',[0 0 250 500]),
im_nrs = [50 49 54 10];

for ii = 1:length(im_nrs)
    % Select an image number:
    inNr = im_nrs(ii);
    thisIm = reshape(stimulus(inNr,:),numOrientations,res,res); 

    % Imagine this pRF
    pp = [res/2 res/2 res/5];
    [~,xx,yy] = makegaussian2d(res,2,2,2,2);
    gaufun1 = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),...
        pp(3),xx,yy,0,0)/(2*pi*pp(3)^2)); % Gaussian or pRF
    imEnergyPrf = imEnergyMean(inNr,:)'.*gaufun1(pp);

    % Make figure for energy for every orientation
    OrientationEnergy = zeros(1,numOrientations);
    for kk = 1:numOrientations
        thisImPrf = thisIm(kk,:,:);
        thisImPrf = thisImPrf(:).*gaufun1(pp);
        % Get energy for this orientation
        OrientationEnergy(kk) = sum(thisImPrf);
    end
    
    orient_order = [3 4 5 6 7 8 1 2];

    % Plot original image and energy across orientations
    subplot(length(im_nrs),2,ii*2-1)
    mx = max(abs(stimuli(:)));
    imagesc(stimuli(:,:,inNr))
    axis image tight off;
    caxis([0 mx]);
    colormap(gray);

    subplot(length(im_nrs),2,ii*2)
    bar(OrientationEnergy(orient_order),'w')
    xlim([0 9]),ylim([0 .2])
    box off
    ylabel('energy')
    title(['im ' int2str(inNr) ' sd=' num2str(var(OrientationEnergy).^.5,3)])
    set(gca,'FontName','Arial','FontSize',8)

end
 
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',fullfile(dataDir,'derivatives','figures',...
    'Fig5'))
print('-depsc','-r300',fullfile(dataDir,'derivatives','figures',...
    'Fig5'))

% [EOF]