function gammaModelDownloadData()
% Download and unzip the data from the Open Science Framework project page
% associated with this paper:
%
%   Hermes D, Petridou N, Kay K, Winawer J. 2019 An image-computable model
%   for the stimulus selectivity of gamma oscillations. eLife 2019;8:e47035.
%   DOI: https://doi.org/10.7554/eLife.47035
% 
% Alternatively, the data can be downloaded manually from this site:
% https://osf.io/eqjxb/
%
% The code downloads a single zip file (2.03GB), places it in the root
% directory of the project, and unzips it into the folder named 'data'


url = 'https://osf.io/q4yad/download';

codeDir = gammaModelPath;
dataDir = fullfile(codeDir,'data');

if ~exist(dataDir,'dir')
    disp(['making data directory: ' dataDir])
    mkdir(dataDir);
end

pth = fullfile(dataDir, 'Hermes2019eLifeData.zip');

disp('start downloading data from OSF')
fname = websave(pth, url);

upzip(pth)

end