function rootPath = gammaModelPath()
% Set directory path for Gamma Model project
%
%     gammaModelPath
%
% Set up the path to the functions called by Gamma Model functions.
% Put data from OSF in rootPath/data/

rootPath = fileparts(which(mfilename));
fprintf('gamma model code directory: %s\n',rootPath)

% Adds the root directory of the Gamma Model tree to the user's path
addpath(rootPath);

% Generates a list of the directories below the Gamma Model tree.
addpath(fullfile(rootPath, 'compute'))
addpath(fullfile(rootPath, 'data_info'))
addpath(fullfile(rootPath, 'render'))

return;

