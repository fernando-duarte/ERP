% WRDS_INSTALL Install the WRDS API

% Add path of this folder
p = fileparts(mfilename('fullpath'));
addpath(p)

% Add path to ssh2
addpath(fullfile(p,'external','ssh2'))
addpath(fullfile(p,'external','passfield'))

% Save path
% savepath

% Make data path
if ~exist(fullfile(p,'data'),'dir')
    mkdir(fullfile(p,'data'))
end