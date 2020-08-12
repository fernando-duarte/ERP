%% run_all

%%% This file runs the Equity Risk Premium (ERP) project. It automatically 
%%% downloads and updates the data, computes specific ERP measures as 
%%% presented in the literature uses those to compute the ERP at 
%%% different horizons. 

clear; clc;

%% set up directories 

% main directory
% pick the folder where this file is located
root_dir = [filesep 'home' filesep 'rcerxr21' filesep 'ERP']; 
cd(root_dir) % do not change directory after running this line or 
             % code will fail when running wrds downloads in download_data

% adding paths
addpath([root_dir filesep 'Code'])
addpath([root_dir filesep 'Input'])
addpath([root_dir filesep 'Temp'])
addpath([root_dir filesep 'Output'])
addpath(genpath([root_dir filesep 'Code' filesep 'functions']))
addpath(genpath([root_dir filesep 'Code' filesep 'Matlab subfunctions']))
addpath(genpath([root_dir filesep 'Code' filesep 'Matlab subfunctions' ...
    filesep 'Factor models']))
addpath(genpath([root_dir filesep 'Code' filesep 'Matlab subfunctions' ...
    filesep 'FredFetch-master']))
addpath(genpath([root_dir filesep 'Code' filesep 'Matlab subfunctions' ...
    filesep 'lightspeed']))
addpath(genpath([root_dir filesep 'Code' filesep 'Matlab subfunctions' ...
    filesep 'WRDS API']))
addpath(genpath([root_dir filesep 'Code' filesep 'NYFed_functions']))
addpath(genpath([root_dir filesep 'Code' filesep 'project_scripts']))
addpath(genpath([root_dir filesep 'Code' filesep 'stata_scripts']))
addpath(genpath([root_dir filesep 'Code' filesep 'wrds_scripts']))


%% creating initialization variables 

% FRED api key 
api = '35ee7b54bc4bf63e6e61903112994800';

% writes to the api.txt file with the provided API key
apiTXT = fopen('Input/api.txt', 'wt');
fprintf(apiTXT,  '%s', api); fclose(apiTXT);

% Downloadable url's and folder specifications

% pathing in Sas Studio WRDS
wrds_code_dir = '/home/frb-ny/fedraj/WRDS/Code';
wrds_out_dir = '/home/frb-ny/fedraj/WRDS/Output';

% Connect to the WRDS server
wrds_username = 'fedraj';
wrds_password = 'Feddata2020$%';

% The Bureau of Economic Analysis website 
beaURL = 'https://apps.bea.gov/histdata/Releases/GDP_and_PI/2020/Q1/Third_June-26-2020/Section2all_xls.xlsx';

% exporting variable declarations
save Temp/Init.mat wrds_code_dir wrds_out_dir wrds_username ...
    wrds_password beaURL

%% running the project code (RAN)

run('download_data.m') 
run('NY_Fed_CM_DAPM.m')  
run('Duke_CFO.m') 
run('cay.m') 
run('collect_data.m')
run('create_variables.m')
run('create_ERP_measures.m')
run('erp_principal_components.m')
run('create_term_structure.m')
run('test_predictability.m')
run('vintages.m')