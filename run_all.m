%% run_all

%%% This file runs the Equity Risk Premium (ERP) project. It automatically 
%%% downloads and updates the data, computes specific ERP measures as 
%%% presented in the literature uses those to compute the ERP at 
%%% different horizons. 

clear 
clc

%% set up directories
%main directory
root_dir = [filesep 'home' filesep 'rcecrs02' filesep 'ERPv7']; % pick the folder where this file is located
cd(root_dir) %%% do not change directory after running this line or code will fail when running wrds downloads in dwnload_data

%%% adding paths
addpath([root_dir filesep 'Code'])
addpath([root_dir filesep 'Input'])
addpath([root_dir filesep 'Temp'])
addpath([root_dir filesep 'Output'])
addpath(genpath([root_dir filesep 'Code' filesep 'functions']))
addpath(genpath([root_dir filesep 'Code' filesep 'Matlab subfunctions']))
addpath(genpath([root_dir filesep 'Code' filesep 'Matlab subfunctions' filesep 'Factor models']))
addpath(genpath([root_dir filesep 'Code' filesep 'Matlab subfunctions' filesep 'FredFetch-master']))
addpath(genpath([root_dir filesep 'Code' filesep 'Matlab subfunctions' filesep 'lightspeed']))
addpath(genpath([root_dir filesep 'Code' filesep 'Matlab subfunctions' filesep 'WRDS API']))
addpath(genpath([root_dir filesep 'Code' filesep 'NYFed_functions']))
addpath(genpath([root_dir filesep 'Code' filesep 'project_scripts']))
addpath(genpath([root_dir filesep 'Code' filesep 'stata_scripts']))
addpath(genpath([root_dir filesep 'Code' filesep 'wrds_scripts']))

%% running the project code
run('download_data.m') %ran
run('NY_Fed_CM_DAPM.m') %ran 
run('Duke_CFO.m') %ran
run('cay.m') %ran 
run('collect_data.m') %ran
run('create_variables.m') %ran
run('create_ERP_measures.m') %ran
run('erp_principal_components.m') %ran
run('create_term_structure.m') %ran
run('test_predictability.m') %ran
run('vintages.m') %ran