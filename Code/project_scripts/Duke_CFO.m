%% Duke_CFO

%%% This file computes the ERP using data from Dke's Survey of CFOs.
%%% Specifically it uses their prediction of S&P 500 returns and prevaling
%%% interest rates. 

%% Set Directory 
input_dir = [root_dir filesep 'Input'];

%% Loading in data

%%% Hand entered duke cfo survey data
file_in = [input_dir filesep 'cfo_data.csv'];
cfo_data = readtable(file_in);

%%% deleting data with nans
ind = isnan(table2array(cfo_data(:,4)));
cfo_data = cfo_data(ind == 0, :);

date_cfo = datenum(cfo_data.SurveyDate_1, 'YYYYQQ');

returns1yr = table2array(cfo_data(:,4));
returns10yr = table2array(cfo_data(:,6));

%%% Bond yields from FRED

fred = fred.latest({'WGS1YR','WGS10YR'});

date_fred = fred.date;
bondyield1yr = fred.value(:,1);
bondyield10yr = fred.value(:,2);

% data = [bondyield1yr, bondyield10yr];

%%% creating quarterly values and matching dates

header = {'WGS1YR', 'WGS10YR'};
fred_ts = timetable(datetime(date_fred, 'ConvertFrom', 'datenum'), ...
    bondyield1yr, bondyield10yr, 'VariableNames', header);
fred_quarterly = retime(fred_ts,'quarterly', 'previous');

% data = [returns1yr returns10yr]; 
header = {'ExpRet1YR', 'ExpRet10YR'};
cfo_ts = timetable(datetime(date_cfo, 'ConvertFrom', 'datenum'), ...
    returns1yr, returns10yr, 'VariableNames', header);
cfo_quarterly = retime(cfo_ts,'quarterly', 'nearest');

%%% creating data to export 

data_out = synchronize(cfo_quarterly, fred_quarterly);
data_out = data_out(~isnan(data_out.ExpRet1YR), :);

%%% exporting data

file_out = [input_dir filesep 'cfo_data_final.mat'];
save(file_out, 'data_out')

clearvars -except root_dir

disp('Finished Duke_CFO.m')
