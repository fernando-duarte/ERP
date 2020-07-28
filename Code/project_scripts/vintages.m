%% Add Vintage


%%% this file adds the results to the ERP_vintages CSV if they are
%%% different than previous results
stateWarning = warning('off','all');

%% set directories

% main directories
excel_output = [root_dir filesep 'Output'];

%todays date
datestring = datestr(today(), 'mmddyyyy');

% Loading in New Data
filename = [excel_output filesep 'ERP_principal_Components_1yr.csv'];
erp = readtable(filename, 'DatetimeType', 'text');
erp = erp(:, 1:2);
erp.Properties.VariableNames = {'date', ['erp_' datestring]};
erp.date = datenum(erp.date);

%%% loading in vintage data
filename = [excel_output filesep 'ERP_vintages.csv'];
erp_vintages = readtable(filename, 'DatetimeType','text');
headers = erp_vintages.Properties.VariableNames;
% headers(1) = {'date'};
erp_vintages.Properties.VariableNames = headers;
erp_vintages.date = datenum(erp_vintages.date);

%%% compare columns of erp_vintages to erp_all
data_vintages = table2array(erp_vintages(:,2:end));
data_erp = erp{:,2};

issame = 0;
if size(data_erp,1) == size(data_vintages,1) %%% comparing number of observations
    for i = 1:size(data_vintages,2) %%% looping through each column in vintage data and seeing if identical
        ind = (data_erp == data_vintages(:,i));
        if sum(ind) == size(ind,1)
            issame = 1;
        end
    end
end

%%% if something has changed, then add new column otherwise keep the same
if issame == 0
    erp_vintages = outerjoin(erp, erp_vintages, 'Keys', 'date');
    erp_vintages.date_erp_vintages = []; %%% deleting added date column
    erp_vintages = erp_vintages(:, [1 3:end 2]); %%% moving new column to the end
    headers(1,end+1) = cellstr(['erp_' datestr(today(), 'MMDDYYYY')]); 
end

%%% replacing NaN with blank
erp_vintages = table2cell(erp_vintages);

for i = 2:size(erp_vintages,2)
    temp = erp_vintages(:,i);
    temp(cellfun(@isnan,temp)) = {[]};
    erp_vintages(:,i) = temp;
end
erp_vintages = cell2table(erp_vintages);
erp_vintages.Properties.VariableNames = headers;
erp_vintages.date = datestr(erp_vintages.date);

%%% exporting to excel
csv_name = [excel_output filesep 'ERP_vintages.csv'];
writetable(erp_vintages, csv_name);

clearvars -except root_dir
disp('Finished vintages.m')