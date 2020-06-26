%% test_predictability

%%% this file finds the out-of-sample (oos) R^2 statistic of a regression
%%% of realized returns on lagged ERP estimates to assess their 
%%% predictability

stateWarning = warning('off','all');

%% set directories

% main directories
input_dir = [root_dir filesep 'Temp' filesep 'Test for predictability' filesep 'Input'];
input_dir2 = [root_dir filesep 'Temp' filesep 'Create ERP measures' filesep 'Input'];
output_dir = [root_dir filesep 'Temp' filesep 'Test for predictability' filesep 'Output'];
excel_output = [root_dir filesep 'Output'];

% make directories if they don't exist
dir_out = {output_dir};
for i=1:length(dir_out)
    if ~exist(dir_out{i},'dir')
        mkdir(dir_out{i});
    end
end
clear i 

%% load data
load([input_dir filesep 'erp_principal_components.mat'])
load([input_dir filesep 'erp_measures.mat'])
load([input_dir2 filesep 'variables.mat'])

%% test predictability of erp models

for geometric_return = [0 1]
    
    % set up labels to name variables
    if geometric_return==1
        return_label = 'geometric';
    else
        return_label = 'arithmetic';
    end
    
    for return_horizon = [1/12 1/4 1/2 1 2 3 4 5]
        
        % set up labels to name variables
        switch return_horizon
            case 1/12
                horizon_label = 'one_month_ahead';
            case 1/4
                horizon_label = 'one_quarter_ahead';
            case 1/2
                horizon_label = 'six_months_ahead';
            case 1
                horizon_label = 'one_year_ahead';
            case 2
                horizon_label = 'two_years_ahead';
            case 3
                horizon_label = 'three_years_ahead';
            case 4
                horizon_label = 'four_years_ahead';
            case 5
                horizon_label = 'five_years_ahead';
        end
        
        % run predictive regressions
        realized_erp_ts = change_horizon(data_ts.realized_erp, data_ts.Time, 1/12,return_horizon,geometric_return, {'realized_erp'}); % change horizon by computing geometric or arithmetic multi-period returns
        predictors = fieldnames(erp_ts.(return_label).(horizon_label))';
        predictors = predictors(1:20);
        
        clear predicted_x dates r2 coefficients residual bint X
        for variable = predictors
            Xvar = timetable(erp_ts.(return_label).(horizon_label).Time, erp_ts.(return_label).(horizon_label).(char(variable)));
            Xvar.Properties.VariableNames = variable;
            XY_ts = synchronize(realized_erp_ts,Xvar);
            X.(variable{:}) = XY_ts.(char(variable)); % said variable is right hand side variable
            Y = XY_ts.realized_erp; % realized erp is left hand side variable
            [predicted_x.(variable{:}), dates.(variable{:}), r2.(variable{:}), coefficients, residual.(variable{:}), bint.(variable{:})] = predictive_regression(Y,X.(variable{:}),XY_ts.Time,12*return_horizon);
        end
        % find which was the best model in real time (using oos R^2)
        cellpredictedx = struct2cell(predicted_x); cellpredictedx_oos = cellfun(@(x) x.oos,cellpredictedx,'UniformOutput',0);
        celldates = struct2cell(dates);
        cellr2 = struct2cell(r2); cellr2_oos = cellfun(@(x) x.oos,cellr2,'UniformOutput',0);
        cell_ts_pred = cellfun(@(ts_date,ts_erp,ts_name) timetable(ts_date,ts_erp),celldates,cellpredictedx_oos,predictors','UniformOutput',0);
        cell_ts_r2 = cellfun(@(ts_date,ts_r2,ts_name) timetable(ts_date,ts_r2),celldates,cellr2_oos,predictors','UniformOutput',0);
        predictor_ts = synchronize(cell_ts_pred{:});
        r2_ts = synchronize(cell_ts_r2{:});
        r2_ts.Properties.VariableNames = predictors;
        
        [~,best_model_index]=max(table2array(r2_ts),[],2);best_model = predictors(best_model_index(end)); best_model_desc = model_description(strcmpi(predictors,best_model));
        [~,best_hist_mean_index]=max(table2array(r2_ts(:,1:2)),[],2);best_hist_mean = hist_mean(best_hist_mean_index(end)); best_hist_mean_desc = model_description(strcmpi(predictors,best_hist_mean));
        [~,best_ddm_index]=max(table2array(r2_ts(:,3:10)),[],2);best_ddm = ddm(best_ddm_index(end)); best_ddm_desc = model_description(strcmpi(predictors,best_ddm));
        [~,best_cross_section_index]=max(table2array(r2_ts(:,11:14)),[],2);best_cross_section = cross_section(best_cross_section_index(end)); best_cross_section_desc = model_description(strcmpi(predictors,best_cross_section));
        [~,best_time_series_index]=max(table2array(r2_ts(:,15:19)),[],2);best_time_series = time_series(best_time_series_index(end)); best_time_series_desc = model_description(strcmpi(predictors,best_time_series));
        [~,best_survey_index]=max(table2array(r2_ts(:,20)),[],2);best_survey = surveys(best_survey_index(end)); best_survey_desc = model_description(strcmpi(predictors,best_survey));
        
        best_r2_desc.(return_label).(horizon_label).all = best_model_desc;
        best_r2_desc.(return_label).(horizon_label).hist_mean = best_hist_mean_desc;
        best_r2_desc.(return_label).(horizon_label).ddm = best_ddm_desc;
        best_r2_desc.(return_label).(horizon_label).cross_section = best_cross_section_desc;
        best_r2_desc.(return_label).(horizon_label).time_series = best_time_series_desc;
        best_r2_desc.(return_label).(horizon_label).surveys = best_survey_desc;
        
        % record oos R^2 for best models in each category
        best_r2.(return_label).(horizon_label).all = r2_ts.(best_model{:})(end);
        best_r2.(return_label).(horizon_label).hist_mean = r2_ts.(best_hist_mean{:})(end);
        best_r2.(return_label).(horizon_label).ddm = r2_ts.(best_ddm{:})(end);
        best_r2.(return_label).(horizon_label).cross_section = r2_ts.(best_cross_section{:})(end);
        best_r2.(return_label).(horizon_label).time_series = r2_ts.(best_time_series{:})(end);
        best_r2.(return_label).(horizon_label).surveys = r2_ts.(best_survey{:})(end);
        
    end
end

%% test predictability of principal components

for geometric_return = [0 1]
    
    % set up labels to name variables
    if geometric_return==1
        return_label = 'geometric';
    else
        return_label = 'arithmetic';
    end
    
    for return_horizon = [1/12 1/4 1/2 1 2 3 4 5]
        
        % set up labels to name variables
        switch return_horizon
            case 1/12
                horizon_label = 'one_month_ahead';
            case 1/4
                horizon_label = 'one_quarter_ahead';
            case 1/2
                horizon_label = 'six_months_ahead';
            case 1
                horizon_label = 'one_year_ahead';
            case 2
                horizon_label = 'two_years_ahead';
            case 3
                horizon_label = 'three_years_ahead';
            case 4
                horizon_label = 'four_years_ahead';
            case 5
                horizon_label = 'five_years_ahead';
        end
        
        % run predictive regressions
        realized_erp_ts = change_horizon(data_ts.realized_erp,data_ts.Time,1/12,return_horizon,geometric_return, {'realized_erp'}); % change horizon by computing geometric or arithmetic multi-period returns
        predictors = fieldnames(erp_pc.(return_label).(horizon_label))';
                
        clear predicted_x dates r2 coefficients residual bint cell*
        for variable = predictors
            Xvar = timetable(erp_pc.(return_label).(horizon_label).(char(variable)).Time, erp_pc.(return_label).(horizon_label).(char(variable)).erp);
            Xvar.Properties.VariableNames = variable;
            XY_ts = synchronize(realized_erp_ts,Xvar);
            X.(variable{:}) = XY_ts.(char(variable)); % said variable is right hand side variable
            Y = XY_ts.realized_erp; % realized erp is left hand side variable
            [predicted_x.(variable{:}), dates.(variable{:}), r2.(variable{:}), coefficients, residual.(variable{:}), bint.(variable{:})] = predictive_regression(Y,X.(variable{:}),XY_ts.Time,12*return_horizon);
        end
               
        % find which was the best model in real time (using oos R^2)
        cellpredictedx = struct2cell(predicted_x); 
        cellpredictedx_oos = cellfun(@(x) x.oos,cellpredictedx,'UniformOutput',0);
        celldates = struct2cell(dates);
        cellr2 = struct2cell(r2); cellr2_oos = cellfun(@(x) x.oos,cellr2,'UniformOutput',0);
        cell_ts_pred = cellfun(@(ts_date,ts_erp,ts_name) timetable(ts_date, ts_erp), celldates,cellpredictedx_oos,predictors','UniformOutput',0);
        cell_ts_r2 = cellfun(@(ts_date,ts_r2,ts_name) timetable(ts_date,ts_r2),celldates,cellr2_oos,predictors','UniformOutput',0);
        predictor_ts = synchronize(cell_ts_pred{:});
        r2_ts = synchronize(cell_ts_r2{:});
        
        % record oos R^2 for best models in each category
        best_r2_pc.(return_label).(horizon_label).all = table2array(r2_ts(:,1));
        best_r2_pc.(return_label).(horizon_label).hist_mean = table2array(r2_ts(:,2));
        best_r2_pc.(return_label).(horizon_label).ddm = table2array(r2_ts(:,3));
        best_r2_pc.(return_label).(horizon_label).cross_section = table2array(r2_ts(:,4));
        best_r2_pc.(return_label).(horizon_label).time_series = table2array(r2_ts(:,5));
        best_r2_pc.(return_label).(horizon_label).surveys = table2array(r2_ts(:,6));
        
    end
end

%% save
save([output_dir filesep 'predictability.mat'],'best_r2','best_r2_pc')

%% write csv

for geometric_return = [0 1]
    
    % set up labels to name variables
    if geometric_return==1
        return_label = 'geometric';
    else
        return_label = 'arithmetic';
    end
    
    % name of the excel file
    filename = [excel_output filesep 'Predictability_Models_' return_label '.csv'];
    
    % set up headers
    columnname ={'One_month','One_quarter','Six_months', 'One_year','Two_years','Three_years','Four_years','Five_years'};
      
    rownames = {'historical mean';'name1';'ddm';'name2';'cross sectional regressions';'name3';'time series regressions';'name4';'surveys';'zeros'};
    
    % save description of best models
    model_type = {'hist_mean','ddm','cross_section','time_series','surveys'};
    horizon_type ={'one_month_ahead','one_quarter_ahead','six_months_ahead','one_year_ahead','two_years_ahead','three_years_ahead','four_years_ahead','five_years_ahead'};
    r2_for_excel = cell(2*length(model_type),length(horizon_type));
    for i = 1:2:size(r2_for_excel,1)
        for j = 1:size(r2_for_excel,2)
            r2_for_excel{i,j} = round(best_r2.(return_label).(horizon_type{j}).(model_type{(i-1)/2+1}),2);
            r2_for_excel{i+1,j} = best_r2_desc.(return_label).(horizon_type{j}).(model_type{(i-1)/2+1}){1};
        end
    end
    
    % changing how NaNs are displayed
    for i = 1:2:size(r2_for_excel,1)
        temp = r2_for_excel(i,:);
        temp(cellfun(@isnan,temp)) = {[]};
        r2_for_excel(i,:) = temp;
    end
    
    % writing csv
    export_to_excel = cell2table(r2_for_excel);
    export_to_excel.Properties.VariableNames = columnname;
    export_to_excel.Properties.RowNames = rownames;
    writetable(export_to_excel, filename, 'WriteRowNames', true);
    
end

%% export R^2 of principal components to excel

for geometric_return = [0 1]
    
    % set up labels to name variables
    if geometric_return==1
        return_label = 'geometric';
    else
        return_label = 'arithmetic';
    end
    
    % name of the excel file
    filename = [excel_output filesep 'Predictability_PC_' return_label '.csv'];

    % set up headers
    columnname ={'One_month','One_quarter','Six_months', 'One_year','Two_years','Three_years','Four_years','Five_years'};
    rowname = {'all';'historical mean';'ddm';'cross sectional regressions';'time series regressions';'surveys'};
    
    % save description of best models
    model_type = {'all','hist_mean','ddm','cross_section','time_series','surveys'};
    horizon_type ={'one_month_ahead','one_quarter_ahead','six_months_ahead','one_year_ahead','two_years_ahead','three_years_ahead','four_years_ahead','five_years_ahead'};
    r2_for_excel = cell(length(model_type),length(horizon_type));
    for i = 1:size(r2_for_excel,1)
        for j = 1:size(r2_for_excel,2)
            r2_for_excel{i,j} = round(best_r2_pc.(return_label).(horizon_type{j}).(model_type{i}),2);
        end
    end
    
    export_to_excel = cell2table(r2_for_excel);
    export_to_excel.Properties.VariableNames = columnname;
    export_to_excel.Properties.RowNames = rowname;
    writetable(export_to_excel, filename, 'WriteRowNames', true);
    
end
warning(stateWarning); % restore warning state

clearvars -except root_dir
disp('Finished test_predictability.m')