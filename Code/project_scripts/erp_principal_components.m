%% erp_principal_components

%%% This file constructs measures of the ERP by taking the principal
%%% components of the ERP estimates of different models.

stateWarning = warning('off','all');

%% set up parameters
start_date = datenum(1960,1,1); % pick initial date to use
end_date = datenum(datetime('today')); % pick final date to use

%% set directories

% main directories
pc_dir = [root_dir filesep 'Temp' filesep 'Principal components'];
input_dir = [root_dir filesep 'Temp' filesep 'Principal components' ...
    filesep 'Input'];
input_dir2 = [root_dir filesep 'Temp' filesep 'Create ERP measures' ...
    filesep 'Input'];
output_dir = [root_dir filesep 'Temp' filesep 'Principal components' ...
    filesep 'Output'];
next_dir = [root_dir filesep 'Temp' filesep 'Create term structure of ERP' ...
    filesep 'Input'];
next_dir2 = [root_dir filesep 'Temp' filesep 'Test for predictability' ...
    filesep 'Input'];
excel_output = [root_dir filesep 'Output'];

% make directories if they don't exist
dir_out = {output_dir,next_dir,next_dir2};

for i=1:length(dir_out)
    if ~exist(dir_out{i},'dir')
        mkdir(dir_out{i});
    end
end

clear i 
addpath(pc_dir)

%% load data
load([input_dir filesep 'erp_measures.mat'])
load([input_dir2 filesep 'variables.mat'])

%% find principal components

for return_horizon = [1/12, 1/4, 1/2, 1, 2, 3, 4, 5]
    
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
    
    % use geometric or arithmetic returns
    for geometric_return = [0 1]
    
        % set up labels to name variables
        if geometric_return==1
            return_label = 'geometric';
        else
            return_label = 'arithmetic';
        end
        
        current_spec = erp_ts.(return_label).(horizon_label); % select current specification
        num_model = sum(~isnan(current_spec{:,:}),2);
        keep1 = num_model > 6;
        keep2 = datenum(current_spec.Time)>=start_date & datenum(current_spec.Time)<=end_date; % select dates
        keep = logical(keep1 .* keep2);
     
        % find principal component of all models
        num_of_PC = 1;
        [pc eigenvectors share_var xs_mean xs_std coeff]  = PCA_unbalanced(table2array(current_spec(keep,:)),num_of_PC,1e-4,1e5,1,1,1);
        erp_pc.(return_label).(horizon_label).all = timetable(current_spec.Time(keep), ...
            pc(:,1)+nanmean(xs_mean), xs_mean, xs_std, pc(:,1), 'VariableNames', {'erp' 'xs_mean' 'xs_std' 'pc1'});
        PCA_info.(return_label).(horizon_label).all.coeff = coeff;
        PCA_info.(return_label).(horizon_label).all.share_var = share_var;
        
        models = table2array(current_spec(keep,:));
        model_loadings_on_PC = nan(size(models,2),num_of_PC);
        for q=1:size(models,2)
            model_loadings_on_PC(q,:) = regress(models(:,q),pc);
        end
        PCA_info.(return_label).(horizon_label).all.model_loadings_on_PC = model_loadings_on_PC;
        
        % find principal component of historical mean
        hist_mean_ts = [];
        for i = 1:length(hist_mean)
            hist_mean_ts = [hist_mean_ts current_spec.(hist_mean{i})];
        end
        [pc eigenvectors share_var xs_mean xs_std coeff]  = PCA_unbalanced(hist_mean_ts(keep,:),num_of_PC,1e-4,1e5,1,1,1);
        erp_pc.(return_label).(horizon_label).hist_mean = timetable(current_spec.Time(keep), ...
            pc(:,1)+nanmean(xs_mean), xs_mean, xs_std, pc(:,1), 'VariableNames', {'erp' 'xs_mean' 'xs_std' 'pc1'});
        PCA_info.(return_label).(horizon_label).hist_mean.coeff = coeff;
        PCA_info.(return_label).(horizon_label).hist_mean.share_var = share_var;
        
        % find principal component of ddm
        ddm_ts = [];
        for i = 1:length(ddm)
            ddm_ts = [ddm_ts current_spec.(ddm{i})];
        end
        [pc eigenvectors share_var xs_mean xs_std coeff]  = PCA_unbalanced(ddm_ts(keep,:),num_of_PC,1e-4,1e5,1,1,1);
        erp_pc.(return_label).(horizon_label).ddm = timetable(current_spec.Time(keep), ...
            pc(:,1)+nanmean(xs_mean), xs_mean, xs_std, pc(:,1), 'VariableNames', {'erp' 'xs_mean' 'xs_std' 'pc1'});
        PCA_info.(return_label).(horizon_label).ddm.coeff = coeff;
        PCA_info.(return_label).(horizon_label).ddm.share_var = share_var;
        
        % find principal component of cross-sectional regressions
        cross_section_ts = [];
        for i = 1:length(cross_section)
            cross_section_ts = [cross_section_ts current_spec.(cross_section{i})];
        end
       [pc eigenvectors share_var xs_mean xs_std coeff]  = PCA_unbalanced(cross_section_ts(keep,:),num_of_PC,1e-4,1e5,1,1,1);
        erp_pc.(return_label).(horizon_label).cross_section = timetable(current_spec.Time(keep), ...
            pc(:,1)+nanmean(xs_mean), xs_mean, xs_std, pc(:,1), 'VariableNames', {'erp' 'xs_mean' 'xs_std' 'pc1'});
        PCA_info.(return_label).(horizon_label).cross_section.coeff = coeff;
        PCA_info.(return_label).(horizon_label).cross_section.share_var = share_var;
        
        % find principal component of time-series regressions
        time_series_ts = [];
        for i = 1:length(time_series)
            time_series_ts = [time_series_ts current_spec.(time_series{i})];
        end
        [pc eigenvectors share_var xs_mean xs_std coeff]  = PCA_unbalanced(time_series_ts(keep,:),num_of_PC,1e-4,1e5,1,1,1);
        erp_pc.(return_label).(horizon_label).time_series = timetable(current_spec.Time(keep), ...
            pc(:,1)+nanmean(xs_mean), xs_mean, xs_std, pc(:,1), 'VariableNames', {'erp' 'xs_mean' 'xs_std' 'pc1'});
        PCA_info.(return_label).(horizon_label).time_series.coeff = coeff;
        PCA_info.(return_label).(horizon_label).time_series.share_var = share_var;
        
        % find principal component surveys
        surveys_ts = [];
        for i = 1:length(surveys)
            surveys_ts = [surveys_ts current_spec.(surveys{i})];
        end
        [pc eigenvectors share_var xs_mean xs_std]  = PCA_unbalanced(surveys_ts(keep,:),num_of_PC,1e-4,1e5,1,1,1);
        erp_pc.(return_label).(horizon_label).surveys = timetable(current_spec.Time(keep), ...
            pc(:,1)+nanmean(xs_mean), xs_mean, xs_std, pc(:,1), 'VariableNames', {'erp' 'xs_mean' 'xs_std' 'pc1'});
        PCA_info.(return_label).(horizon_label).surveys.coeff = 1;
        PCA_info.(return_label).(horizon_label).surveys.share_var = 1;
        
        % compute risk free rates
        if return_horizon < 1.5
            rf_ts.(return_label).(horizon_label) = change_horizon(data_ts.RF, data_ts.Time, 1/12, return_horizon, geometric_return, {'RF'});
        else 
            label = ['SVENY0' num2str(return_horizon)];
            rf_ts.(return_label).(horizon_label) = timetable(data_ts.Time, data_ts.(label));
        end
    end
end

%% save
save([output_dir filesep 'erp_principal_components.mat'],'erp_pc','PCA_info')
save([next_dir filesep 'erp_principal_components.mat'],'erp_pc')
save([next_dir2 filesep 'erp_principal_components.mat'],'erp_pc')

%% write csv

for return_horizon = [1/12 1/4 1/2 1 2 3 4 5]
    
    % set up labels to name variables
    switch return_horizon
        case 1/12
            horizon_label = 'one_month_ahead';
            horizon_label_excel = '1m';
            horizon_label_excel2 = 'One month ahead';
        case 1/4
            horizon_label = 'one_quarter_ahead';
            horizon_label_excel = '1q';
            horizon_label_excel2 = 'One quarter ahead';
        case 1/2
            horizon_label = 'six_months_ahead';
            horizon_label_excel = '6m';
            horizon_label_excel2 = 'Six months ahead';
        case 1
            horizon_label = 'one_year_ahead';
            horizon_label_excel = '1yr';
            horizon_label_excel2 = 'One year ahead';
        case 2
            horizon_label = 'two_years_ahead';
            horizon_label_excel = '2yr';
            horizon_label_excel2 = 'Two years ahead';
        case 3
            horizon_label = 'three_years_ahead';
            horizon_label_excel = '3yr';
            horizon_label_excel2 = 'Three years ahead';
        case 4
            horizon_label = 'four_years_ahead';
            horizon_label_excel = '4yr';
            horizon_label_excel2 = 'Four years ahead';
        case 5
            horizon_label = 'five_years_ahead';
            horizon_label_excel = '5yr';
            horizon_label_excel2 = 'Five years ahead';
    end
    
    % name of the excel file
    filename = [excel_output filesep 'ERP_principal_Components_' ...
        horizon_label_excel '.csv'];

    % set up headers
    headers = {horizon_label_excel2, 'historical mean', 'ddm', ...
        'cross sectional regressions', 'time series regressions', ...
        'surveys', 'risk free rate', 'recession', 'zeros'};
    
    % save estimates of ERP found using principal components
    % plus a column with a recession dummy and a column of zeros
    recession_ts = 200*data_ts.recession-100; % +100 means recession, -100 means not recession (useful for excel graphs)
    save_ts = synchronize(timetable(erp_pc.(return_label).(horizon_label).all.Time, erp_pc.(return_label).(horizon_label).all.erp), timetable(data_ts.Time,recession_ts));
    save_ts = synchronize(save_ts, rf_ts.(return_label).(horizon_label));
    model_type = {'all','hist_mean','ddm','cross_section','time_series','surveys'};
    for j=2:length(model_type)
        save_ts = synchronize(save_ts,timetable(erp_pc.(return_label).(horizon_label).(model_type{j}).Time, erp_pc.(return_label).(horizon_label).(model_type{j}).erp));
    end
    save_ts = save_ts(save_ts.Time>=datetime(1960,1,1),:);
    save_ts = save_ts(:, [1 4 5 6 7 8 3 2]);
    models = table2array(save_ts); % put models in the right order and append recession dummy column at the end
    dates = datenum(save_ts.Time);

    export_to_excel = array2table([dates round(models,2)]);
    ind = sum(isnan(table2array(export_to_excel)), 2);
    export_to_excel = export_to_excel(ind < 6,:);
    export_to_excel = table2cell(export_to_excel);
    
    for i = 1:size(export_to_excel,2)
        temp = export_to_excel(:,i);
        temp(cellfun(@isnan,temp)) = {[]};
        export_to_excel(:,i) = temp;
    end

    export_to_excel = cell2table(export_to_excel);
    export_to_excel.Properties.VariableNames = ['date' model_type 'RF' 'recession'];
    
    export_to_excel.date = datestr(export_to_excel.date);
    writetable(export_to_excel, filename);
        
end

for return_horizon = [1/12 1/4 1/2 1 2 3 4 5]
    
    % set up labels to name variables
    switch return_horizon
        case 1/12
            horizon_label = 'one_month_ahead';
            horizon_label_excel = '1m';
        case 1/4
            horizon_label = 'one_quarter_ahead';
            horizon_label_excel = '1q';
        case 1/2
            horizon_label = 'six_months_ahead';
            horizon_label_excel = '6m';
        case 1
            horizon_label = 'one_year_ahead';
            horizon_label_excel = '1yr';
        case 2
            horizon_label = 'two_years_ahead';
            horizon_label_excel = '2yr';
        case 3
            horizon_label = 'three_years_ahead';
            horizon_label_excel = '3yr';
        case 4
            horizon_label = 'four_years_ahead';
            horizon_label_excel = '4yr';
        case 5
            horizon_label = 'five_years_ahead';
            horizon_label_excel = '5yr';
    end
    
      % name of the excel file
        filename = [excel_output filesep 'Dispersion_' horizon_label_excel '.csv'];

        
        % set up headers
        headers = {'XS_Mean','XS_Std_dev','perc_25','perc_75','IQR','Max','Min','RF', 'recession'};
        
        % save estimates of 1-year ahead ERP found using principal components
        % plus a column with a recession dummy and a column of zeros
        xs_mean = timetable(erp_pc.(return_label).(horizon_label).all.Time,erp_pc.(return_label).(horizon_label).all.xs_mean, 'VariableNames', {'xs_mean'});
        xs_std = timetable(erp_pc.(return_label).(horizon_label).all.Time,erp_pc.(return_label).(horizon_label).all.xs_std,'VariableNames', {'xs_std'});
        percentile25 = timetable(erp_ts.(return_label).(horizon_label).Time,prctile(table2array(erp_ts.(return_label).(horizon_label)),25,2),'VariableNames', {'percentile25'});
        percentile75 = timetable(erp_ts.(return_label).(horizon_label).Time,prctile(table2array(erp_ts.(return_label).(horizon_label)),75,2),'VariableNames',{'percentile75'});
        iqr = timetable(erp_ts.(return_label).(horizon_label).Time,table2array(percentile75)-table2array(percentile25),'VariableNames', {'iqr'});
        xs_max = timetable(erp_ts.(return_label).(horizon_label).Time,max(table2array(erp_ts.(return_label).(horizon_label)),[],2),'VariableNames',{'xs_max'});
        xs_min = timetable(erp_ts.(return_label).(horizon_label).Time,min(table2array(erp_ts.(return_label).(horizon_label)),[],2),'VariableNames', {'xs_min'});
        save_ts_fields = {'xs_mean','xs_std','percentile25','percentile75','iqr','xs_max','xs_min'};
        save_ts = synchronize(xs_mean,xs_std,percentile25,percentile75,iqr,xs_max,xs_min);
        save_ts = synchronize(save_ts, rf_ts.(return_label).(horizon_label));
        recession_ts = 200*data_ts.recession-100; % +100 means recession, -100 means not recession (useful for excel graphs)
        save_ts = synchronize(save_ts, timetable(data_ts.Time, recession_ts,'VariableNames', {'recession'}));
        save_ts = save_ts(save_ts.Time>=datetime(1960,1,1),:);
        data = table2array(save_ts); % put models in the right order and append recession dummy column at the end
        dates = datenum(save_ts.Time);
        
        export_to_excel = array2table([dates round(data,2)]);
        ind = sum(isnan(table2array(export_to_excel)), 2);
        export_to_excel = export_to_excel(ind < 6,:);
        export_to_excel = table2cell(export_to_excel);
    
        for i = 1:size(export_to_excel,2)
            temp = export_to_excel(:,i);
            temp(cellfun(@isnan,temp)) = {[]};
            export_to_excel(:,i) = temp;
        end

        export_to_excel = cell2table(export_to_excel);
        export_to_excel.Properties.VariableNames = ['date' headers];
                
        export_to_excel.date = datestr(export_to_excel.date);
        writetable(export_to_excel, filename);
        
end

warning(stateWarning); % restore warning state

clearvars -except root_dir

disp('Finished erp_principal_components.m')