%% create_term_structure

%%% This file constructs the term structure of the ERP.

stateWarning = warning('off','all');

%% set directories

% main directories
input_dir = [root_dir filesep 'Temp' filesep 'Create term structure of ERP' filesep ...
    'Input'];
input_dir2 = [root_dir filesep 'Temp' filesep 'Create ERP measures' filesep 'Input'];
output_dir = [root_dir filesep 'Temp' filesep 'Create term structure of ERP' filesep ...
    'Output'];
next_dir = [root_dir filesep 'Temp' filesep 'Test for predictability' filesep 'Input'];
excel_output = [root_dir filesep 'Output'];

% make directories if they don't exist
dir_out = {output_dir,next_dir};
for i=1:length(dir_out)
    if ~exist(dir_out{i},'dir')
        mkdir(dir_out{i});
    end
end

clear i 

%% load data
load([input_dir filesep 'erp_principal_components.mat'])
load([input_dir2 filesep 'variables.mat'])

%% create term structure of ERP

for geometric_return = [0 1]
    
    % set up labels to name variables
    if geometric_return==1
        return_label = 'geometric';
    else
        return_label = 'arithmetic';
    end
    
    % initialize variable with first specification
    term_structure.(return_label) = timetable(erp_pc.(return_label).('one_month_ahead').all.Time, ...
        erp_pc.(return_label).('one_month_ahead').all.erp);
    
    % loop through other specifications
    for return_horizon = [1/4 1/2 1 2 3 4 5]
        
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
        
        term_structure.(return_label) = synchronize(term_structure.(return_label),...
            timetable(erp_pc.(return_label).(horizon_label).all.Time, erp_pc.(return_label).(horizon_label).all.erp));
    end
    term_structure.(return_label).Properties.VariableNames = {'erp_one_month_ahead',...
            'erp_one_quarter_ahead', 'erp_six_months_ahead', 'erp_one_year_ahead', ...
            'erp_two_years_ahead', 'erp_three_years_ahead', 'erp_four_years_ahead', ...
            'erp_five_years_ahead'};
end

%% create "compact" term structure for periods in which ERP is really high or really low

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
        
        % find highest and lowest ERP dates
        num_points = 6; % pick how many high and low points are selected
        one_erp = erp_pc.(return_label).(horizon_label).all.erp;
        one_year_erp_dates = erp_pc.(return_label).(horizon_label).all.Time;
        [~,minmax_index] = sort(one_erp,1,'ascend');
        dates = datenum([one_year_erp_dates(minmax_index(1:num_points)); one_year_erp_dates(minmax_index(end-num_points+1:end))]);
        dates = [dates; datenum(1982,12,31); datenum(1987,9,30)]; % add december 1982 and september 1987 because there was an elevated ERP compared to contiguous periods
        
        [~,dates_index]= intersect(datenum(term_structure.(return_label).Time),dates);
        
        % create "compact" term structure with high/low ERP point, plus the
        % latest point
       val = strcat('erp_',horizon_label);
        compact_term_structure.(return_label).(horizon_label) = timetable(term_structure.(return_label).Time(dates_index), term_structure.(return_label).(val)(dates_index));
        
    end
end

%% create "counterfactual" term structure

% construct yield curve
yield_curve = timetable(data_ts.Time,data_ts.('RF'));
for return_horizon = [1/12 1/4 1/2 1 2 3 4 5]
    
    % set up labels to name variables
    switch return_horizon
        case 1/12
            horizon_label = {'one_month_ahead'};
            data_label = 'RF'; % Fama-French one month risk-free rate
        case 1/4
            horizon_label = {'one_quarter_ahead'};
            data_label ='TB3MS'; % FRED three month yield
        case 1/2
            horizon_label = {'six_months_ahead'};
            data_label = 'TB6MS'; % FRED six month yield
        case 1
            horizon_label = {'one_year_ahead'};
            data_label = 'SVENY01'; % Fed Board 1 year yield
        case 2
            horizon_label = {'two_years_ahead'};
            data_label = 'SVENY02'; % Fed Board 2 year yield
        case 3
            horizon_label = {'three_years_ahead'};
            data_label = 'SVENY03'; % Fed Board 3 year yield
        case 4
            horizon_label = {'four_years_ahead'};
            data_label = 'SVENY04'; % Fed Board 4 year yield
        case 5
            horizon_label = {'five_years_ahead'};
            data_label = 'SVENY05'; % Fed Board 5 year yield
    end
    
    yield_curve = synchronize(yield_curve,timetable(data_ts.Time,data_ts.(data_label),'VariableNames', horizon_label));
     
end
yield_curve = yield_curve(:,2:end);

% eliminate the last observations for yield_curve and
% compact_term_structure if they are NaN
any_yield_nan = any(isnan(table2array(yield_curve)),2); % finds if there is at least one NaN in each row
number_to_remove = find(diff(cumsum(flipud(any_yield_nan)))==0,1,'first'); % finds the number of NaN at the end of the time-series
last_date_keep = yield_curve.Time(end-number_to_remove); % find last day to keep

% compute mean of yield curve over time
header_yields ={'one_month_ahead';'one_quarter_ahead';'six_months_ahead'; 'one_year_ahead';'two_years_ahead';'three_years_ahead'};
yield_curve_mat = [];
for i = 1:length(header_yields)
    yield_curve_mat = [yield_curve_mat yield_curve.(header_yields{i})]; % convert financial time series to matrix and order maturities by row
end
mean_yield_curve = mean(yield_curve_mat(~any_yield_nan,:)); % compute mean of the yield curve

% compute counterfactual ERP term structure
headers_erp ={'erp_one_month_ahead';'erp_one_quarter_ahead';'erp_six_months_ahead'; 'erp_one_year_ahead';'erp_two_years_ahead';'erp_three_years_ahead'};

for geometric_return = [0 1]
    
    % set up labels to name variables
    if geometric_return==1
        return_label = 'geometric';
    else
        return_label = 'arithmetic';
    end
    
    % compute counterfactual ERP term structure if yield is at its mean level
    erp_term_structure_today = term_structure.(return_label)(end-number_to_remove,:);
    %compact_term_structure.(return_label).('one_month_ahead')(end); % horizon is irrelevant, last point is always last point in the time series
    yield_curve_today = yield_curve(yield_curve.Time == erp_term_structure_today.Time,:);
    yield_curve_diff = [];
    for i = 1:length(header_yields)
        yield_curve_diff = [yield_curve_diff yield_curve_today.((header_yields{i}))-mean_yield_curve(i)];
    end
    for i =1:length(headers_erp)
        counterfactual_term_structure.(return_label) = erp_term_structure_today.(headers_erp{i})+yield_curve_diff(i);
        today_term_structure_mat.(return_label) = erp_term_structure_today.(headers_erp{i});
    end
    
end

% compute counterfactual ERP term structure if yield is at its mean level,
% taking into account historical relation between erp and 10 year yield

for geometric_return = [0 1]
    
    % set up labels to name variables
    if geometric_return==1
        return_label = 'geometric';
    else
        return_label = 'arithmetic';
    end
    
    regression_ts = synchronize(timetable(data_ts.Time, data_ts.SVENY05, data_ts.SVENY10, table2array(yield_curve)), term_structure.(return_label));
    regression_ts = regression_ts(~any(isnan(table2array(regression_ts)),2),:); % eliminate missing observations
    Y_all = table2array(regression_ts(:,4:end));
    X_all = [regression_ts.Var3(:,4) regression_ts.Var1 regression_ts.Var2];
    counterfactual_date = regression_ts.Time(end);
    yield_curve_diff2 = -[yield_curve_diff(4) regression_ts.Var1(end)-nanmean(regression_ts.Var1) regression_ts.Var2(end)-nanmean(regression_ts.Var2)];
    
    j=1;
    for return_horizon = [1/12 1/4 1/2 1 2 3 4 5]
        
        % run a regression of ERP on bond yields
        window = size(Y_all,1); % run regression over last window months
        Y = diff(Y_all(end-window+1:end,j));
        X = [ones(size(Y)) diff(X_all(end-window+1:end,:))];      
        [b,bint,r,rint,stats] = regress(Y,X);
        counterfactual_regression.(return_label)(1,j) = Y(end)+[1 yield_curve_diff2]*b;
        
        j = j+1;
        
    end
end

%% save
save([output_dir filesep 'term_structure_erp.mat'],'term_structure')
save([next_dir filesep 'term_structure_erp.mat'],'term_structure')

%% write csv

for geometric_return = [0 1]
    
    % set up labels to name variables
    if geometric_return==1
        return_label = 'geometric';
    else
        return_label = 'arithmetic';
    end
    
    % name of the excel file
    filename = [excel_output filesep 'Term_Structure_' return_label '.csv'];
    
    % headers 
    headers ={'erp_one_month_ahead';'erp_one_quarter_ahead'; ...
        'erp_six_months_ahead'; 'erp_one_year_ahead';...
        'erp_two_years_ahead';'erp_three_years_ahead';...
        'erp_four_years_ahead'; 'erp_five_years_ahead'};
    
    % saving values
    save_ts = term_structure.(return_label);
    recession_ts = 200*data_ts.recession-100; % +100 means recession, -100 means not recession (useful for excel graphs)
    save_ts = synchronize(save_ts, timetable(data_ts.Time, recession_ts, 'VariableNames', {'recession'}));
    save_ts = synchronize(save_ts, timetable(save_ts.Time, zeros(size(save_ts.Time)),'VariableNames', {'zeros'}));
    save_ts = save_ts(save_ts.Time>=datetime(1959,12,31),:);
    mean_term_structure.(return_label) = nan(1,size(save_ts,2));
    mean_term_structure.(return_label)(:,1:8) = nanmean(table2array(save_ts(:,1:8)));
    models = [table2array(save_ts); mean_term_structure.(return_label); zeros(1,size(save_ts,2))]; 
    models = round(models,2);
    dates = datenum(save_ts.Time);
    
    rownames = cellstr(datestr(save_ts.Time));
    rownames{length(rownames)+1,1} = 'mean';
    rownames{length(rownames)+1,1} = 'zeros';
    %write csv
    export_to_excel = array2table([round(models,2)]);
    export_to_excel.Properties.VariableNames = [headers' 'recession' 'zeros'];
    export_to_excel.Properties.RowNames = rownames;
       
    ind = sum(isnan(table2array(export_to_excel(:,2:end))),2);
    export_to_excel = export_to_excel(ind < 6,:);
    
    writetable(export_to_excel, filename);
    
end

for geometric_return = [0 1]
    
    % set up labels to name variables
    if geometric_return==1
        return_label = 'geometric';
    else
        return_label = 'arithmetic';
    end
    
    % name of the excel file
    filename = [excel_output filesep 'Compact_Term_Structure_' return_label '.csv'];
    
 
    % headers2 ={'One month','One quarter','Six months', 'One year','Two years','Three years','Four years','Five years'};
    headers2 ={'One month','One quarter','Six months', 'One year','Two years','Three years'};
        
    % save term structure of ERP
    horizon_label = 'one_month_ahead';
    headers ={'one_month_ahead';'one_quarter_ahead';'six_months_ahead'; 'one_year_ahead';'two_years_ahead';'three_years_ahead';...
        'four_years_ahead';'five_years_ahead'};
    save_ts = compact_term_structure.(return_label).(char(headers(1)));
    for i = 2:length(headers)
        save_ts = synchronize(save_ts, compact_term_structure.(return_label).(char(headers(i))));
    end
    save_ts.Properties.VariableNames = headers; 
    mean_term_structure.(return_label) = nan(1,size(save_ts,2));
    mean_term_structure.(return_label)(:,1:8) = nanmean(table2array(save_ts(:,1:8)));
    models = [table2array(save_ts); mean_term_structure.(return_label); zeros(1,size(save_ts,2))]; % add two rows at the end, one with the average term structure and one with zeros
    rownames = cellstr(datestr(save_ts.Time));
    rownames{length(rownames)+1,1} = 'mean';
    rownames{length(rownames)+1,1} = 'zeros';
    
    % write csv
    export_to_excel = array2table(round(models,2));
    export_to_excel.Properties.VariableNames = headers';
    export_to_excel.Properties.RowNames = rownames;
    
    ind = sum(isnan(table2array(export_to_excel(:,2:end))),2);
    export_to_excel = export_to_excel(ind < 3,:);
    
    writetable(export_to_excel, filename, 'WriteRowNames', true);
  
end

for geometric_return = [0 1]
    
    % set up labels to name variables
    if geometric_return==1
        return_label = 'geometric';
    else
        return_label = 'arithmetic';
    end
    
    % name of the excel file
    filename = [excel_output filesep 'Earnings_' return_label '.csv'];
    
    % set up headers
    headers1 = {'date', 'Real_earnings_growth','YoY_change_in_expectations_of_EPS','zeros'};
    
    % save term structure of ERP
    dates_export = datenum(data_ts.Time(end-window+1:end));
    real_E_growth = data_ts.real_E_growth(end-window+1:end);
    real_yoy_expectations = data_ts.eps_1yr_ahead(end-window+1:end)-data_ts.eps_1yr_ahead(end-window+1-12:end-12);
    data = [dates_export round(real_E_growth,2) round(real_yoy_expectations,2) zeros(size(real_yoy_expectations,1),1)];
    
    % write csv
    export_to_excel = array2table(round(data,2));
    ind = sum(isnan(table2array(export_to_excel(:,2:end))),2);
    export_to_excel = export_to_excel(ind < 2,:);
    export_to_excel.Properties.VariableNames = headers1;
    export_to_excel.date = datestr(export_to_excel.date);
    export_to_excel = table2cell(export_to_excel);
    
    for i = 2:size(export_to_excel,2)
        temp = export_to_excel(:,i);
        temp(cellfun(@isnan,temp)) = {[]};
        export_to_excel(:,i) = temp;
    end
    
    export_to_excel = cell2table(export_to_excel);
    export_to_excel.Properties.VariableNames = headers1;
    
    
    writetable(export_to_excel, filename); 
    
end

warning(stateWarning); % restore warning state

clearvars -except root_dir

disp('Finished create_term_structure.m')

