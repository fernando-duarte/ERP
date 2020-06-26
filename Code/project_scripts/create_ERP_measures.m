%% create_ERP_measures

%%% This file constructs several measures of the ERP at 
%%% different horizons using different models.

stateWarning = warning('off','all');

model_ratio = .5;

%% set directories

% main directories
input_dir = [root_dir filesep 'Temp' filesep 'Create ERP measures' filesep 'Input'];
output_dir = [root_dir filesep 'Temp' filesep 'Create ERP measures' filesep 'Output'];
next_dir = [root_dir filesep 'Temp' filesep 'Principal components' filesep 'Input'];
next_dir2 = [root_dir filesep 'Temp' filesep 'Test for predictability' filesep 'Input'];
Matlab_subfunctions = [root_dir filesep 'Code' filesep 'Matlab subfunctions'];
factor_models = [Matlab_subfunctions filesep 'Factor models' filesep];
lightspeed = [Matlab_subfunctions filesep 'lightspeed' filesep];
excel_output = [root_dir filesep 'Output'];

% make directories if they don't exist
dir_out = {output_dir,next_dir,next_dir2};
for i=1:length(dir_out)
    if ~exist(dir_out{i},'dir')
        mkdir(dir_out{i});
    end
end
cd(root_dir)

%adding paths for functions
addpath(Matlab_subfunctions)
addpath(factor_models)
addpath(lightspeed)
clear i

%% load data
load([input_dir filesep 'variables.mat'])

%% construct ERP measures

% loop through different horizons (monthly, quarterly, semi-annual, annual,
% 2, 3, 4 and 5 years
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
    
    % use geometric or arithmetic returns
    for geometric_return = [0 1]
        
        % set up labels to name variables
        if geometric_return==1
            return_label = 'geometric';
        else
            return_label = 'arithmetic';
        end
        
        % counter for model number
        model_no = 1;
        
        % create placeholders to keep track of frequency of models
        monthly = {};
        quarterly = {};
        annual = {};
        
        % create placeholders to keep track of model description
        model_description = {};
        
        %% Historical mean models
        hist_mean = {}; % will be used to keep track of what models belong to this category
        
        % model: historical mean as far back as data allows
        summary = 'historical mean as far back as data allows';
        realized_erp_ts = change_horizon(data_ts.realized_erp,data_ts.Time, 1/12,return_horizon,geometric_return,{'realized_erp'}); % change horizon by computing geometric or arithmetic multi-period returns
        realized_erp = lag(realized_erp_ts); % shift time-series so that you only use past information when taking the mean
        keep = ~isnan(realized_erp.realized_erp);
        model = nancumsum(realized_erp.realized_erp,1,2)./nancumsum(keep,1,2); % historical mean of realized ERP
        erp_ts.(return_label).(horizon_label) = timetable(realized_erp_ts.Time, model);
        
        % record characteristics
        hist_mean = [hist_mean ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        
        monthly = [monthly ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        % model: historical mean using last five years
        summary = 'historical mean using last five years';
        realized_erp_clean = realized_erp(keep,:);
        window_size = 5*12;
        model = nan(size(realized_erp));
        model(keep) = filter(ones(1,window_size)/window_size,1,realized_erp_clean.realized_erp); % average over last 60 months
        erp_ts.(return_label).(horizon_label) = synchronize(erp_ts.(return_label).(horizon_label), timetable(realized_erp_ts.Time, model));
       
        % record characteristics
        hist_mean = [hist_mean ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        monthly = [monthly ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        %% Dividend discount mean models
        ddm = {}; % keep track of what models belong to this category
        
        % model: E/P minus nominal 10yr yield
        summary = 'E/P minus nominal 10yr yield';
        model = data_ts.ep-data_ts.RLONG;
        erp_ts.(return_label).(horizon_label) = join(erp_ts.(return_label).(horizon_label), timetable(data_ts.Time, model));
       
        % record characteristics
        ddm = [ddm ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        monthly = [monthly ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        % model: Shiller: 1/CAPE minus nominal 10yr yield
        summary = 'Shiller: 1/CAPE minus nominal 10yr yield';
        model = (100./data_ts.CAPE)-data_ts.RLONG;
        erp_ts.(return_label).(horizon_label) =  join(erp_ts.(return_label).(horizon_label), timetable(data_ts.Time, model));
       
        % record characteristics
        ddm = [ddm ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        monthly = [monthly ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        % model: E/P minus real 10yr yield
        summary = 'E/P minus real 10yr yield';
        model = data_ts.ep-data_ts.FII10;
        erp_ts.(return_label).(horizon_label) =  join(erp_ts.(return_label).(horizon_label), timetable(data_ts.Time, model));
       
        % record characteristics
        ddm = [ddm ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        monthly = [monthly ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        % model: expected E/P minus real 10yr yield
        summary = 'expected E/P minus real 10yr yield';
        model = 100*data_ts.eps_1yr_ahead./data_ts.sp500_price-data_ts.FII10;
        erp_ts.(return_label).(horizon_label) =  join(erp_ts.(return_label).(horizon_label), timetable(data_ts.Time, model));
       
        % record characteristics
        ddm = [ddm ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        monthly = [monthly ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        % model: expected E/P minus nominal 10yr yield
        summary = 'expected E/P minus nominal 10yr yield';
        model = 100*data_ts.eps_1yr_ahead./data_ts.sp500_price-data_ts.RLONG;
        erp_ts.(return_label).(horizon_label) =  join(erp_ts.(return_label).(horizon_label), timetable(data_ts.Time, model));
       
        % record characteristics
        ddm = [ddm ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        monthly = [monthly ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        % model: Two-stage DDM
        summary = 'Two-stage DDM';
        
        % estimate short-run (5 years) growth rate of earnings
        lagged = lag(data_ts,-1);
        E_growth_full = 100*12*(data_ts.E(2:end)./data_ts.E(1:end-1) - 1); % growth rate of earnings
        window_size = 5*12;
        E_growth_5yr = filter(ones(1,window_size)/window_size,1,E_growth_full); % average growth rate over the last 5 years
        earnings_growth = E_growth_5yr(2:end);
        lag_earnings_growth = E_growth_5yr(1:end-1);
        ep = data_ts.ep(2:end); % earnings-price ratio
        lag_ep = ep(1:end-1);% lagged earnings-price ratio
        keep = ~isnan(earnings_growth) & ~isnan(lag_earnings_growth) & ~isnan(lag_ep);
        X = [ones(size(earnings_growth(keep))) lag_earnings_growth(keep) lag_ep(keep)];
        Y = earnings_growth(keep);
        g_SR = nan(size(E_growth_full));
        g_SR(keep) = X*regress(Y,X);
        g_SR(end+1) = NaN;
        % estimate for long-run growth rate of earnings
        g_LR = 2.2; % by assumption 2.2%
        
        % payout ratio
        payout = 0.4; % by assumption 40%
        
        % construct ERP from 2-stage DDM model
        model=(100+g_LR+5*(g_SR-g_LR)).*(data_ts.E*payout./data_ts.P)-data_ts.RLONG+g_LR;
        erp_ts.(return_label).(horizon_label) =  join(erp_ts.(return_label).(horizon_label), timetable(data_ts.Time, model));
       
        % record characteristics
        ddm = [ddm ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        monthly = [monthly ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        % model: 6 stage DDM by Damodaran
        summary = '6 stage DDM by Damodaran';
        model = data_ts.damodaran_ddm;
        erp_ts.(return_label).(horizon_label) = join(erp_ts.(return_label).(horizon_label), timetable(data_ts.Time, model));
       
        % record characteristics
        ddm = [ddm ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        annual = [annual ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        % model: 6 stage free cash flow DDM by Damodaran
        summary = '6 stage free cash flow DDM by Damodaran';
        model = data_ts.damodaran_fcfe;
        erp_ts.(return_label).(horizon_label) = join(erp_ts.(return_label).(horizon_label), timetable(data_ts.Time, model));
       
        % record characteristics
        ddm = [ddm ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        annual = [annual ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        %% Cross-sectional regressions
        cross_section = {}; % keep track of what models belong to this category
        
        % model: fama-french
        summary = 'fama-french';
        portfolios = cellfun(@(x) ['ff_p' num2str(x)],mat2cell((1:25)',ones(25,1),1),'UniformOutput',0); % pick what portfolios to use for the cross-sectional regression
        variables = [];
        for i = 1:length(portfolios)
            variables = [variables data_ts.(char(portfolios{i}))];
        end
        names = {'MktRF'; 'HML'; 'SMB'};
        names = vertcat(names,portfolios);
        ff_ts = change_horizon([data_ts.MktRF,data_ts.HML,data_ts.SMB,variables], data_ts.Time, 1/12,return_horizon,geometric_return,names); % change horizon by computing geometric or arithmetic multi-period returns
        X = []; Y = [];
        for i = 1:3
            X = [X ff_ts.(char(names{i}))];
        end
        for i = 1:length(portfolios)
            Y = [Y ff_ts.(char(portfolios{i}))];
        end
        keep = ~any(isnan(Y),2) | ~any(isnan(X),2);
        window = 60; % 5 year window
        [alpha, beta, stats] = rolling_regress(Y(keep,:),X(keep,:),'one_sided_uniform',(window-1)/size(X(keep,:),1));
        [lambda, errors, lambda_t]= fama_macbeth(Y(keep,:),beta);
        model = nancumsum(lambda_t(:,1),1,2)./nancumsum(~isnan(lambda_t(:,1)),1, 2);
        erp_ts.(return_label).(horizon_label) = join(erp_ts.(return_label).(horizon_label), timetable(ff_ts.Time, model));
       
        % record characteristics
        cross_section = [cross_section ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        monthly = [monthly ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        clear portfolios
        
        % model: fama-french and momentum
        summary = 'fama-french and momentum';
        portfolios = [cellfun(@(x) ['ff_p' num2str(x)],mat2cell((1:25)',ones(25,1),1),'UniformOutput',0);cellfun(@(x) ['mom_p' num2str(x)],mat2cell((1:10)',ones(10,1),1),'UniformOutput',0)]; % pick what portfolios to use for the cross-sectional regression
        variables = [];
        for i = 1:length(portfolios)
            variables = [variables data_ts.(char(portfolios{i}))];
        end
        names = {'MktRF'; 'HML'; 'SMB';'mom'};
        names = vertcat(names,portfolios);
        ff_ts = change_horizon([data_ts.MktRF,data_ts.HML,data_ts.SMB,data_ts.mom,variables],data_ts.Time,1/12,return_horizon,geometric_return,names); % change horizon by computing geometric or arithmetic multi-period returns
        X = []; Y = [];
        for i = 1:4
            X = [X ff_ts.(char(names{i}))];
        end
        for i = 1:length(portfolios)
            Y = [Y ff_ts.(char(portfolios{i}))];
        end
        keep = ~any(isnan(Y),2) & ~any(isnan(X),2);
        ff_ts = ff_ts(keep,:);
        window = 60; % 5 year window
        [alpha, beta, stats] = rolling_regress(Y(keep,:),X(keep,:),'one_sided_uniform',(window-1)/size(X(keep,:),1));
        [lambda, errors, lambda_t]= fama_macbeth(Y(keep,:),beta);
        model = nancumsum(lambda_t(:,1),1,2)./nancumsum(~isnan(lambda_t(:,1)),1, 2);
        erp_ts.(return_label).(horizon_label) =  synchronize(erp_ts.(return_label).(horizon_label), timetable(ff_ts.Time, model));
       
        clear portfolios
        % record characteristics
        cross_section = [cross_section ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        monthly = [monthly ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        % model: fama-french, momentum and inflation
        summary = 'fama-french, momentum and inflation';
        portfolios = [cellfun(@(x) ['ff_p' num2str(x)],mat2cell((1:25)',ones(25,1),1),'UniformOutput',0);cellfun(@(x) ['mom_p' num2str(x)],mat2cell((1:10)',ones(10,1),1),'UniformOutput',0)]; % pick what portfolios to use for the cross-sectional regression
        variables = [];
        for i = 1:length(portfolios)
            variables = [variables data_ts.(char(portfolios{i}))];
        end
        names = {'MktRF'; 'HML'; 'SMB';'mom';'inflation'};
        names = vertcat(names,portfolios);
        ff_ts = change_horizon([data_ts.MktRF,data_ts.HML,data_ts.SMB,data_ts.mom,data_ts.inflation variables],data_ts.Time,1/12,return_horizon,geometric_return,names); % change horizon by computing geometric or arithmetic multi-period returns
         X = []; Y = [];
        for i = 1:5
            X = [X ff_ts.(char(names{i}))];
        end
        for i = 1:length(portfolios)
            Y = [Y ff_ts.(char(portfolios{i}))];
        end
        keep = ~any(isnan(Y),2) & ~any(isnan(X),2);
        ff_ts = ff_ts(keep,:);
        window = 60; % 5 year window
        [alpha beta stats] = rolling_regress(Y(keep,:),X(keep,:),'one_sided_uniform',(window-1)/size(X(keep,:),1));
        [lambda errors lambda_t]= fama_macbeth(Y(keep,:),beta);
        model = nancumsum(lambda_t(:,1),1,2)./nancumsum(~isnan(lambda_t(:,1)),1, 2);
        erp_ts.(return_label).(horizon_label) =  synchronize(erp_ts.(return_label).(horizon_label), timetable(ff_ts.Time, model));
       
        clear c portfolios
        % record characteristics
        cross_section = [cross_section ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        monthly = [monthly ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        % model: Adrian, Crump, and Moench (2012)
        summary = 'Adrian, Crump, and Moench';
        model = data_ts.NYFed_ERP;
        model = change_horizon(model,data_ts.Time, 1/12,return_horizon,geometric_return, {'NYFed_ERP'}); % change horizon by computing geometric or arithmetic multi-period returns
        erp_ts.(return_label).(horizon_label) = synchronize(erp_ts.(return_label).(horizon_label), model);
        
        % record characteristics
        cross_section = [cross_section ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        monthly = [monthly ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        %% Time-series regressions
        time_series = {}; % keep track of what models belong to this category
        
        % model: Predictor is D/P
        summary = 'Time-series predictor is D/P';
        realized_erp_ts = change_horizon(data_ts.realized_erp,data_ts.Time, 1/12,return_horizon,geometric_return, {'realized_erp'}); % change horizon by computing geometric or arithmetic multi-period returns
        Xvar = timetable(data_ts.Time,data_ts.dp);
        Xvar.Properties.VariableNames = {'dp'};
        XY_ts = synchronize(realized_erp_ts,Xvar); % make sure dates for dependent and independent variables agree
        X = XY_ts.dp; % dividend price ratio is right hand side variable
        Y = XY_ts.realized_erp; % realized erp is left hand side variable
        dates = XY_ts.Time;
        [model, dates] = predictive_regression(Y,X,dates,12*return_horizon);
        erp_ts.(return_label).(horizon_label) =  synchronize(erp_ts.(return_label).(horizon_label), timetable(dates, model.oos));
        
        % record characteristics
        time_series = [time_series ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        monthly = [monthly ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        % model: Best predictor of Goyal and Welch
        summary = 'Best predictor in Goyal and Welch';
        predictors = {'dp','ep','CAPE','bm','roe','term_spread','BAA_AAA','share_equity_issuance','cay','inflation'}; 
        
        % run regressions
        clear predicted_x dates r2 coefficients residual bint X
        for variable = predictors
            Xvar = timetable(data_ts.Time, data_ts.(char(variable)));
            Xvar.Properties.VariableNames = variable;
            XY_ts = synchronize(realized_erp_ts,Xvar);
            X.(variable{:}) = XY_ts.(char(variable)); % said variable is right hand side variable
            Y = XY_ts.realized_erp; % realized erp is left hand side variable
            [predicted_x.(variable{:}), dates.(variable{:}), r2.(variable{:}), coefficients, residual.(variable{:}), bint.(variable{:})] = predictive_regression(Y,X.(variable{:}),XY_ts.Time,12*return_horizon);
            clear Xvar XY_ts X Y 
        end
        
        % find which was the best model in real time (using oos R^2)
        cellpredictedx = struct2cell(predicted_x); 
        cellpredictedx_oos = cellfun(@(x) x.oos,cellpredictedx,'UniformOutput',0);
        celldates = struct2cell(dates);
        cellr2 = struct2cell(r2); 
        cellr2_oos = cellfun(@(x) x.oos,cellr2,'UniformOutput',0);
        cell_ts_pred = cellfun(@(ts_date, ts_erp, ts_name) timetable(ts_date,ts_erp),celldates, cellpredictedx_oos,predictors','UniformOutput',0);
        cell_ts_r2 = cellfun(@(ts_date, ts_r2, ts_name) timetable(ts_date,ts_r2),celldates, cellr2_oos,predictors','UniformOutput',0);
        
        predictor_ts = synchronize(cell_ts_pred{:});
        r2_ts = synchronize(cell_ts_r2{:});
        [~,keep]=max(table2array(r2_ts),[],2);
        predictor_mat = table2array(predictor_ts);
        best_predictors = predictors(keep);
        
        model = predictor_mat(sub2ind(size(predictor_mat),(1:size(predictor_mat,1))',keep));
        erp_ts.(return_label).(horizon_label) =  synchronize(erp_ts.(return_label).(horizon_label), timetable(predictor_ts.ts_date, model));
        
        % record characteristics
        time_series = [time_series ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        monthly = [monthly ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        % model: Best predictor of Campbell and Thompson
        summary = 'Best predictor in Campbell and Thompson';
        predictors = {'dp','ep','CAPE','bm','roe','term_spread','BAA_AAA','inflation'};
        
        % run regressions
        clear predicted_x dates r2 coefficients residual bint
        for variable = predictors
               Xvar = timetable(data_ts.Time, data_ts.(char(variable)));
            Xvar.Properties.VariableNames = variable;
            XY_ts = synchronize(realized_erp_ts,Xvar);
            X.(variable{:}) = XY_ts.(char(variable)); % said variable is right hand side variable
            Y = XY_ts.realized_erp; % realized erp is left hand side variable
            [predicted_x.(variable{:}), dates.(variable{:}), r2.(variable{:}), coefficients, residual.(variable{:}), bint.(variable{:})] = predictive_regression(Y,X.(variable{:}),XY_ts.Time,12*return_horizon);
       end
        
        % find which was the best model in real time (using oos R^2)
        celldates = struct2cell(dates);
        
        cellr2 = struct2cell(r2);
        cellr2_oos = cellfun(@(x) x.oos,cellr2,'UniformOutput',0);
        cellr2_oos_slope = cellfun(@(x) x.oos_slope,cellr2,'UniformOutput',0);
        cellr2_oos_nonneg = cellfun(@(x) x.oos_nonneg,cellr2,'UniformOutput',0);
        cellr2_oos_both = cellfun(@(x) x.oos_both,cellr2,'UniformOutput',0);
        cell_ts_r2  = cellfun(@(ts_date, r2, ts_name) timetable(ts_date,r2),celldates, cellr2_oos,predictors','UniformOutput',0);
        cell_ts_r2_slope = cellfun(@(ts_date, slope, ts_name) timetable(ts_date,slope),celldates, cellr2_oos_slope,predictors','UniformOutput',0);
        cell_ts_r2_nonneg = cellfun(@(ts_date, nonneg, ts_name) timetable(ts_date,nonneg),celldates, cellr2_oos_nonneg,predictors','UniformOutput',0);
        cell_ts_r2_both = cellfun(@(ts_date, both, ts_name) timetable(ts_date,both),celldates, cellr2_oos_both,predictors','UniformOutput',0);
        r2_ts = cellfun(@(date,r2,slope,nonneg,both,name) timetable(date,r2,slope,nonneg,both),celldates, cell_ts_r2, cell_ts_r2_slope , cell_ts_r2_nonneg , cell_ts_r2_both, predictors','UniformOutput', 0); 
       
        for i = 1:length(r2_ts)
            r2_ts{i,1} = timetable(r2_ts{i,1}.date, r2_ts{i,1}.r2(:,1).r2, r2_ts{i,1}.slope(:,1).slope, r2_ts{i,1}.nonneg(:,1).nonneg, r2_ts{i,1}.both(:,1).both, 'VariableNames', {'r2', 'slope', 'nonneg', 'both'});
        end
        
        cellpredictedx = struct2cell(predicted_x);
        cellpredictedx_oos = cellfun(@(x) x.oos,cellpredictedx,'UniformOutput',0);
        cellpredictedx_oos_slope = cellfun(@(x) x.oos_slope,cellpredictedx,'UniformOutput',0);
        cellpredictedx_oos_nonneg = cellfun(@(x) x.oos_nonneg,cellpredictedx,'UniformOutput',0);
        cellpredictedx_oos_both = cellfun(@(x) x.oos_both,cellpredictedx,'UniformOutput',0);
        cell_ts_pred = cellfun(@(date,r2,slope,nonneg,both,name) timetable(date,r2,slope,nonneg,both),celldates, cellpredictedx_oos,cellpredictedx_oos_slope,cellpredictedx_oos_nonneg,cellpredictedx_oos_both,predictors','UniformOutput', 0); 
        
        predictor_ts = synchronize(cell_ts_pred{:});
        r2_ts = synchronize(r2_ts{:});
        
        [~,keep]=max(table2array(r2_ts),[],2);
        predictor_mat = table2array(predictor_ts);
        
        model = predictor_mat(sub2ind(size(predictor_mat),(1:size(predictor_mat,1))',keep));
        erp_ts.(return_label).(horizon_label) =  synchronize(erp_ts.(return_label).(horizon_label), timetable(predictor_ts.date, model));
        
        % record characteristics
        time_series = [time_series ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        monthly = [monthly ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        % model: Fama and French (2002)
        summary = 'Best predictor in Fama French 2002';
        ff_erp1 = data_ts.dividend_yield+data_ts.real_D_growth-data_ts.RLONG;
        ff_erp2 = data_ts.dividend_yield+data_ts.real_E_growth-data_ts.RLONG;
        ff_erp3 = data_ts.dividend_yield+data_ts.real_P_growth-data_ts.RLONG;
        
        predictor1 = nancumsum(ff_erp1,1, 2)./nancumsum(~isnan(ff_erp1),1, 2);
        predictor2 = nancumsum(ff_erp2,1, 2)./nancumsum(~isnan(ff_erp2),1, 2);
        predictor3 = nancumsum(ff_erp3,1, 2)./nancumsum(~isnan(ff_erp3),1, 2);
        names = {'predictor1','predictor2','predictor3'};
    
        ff_ts = change_horizon([predictor1 predictor2 predictor3], data_ts.Time ,1/12,return_horizon,geometric_return, names); % change horizon by computing geometric or arithmetic multi-period returns
        
        predictors = {'predictor1','predictor2','predictor3'};
        
        % run regressions
       
        clear predicted_x dates r2 coefficients residual bint
        for variable = predictors
            Xvar = timetable(ff_ts.Time, ff_ts.predictor1);
            Xvar.Properties.VariableNames = variable;
            XY_ts = synchronize(realized_erp_ts,Xvar);
            X.(variable{:}) = XY_ts.(char(variable)); % said variable is right hand side variable
            Y = XY_ts.realized_erp; % realized erp is left hand side variable
            [predicted_x.(variable{:}), dates.(variable{:}), r2.(variable{:}), coefficients, residual.(variable{:}), bint.(variable{:})] = predictive_regression(Y,X.(variable{:}),XY_ts.Time,12*return_horizon);
        end
        
        % find which was the best model in real time (using oos R^2)
        celldates = struct2cell(dates);
        
        cellr2 = struct2cell(r2);
        cellr2_oos = cellfun(@(x) x.oos,cellr2,'UniformOutput',0);
        cellr2_oos_slope = cellfun(@(x) x.oos_slope,cellr2,'UniformOutput',0);
        cellr2_oos_nonneg = cellfun(@(x) x.oos_nonneg,cellr2,'UniformOutput',0);
        cellr2_oos_both = cellfun(@(x) x.oos_both,cellr2,'UniformOutput',0);
        cell_ts_r2  = cellfun(@(ts_date, r2, ts_name) timetable(ts_date,r2),celldates, cellr2_oos,predictors','UniformOutput',0);
        cell_ts_r2_slope = cellfun(@(ts_date, slope, ts_name) timetable(ts_date,slope),celldates, cellr2_oos_slope,predictors','UniformOutput',0);
        cell_ts_r2_nonneg = cellfun(@(ts_date, nonneg, ts_name) timetable(ts_date,nonneg),celldates, cellr2_oos_nonneg,predictors','UniformOutput',0);
        cell_ts_r2_both = cellfun(@(ts_date, both, ts_name) timetable(ts_date,both),celldates, cellr2_oos_both,predictors','UniformOutput',0);
        r2_ts = cellfun(@(date,r2,slope,nonneg,both,name) timetable(date,r2,slope,nonneg,both),celldates, cell_ts_r2, cell_ts_r2_slope , cell_ts_r2_nonneg , cell_ts_r2_both, predictors','UniformOutput', 0); 
       
        for i = 1:length(r2_ts)
            r2_ts{i,1} = timetable(r2_ts{i,1}.date, r2_ts{i,1}.r2(:,1).r2, r2_ts{i,1}.slope(:,1).slope, r2_ts{i,1}.nonneg(:,1).nonneg, r2_ts{i,1}.both(:,1).both, 'VariableNames', {'r2', 'slope', 'nonneg', 'both'});
        end
        
        cellpredictedx = struct2cell(predicted_x);
        cellpredictedx_oos = cellfun(@(x) x.oos,cellpredictedx,'UniformOutput',0);
        cellpredictedx_oos_slope = cellfun(@(x) x.oos_slope,cellpredictedx,'UniformOutput',0);
        cellpredictedx_oos_nonneg = cellfun(@(x) x.oos_nonneg,cellpredictedx,'UniformOutput',0);
        cellpredictedx_oos_both = cellfun(@(x) x.oos_both,cellpredictedx,'UniformOutput',0);
        cell_ts_pred = cellfun(@(date,r2,slope,nonneg,both,name) timetable(date,r2,slope,nonneg,both),celldates, cellpredictedx_oos,cellpredictedx_oos_slope,cellpredictedx_oos_nonneg,cellpredictedx_oos_both,predictors','UniformOutput', 0); 
        
        predictor_ts = synchronize(cell_ts_pred{:});
        r2_ts = synchronize(r2_ts{:});
       
        [~,keep]=max(table2array(r2_ts),[],2);
        predictor_mat = table2array(predictor_ts);
        
        model = predictor_mat(sub2ind(size(predictor_mat),(1:size(predictor_mat,1))',keep));
        erp_ts.(return_label).(horizon_label) =  synchronize(erp_ts.(return_label).(horizon_label), timetable(predictor_ts.date, model));
         
        % record characteristics
        time_series = [time_series ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        monthly = [monthly ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        % model: Baker Wurgler sentiment measure
        summary = 'Baker Wurgler sentiment measure';
        Xvar = timetable(data_ts.Time, data_ts.sentiment, 'VariableNames', {'sentiment'});
        XY_ts = synchronize(realized_erp_ts,Xvar); % make sure dates for dependent and independent variables agree
        X = XY_ts.sentiment; % sentiment measure
        Y = XY_ts.realized_erp; % realized erp is left hand side variable
        dates = XY_ts.Time;
        [model, dates] = predictive_regression(Y,X,dates,12*return_horizon);
        erp_ts.(return_label).(horizon_label) =  synchronize(erp_ts.(return_label).(horizon_label), timetable(dates, model.oos));
        
        % record characteristics
        time_series = [time_series ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        monthly = [monthly ['model' num2str(model_no,'%02u')]];
        model_no = model_no+1;
        
        %% Surveys
        surveys = {}; % keep track of what models belong to this category
        
        % Duke survey of CFOs
        summary = 'Duke survey of CFOs';
        if return_horizon < 1
            % for horizons less than a year, there are no
            % forecats
            erp_ts.(return_label).(horizon_label) =  synchronize(erp_ts.(return_label).(horizon_label), timetable(erp_ts.(return_label).(horizon_label).Time, nan(size(erp_ts.(return_label).(horizon_label).Time,1),1)));
        else
            % for longer horizons, interpolate between the 1-year and
            % 10-year survey estimates
            one_year_erp =  data_ts.premium1yr;
            ten_years_erp = data_ts.premium10yr;
            alpha = (10-return_horizon)/9; % weight for the one-year survey is alpha, weight for the 10 year is 1-alpha
            if geometric_return==1
                % geometric interpolation
                model = (one_year_erp.^alpha).*ten_years_erp.^(1-alpha);
                model = real(model);
            else
                % arithmetic interpolation
                model = alpha*one_year_erp+(1-alpha)*ten_years_erp;
            end
            erp_ts.(return_label).(horizon_label) =   synchronize(erp_ts.(return_label).(horizon_label), timetable(data_ts.Time, model));
            
        end
        % record characteristics
        surveys = [surveys ['model' num2str(model_no,'%02u')]];
        model_description = [model_description summary];
        quarterly = [quarterly ['model' num2str(model_no,'%02u')]];
        
        
        %% Naming Variables
        number_of_models = size(model_description,2);
        headers1 = cellfun(@(x) ['model' num2str(x,'%02u')],mat2cell((1:number_of_models)',ones(number_of_models,1),1),'UniformOutput',0)';
        erp_ts.(return_label).(horizon_label).Properties.VariableNames = headers1;
    end
end

%% save
save([output_dir filesep 'erp_measures.mat'],'erp_ts','hist_mean','ddm','cross_section','time_series','surveys','monthly','quarterly','annual','model_description')
save([next_dir filesep 'erp_measures.mat'],'erp_ts','hist_mean','ddm','cross_section','time_series','surveys','monthly','quarterly','annual','model_description')
save([next_dir2 filesep 'erp_measures.mat'],'erp_ts','hist_mean','ddm','cross_section','time_series','surveys','monthly','quarterly','annual','model_description')


%% write to csv

% set up headers
headers1 = cellfun(@(x) ['model' num2str(x,'%02u')],mat2cell((1:number_of_models)',ones(number_of_models,1),1),'UniformOutput',0)';
headers2 = model_description;

% save all estimates of 1-year ahead ERP, one column with a recession
% dummy and one column with zeros
recession_ts = 200*data_ts.recession-100; % +100 means recession, -100 means not recession (useful for excel graphs)
save_ts = synchronize(erp_ts.(return_label).('one_year_ahead'),timetable(data_ts.Time,data_ts.SVENY01));

save_ts = synchronize(save_ts,timetable(data_ts.Time,recession_ts)); % select right ERP measure and add a column for recession dummies
save_ts = save_ts(save_ts.Time>=datetime(1959,12,1),:); % pick dates
models = table2array(save_ts);
models = [models zeros(size(models,1),1)];
modelnum = sum(~isnan(models(:,1:end-2)),2);
dates = datenum(save_ts.Time);

% creating output table
output_for_excel = table2cell(array2table(round([dates models modelnum],2)));

% replacing nans with blanks
for i = 1:size(output_for_excel,2)
    temp = output_for_excel(:,i);
    temp(cellfun(@isnan,temp)) = {[]};
    output_for_excel(:,i) = temp;
end

% creating second header row
output_for_excel(1,:) = cell(1,size(output_for_excel,2));
output_for_excel(1,2:end-4) = headers2;

% restricting to output with at least a few models
output_for_excel = output_for_excel(modelnum > 5,:);
output_for_excel = cell2table(output_for_excel);

% creating first header row
output_for_excel.Properties.VariableNames = ['date' headers1 'RF' 'recession' 'zeros' 'ModelNum'];

% formating dates
x = output_for_excel.date(2:end,:);
date = cell(size(output_for_excel,1),1);
date(1,:) = cellstr('');
date(2:end,:) = cellstr(datestr(cell2mat(x)));
output_for_excel.date = date;

%writing to csv
file_name = [excel_output filesep 'ERP_1yr_ahead.csv'];
writetable(output_for_excel, file_name);

warning(stateWarning); % restore warning state

clearvars -except root_dir

disp('Finished create_ERP_measures.m')
