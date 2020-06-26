function [predicted_x, dates, r2, coefficients, residual, bint] = predictive_regression(Y,X,dates,horizon)
%PREDICTIVE_REGRESSION Runs linear regressions predicting Y(t+k) using X(t) with k>0
%
%   PREDICTED_X = PREDICTIVE_REGRESSION(Y,X,DATES,HORIZON) takes as inputs a T-by-1 vector Y,
%   a T-by-1 vector X, a T-by-1 vector DATES and a positive integer HORIZON (smaller than T). The
%   function returns a stucture PREDICTED_X with 6 fields that give different ways of
%   predicting Y(t+HORIZON) using X(t):
%       historical_average: the average of X upto the current observation
%       is: the predicted values of Y in a regression of Y(t+HORIZON) on X(t) using a
%       full-sample OLS regression
%       oos: the predicted values of Y in a regression of Y(t+HORIZON) on X(t) using all
%       observations up to time t only
%       oos_slope: same as oos but replaces the slope coefficient by zero
%       if it is estimated to be negative in any of the regressions
%       oos_nonneg: same as oos but replaces the estimate of Y(t+HORIZON) by zero
%       if it is estimated to be negative in any of the regressions
%       oos_nonneg: same as oos but replaces the the slope coefficient by zero
%       if it is estimated to be negative in any of the regressions, and
%       then replaces the estimate of Y(t+HORIZON) by zero
%       if it is estimated to be negative
%   Note: if any observation of X or Y is NaN, then the observation is
%   removed before any regressions are run
%
%   [PREDICTED_X DATES]= PREDICTIVE_REGRESSION(...) also returns the input DATES
%   but with observations for which X or Y have NaN's removed
%
%   [PREDICTED_X DATES R2 COEFFICIENTS]= PREDICTIVE_REGRESSION(...) also returns a
%   structure R2 with the same 6 fields of PREDICTED_X  that contain the
%   R-square of the regression (in this context, the R^2 is defined
%   relative to the historical average, ie R^2 = 1 -
%   sse_predict/sse_historical, where sse_predict is the sum of squared
%   errors of the predictive model and sse_historical is the sum of squared
%   errors using only the historical mean of X as a predictor
%
%   [PREDICTED_X DATES R2 COEFFICIENTS]= PREDICTIVE_REGRESSION(...) also returns a
%   structure COEFFICIENTS with the same 6 fields of PREDICTED_X that contain the
%   regression coefficients of the predictive regressions
%
%   [PREDICTED_X DATES R2 COEFFICIENTS RESIDUALS]= PREDICTIVE_REGRESSION(...) also returns a
%   structure RESIDUALS with the same 6 fields of PREDICTED_X that contain the
%   residuals of the predictive regressions
%
%   [PREDICTED_X DATES R2 COEFFICIENTS RESIDUALS BINT]= PREDICTIVE_REGRESSION(...) also returns a
%   structure BINT with confidence intervals for COEFFICIENTS

%% eliminate nan's
nanXY = any(isnan(X),2) | isnan(Y);
X = X(~nanXY);Y = Y(~nanXY);

% if not enough observations or horizon is smaller than a single period, then exit
if length(Y)-horizon<2  || horizon<1 || isempty(dates(~nanXY))
    dates = dates(end);
    predicted_x.historical_average = NaN;
    predicted_x.is =  NaN;
    predicted_x.oos =  NaN;
    predicted_x.oos_slope =  NaN;
    predicted_x.oos_nonneg =  NaN;
    predicted_x.oos_both = NaN;
    
    residual.historical_average = NaN;
    residual.is =  NaN;
    residual.oos =  NaN;
    residual.oos_slope =  NaN;
    residual.oos_nonneg =  NaN;
    residual.oos_both = NaN;
    
    r2.historical_average =  NaN;
    r2.is =  NaN;
    r2.oos =  NaN;
    r2.oos_slope =  NaN;
    r2.oos_nonneg =  NaN;
    r2.oos_both =  NaN;
    coefficients = NaN;
    residual = NaN;
    
    bint = NaN;
    return
    
else % run the regressions

dates = dates(~nanXY);
last_date = dates(end);

%% initialize variables
coefficients = cell(length(X),1);
bint = cell(length(X),1);
predicted_erp = nan(length(X),1);
predicted_erp_slope = nan(length(X),1);
predicted_erp_nonneg = nan(length(X),1);
predicted_erp_both = nan(length(X),1);
sqr_err_historical_mean = nan(length(X),1);
mean_sqr_err_historical_mean = nan(length(X),1);
predicted_historical_average = nan(length(X),1);
residual.historical_average  = nan(length(X),1);
sqr_err_full_sample = nan(length(X),1);
mean_sqr_err_full_sample = nan(length(X),1);
residual.is = nan(length(X),1);
sqr_err_oos = nan(length(X),1);
mean_sqr_err_oos = nan(length(X),1);
residual.oos = nan(length(X),1);
sqr_err_oos_slope = nan(length(X),1);
mean_sqr_err_oos_slope = nan(length(X),1);
residual.oos_slope = nan(length(X),1);
sqr_err_oos_nonneg = nan(length(X),1);
mean_sqr_err_oos_nonneg = nan(length(X),1);
residual.oos_nonneg = nan(length(X),1);
sqr_err_oos_both = nan(length(X),1);
mean_sqr_err_oos_both = nan(length(X),1);
residual.oos_both  = nan(length(X),1);

%% run predictive regressions


    
    for t=2:length(Y)-horizon % start at t=2 to have at least 2 observations
        [coefficients{t+horizon}, bint{t+horizon}]=regress(Y(1+horizon:t+horizon),[ones(size(Y(1+horizon:t+horizon))) X(1:t,:)]);
        
        % unrestricted regression
        predicted_erp(t+horizon+1) = [1 X(t)]*coefficients{t+horizon};
        
        % restricted regression
        predicted_erp_slope(t+horizon+1) = coefficients{t+horizon}(1)+X(t)*max(0,coefficients{t+horizon}(2));  % restrict to positive slope
        predicted_erp_nonneg(t+horizon+1) = max(0,[1 X(t)]*coefficients{t+horizon}); % positive estimate
        predicted_erp_both(t+horizon+1) = max(0,predicted_erp_slope(t+horizon+1)); % first correct sign of slope, then positive estimate
        
    end
end
%% get R^2

% historical average
full_historical_average = cumsum(Y)./(1:length(Y))';
predicted_historical_average(horizon+1:end)=full_historical_average(1:end-horizon);
residual.historical_average(horizon+1:end) = Y(horizon+1:end)-full_historical_average(1:end-horizon);
sqr_err_historical_mean(horizon+1:end) = ((residual.historical_average(horizon+1:end))/100).^2;
%mean_sqr_err_historical_mean(horizon+1:end) = cumsum(sqr_err_historical_mean(horizon+1:end))./(1:length(sqr_err_historical_mean(horizon+1:end)))';
mean_sqr_err_historical_mean(3+horizon:end) = cumsum(sqr_err_historical_mean(3+horizon:end))./(1:length(sqr_err_historical_mean(3+horizon:end)))';

% full sample regression
predicted_full_ols = [ones(size(X)) X]*coefficients{end};
residual.is(horizon+1:end) = Y(horizon+1:end)-predicted_full_ols(1:end-horizon);
sqr_err_full_sample(horizon+1:end) = ((residual.is(horizon+1:end))/100).^2;
%mean_sqr_err_full_sample(horizon+1:end) = cumsum(sqr_err_full_sample(horizon+1:end))./(1:length(sqr_err_full_sample(horizon+1:end)))';
mean_sqr_err_full_sample(3+horizon:end) = cumsum(sqr_err_full_sample(3+horizon:end))./(1:length(sqr_err_full_sample(3+horizon:end)))';

% out-of-sample regression, unconstrained
residual.oos(3+horizon:end)  = Y(3+horizon:end)-predicted_erp(3+horizon:end-1);
sqr_err_oos(3+horizon:end)  = ((residual.oos(3+horizon:end))/100).^2;
mean_sqr_err_oos(3+horizon:end) = cumsum(sqr_err_oos(3+horizon:end))./(1:length(sqr_err_oos(3+horizon:end)))';

% out-of-sample regression, positive slope
residual.oos_slope(3+horizon:end) = Y(3+horizon:end)-predicted_erp_slope(3+horizon:end-1);
sqr_err_oos_slope(3+horizon:end)  = ((residual.oos_slope(3+horizon:end))/100).^2;
mean_sqr_err_oos_slope(3+horizon:end) = cumsum(sqr_err_oos_slope(3+horizon:end))./(1:length(sqr_err_oos_slope(3+horizon:end)))';

% out-of-sample regression, positive erp
residual.oos_nonneg(3+horizon:end)  = Y(3+horizon:end)-predicted_erp_nonneg(3+horizon:end-1);
sqr_err_oos_nonneg(3+horizon:end)  = ((residual.oos_nonneg(3+horizon:end) )/100).^2;
mean_sqr_err_oos_nonneg(3+horizon:end) = cumsum(sqr_err_oos_nonneg(3+horizon:end))./(1:length(sqr_err_oos_nonneg(3+horizon:end)))';

% out-of-sample regression, positive slope then positive erp
residual.oos_both(3+horizon:end) = Y(3+horizon:end)-predicted_erp_both(3+horizon:end-1);
sqr_err_oos_both(3+horizon:end)  = ((residual.oos_both(3+horizon:end))/100).^2;
mean_sqr_err_oos_both(3+horizon:end) = cumsum(sqr_err_oos_both(3+horizon:end))./(1:length(sqr_err_oos_both(3+horizon:end)))';

% find R-squared relative to historical mean
is_R2 = 1-mean_sqr_err_full_sample./mean_sqr_err_historical_mean;
oos_R2 = 1-mean_sqr_err_oos./mean_sqr_err_historical_mean;
oos_R2_slope = 1-mean_sqr_err_oos_slope./mean_sqr_err_historical_mean;
oos_R2_nonneg = 1-mean_sqr_err_oos_nonneg./mean_sqr_err_historical_mean;
oos_R2_both = 1-mean_sqr_err_oos_both./mean_sqr_err_historical_mean;

%% put in desired form for output

nan_all = isnan(is_R2) | isnan(oos_R2) | isnan(oos_R2_slope)| isnan(oos_R2_nonneg)| isnan(oos_R2_both);
dates = dates(~nan_all);

predicted_x.historical_average = predicted_historical_average(~nan_all);
predicted_x.is = predicted_full_ols(~nan_all);
predicted_x.oos = predicted_erp(~nan_all);
predicted_x.oos_slope = predicted_erp_slope(~nan_all);
predicted_x.oos_nonneg = predicted_erp_nonneg(~nan_all);
predicted_x.oos_both = predicted_erp_both(~nan_all);

residual.historical_average =residual.historical_average(~nan_all);
residual.is = residual.is(~nan_all);
residual.oos = residual.oos(~nan_all);
residual.oos_slope = residual.oos_slope(~nan_all);
residual.oos_nonneg = residual.oos_nonneg(~nan_all);
residual.oos_both = residual.oos_both(~nan_all);

r2.historical_average = zeros(size(dates));
r2.is = is_R2(~nan_all);
r2.oos = oos_R2(~nan_all);
r2.oos_slope = oos_R2_slope(~nan_all);
r2.oos_nonneg = oos_R2_nonneg(~nan_all);
r2.oos_both = oos_R2_both(~nan_all);

if isempty(dates)
    dates = last_date;
    predicted_x.historical_average = NaN;
    predicted_x.is =  NaN;
    predicted_x.oos =  NaN;
    predicted_x.oos_slope =  NaN;
    predicted_x.oos_nonneg =  NaN;
    predicted_x.oos_both = NaN;
    
    residual.historical_average = NaN;
    residual.is =  NaN;
    residual.oos =  NaN;
    residual.oos_slope =  NaN;
    residual.oos_nonneg =  NaN;
    residual.oos_both = NaN;
    
    r2.historical_average =  NaN;
    r2.is =  NaN;
    r2.oos =  NaN;
    r2.oos_slope =  NaN;
    r2.oos_nonneg =  NaN;
    r2.oos_both =  NaN;
    
    bint = NaN;
    
end
