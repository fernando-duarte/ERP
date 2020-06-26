function y = change_horizon(data,times,new_frequency,return_horizon,geometric_return,names)
%CHANGE_HORIZON Computes returns at different horizons using different
%aggregation methods
%
%   Y = CHANGE_HORIZON(data,new_frequency) returns the
%   time-series object Y_FS at frequency new_frequency. The time-series
%   data must contain annualized returns expressed in percentage points.
%   The number new_frequency is in units of years, so 1/12 means monthly,
%   1/4 means quarterly, 1 means annual. new_frequency must be larger than
%   or equal to the frequency of data. data is assumed to have no gaps
%   (i.e. there can be NaN's at the beginning and the end but no NaN's
%   between two non-NaN numbers). Allowed frequencies are monthly,
%   quarterly, semi-annual, annual and 2-5 years.
%
%   Y = CHANGE_HORIZON(data,new_frequency,return_horizon) specifies
%   the horizon at which returns are calculated. return_horizon must be
%   larger than or equal to the frequency of data. For example, if data has
%   frequency 1/12 (monthly), new_frequency is 1/4 (quarterly) and
%   return_horizon is 1 (annual), then Y_TS is a time series that has one
%   observation per quarter, and each observation is the returns over the
%   previous year (i.e. over 12 months in data). If left unspecified,
%   return_horizon is set to new_frequency.
%
%   Y = CHANGE_HORIZON(data,new_frequency,return_horizon,geometric_return)
%   computes returns using geometric averages if geometric_average is 0, and
%   geometric returns if geometric_average is 1. Default is 0.

%% check inputs and set defaults

if ~exist('return_horizon','var')
    return_horizon = new_frequency;
elseif isempty(return_horizon)
    return_horizon = new_frequency;
end

if ~exist('geometric_return','var')
    geometric_return = 0;
elseif isempty(geometric_return)
    geometric_return = 0;
end

if return_horizon<new_frequency
    error('return_horizon must be greater than or equal to new_frequency')
end

%% prepare data
dt = datenum(times(2))-datenum(times(1));

if 28<=dt && dt<=31
    old_frequency = 1/12; % monthly
elseif 31<dt && dt<=93
    old_frequency = 1/4; % quarterly
elseif 93<dt && dt<=186
    old_frequency = 1/2; % semi-annual
elseif 186<dt && dt<=366
    old_frequency = 1; % annual
elseif 366<dt && dt<=2*366
    old_frequency = 2;
elseif 2*366<dt && dt<=3*366
    old_frequency = 3;
elseif 3*366<dt && dt<=4*366
    old_frequency = 4;
elseif 4*366<dt && dt<=5*366
    old_frequency = 5;
else
    error('Frequency of time series is too small')
end

if new_frequency<old_frequency
    error('new_frequency must be greater than or equal to old_frequency')
end

%% compute returns over desired horizon

num_series = size(data,2);

for i=1:num_series
    returns = data(:,i);
    
    if ~isnan(returns(1))
        first_nan=1;
    else
        first_nan = find(~isnan(returns),1,'first'); % find first nan
    end
    if ~isnan(returns(end))
        last_nan=size(returns,1);
    else
        last_nan = find(~isnan(returns),1,'last'); % find last nan
    end
    keep = first_nan:last_nan; % index of non-nan observations
    
    return_dates = datenum(times(keep));
    returns = returns(keep);
    nperiods = return_horizon/old_frequency; % number of periods in original series contained in one period of the new series
    
    if geometric_return==1
        
        val = returns/(100/old_frequency);
        price = zeros(length(val)+1,1);
        price(1) = 1;
        for j = 2:length(val)+1
            price(j) = price(j-1)*(1+val(j-1));
        end  % convert returns to price
        new_returns = (price(1+nperiods:end)./price(1:end-nperiods)-1)*100/return_horizon; % take geometric returns over nperiods and annualize
        new_returns_dates = return_dates(nperiods:end); 
    elseif geometric_return==0
        new_returns = filter(ones(nperiods,1)/nperiods,1,returns); % take trailing moving average over nperiods
        new_returns = new_returns(nperiods:end);
        new_returns_dates = return_dates(nperiods:end);  
    end
    
    if i==1
        y = timetable(datetime(new_returns_dates, 'ConvertFrom', 'datenum'),new_returns);
    else
        z = timetable(datetime(new_returns_dates, 'ConvertFrom', 'datenum'),new_returns);
        y = synchronize(y,z);
    end
end

%% Getting variable names

y.Properties.VariableNames = names;

%% convert to new frequency

keep = ~all(isnan(table2array(y)),2);
y = y(keep,:);
nperiods = new_frequency/old_frequency; % number of periods in original series contained in one period of the new series

if new_frequency==1/12
    first_month = 1:12;
elseif new_frequency==1/4
    first_month = 3;
elseif new_frequency==1/2
    first_month = 6;
else
    first_month = 12;
end
[~,month] = datevec(y.Time);

y = y(find(ismember(month,first_month),1,'first'):nperiods:end,:);

    
    
    
    
    
    
    
    
    

