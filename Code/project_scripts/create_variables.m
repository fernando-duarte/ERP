%% create_variables

%%% this file cleans up, homogenizes and merges all the data. Then, it
%%% creates variables that will be used to construct ERP measures.

stateWarning = warning('off','all');

%% set directories

% main directories

input_dir = [root_dir filesep 'Temp' filesep 'Create variables' filesep 'Input'];
output_dir = [root_dir filesep 'Temp' filesep 'Create variables' filesep 'Output'];
next_dir = [root_dir filesep 'Temp' filesep 'Create ERP measures' filesep 'Input'];
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
load([input_dir filesep 'data.mat'])

%% put all time series into a single one
export_to_excel = {}; % will be filled in with content to be exported to Excel

% damodaran
damodaran_clean = retime(damodaran_ts,'monthly', 'fillwithmissing'); % annualize returns and convert to monthly
damodaran_clean.damodaran_ddm = damodaran_clean.damodaran_ddm * 100;
damodaran_clean.damodaran_fcfe = damodaran_clean.damodaran_fcfe * 100;

% cfo
cfo_clean = retime(cfo_ts,'quarterly', 'previous'); 

% fama-french

ff_clean = retime(ff_ts, 'monthly', 'fillwithmissing');
mom_clean = retime(mom_ts, 'monthly', 'fillwithmissing');
ff_portfolios_clean = retime(ff_portfolios_ts, 'monthly', 'fillwithmissing');
mom_portfolios_clean = retime(mom_portfolios_ts, 'monthly', 'fillwithmissing');

ff_clean{:,:} = ff_clean{:,:}*12; %multiply by 12 to annualize returns
mom_clean{:,:} = mom_clean{:,:}*12; %multiply by 12 to annualize returns 
ff_portfolios_clean{:,:} = ff_portfolios_clean{:,:}*12; %multiply by 12 to annualize returns
mom_portfolios_clean{:,:} = mom_portfolios_clean{:,:}*12; %multiply by 12 to annualize returns


% board
yields_clean = retime(yields_ts, 'Monthly', 'previous');
tips_clear = retime(tips_ts, 'Monthly','previous' );

% ny fed
nyfed_clear = retime(nyfed_ts, 'Monthly','fillwithmissing' );

% shiller
shiller_clean = retime(shiller_ts,'Monthly','fillwithmissing'); % convert to monthly (take last day of month)

% thomson
thomson_monthly = retime(thomson_ts,'Monthly', 'previous'); % convert to monthly (take last day of month)

% linearly interpolate estimates to get 1, 2, 3 and 4 year ahead estimates
% (current estimates are all for end-of-year EPS)
[~,thomson_month] = datevec(thomson_monthly.date);
alpha = thomson_month/12;
eps_1yr_ahead = thomson_monthly.EPS_y1.*(1-alpha)+thomson_monthly.EPS_y2.*alpha;
eps_2yr_ahead = thomson_monthly.EPS_y2.*(1-alpha)+thomson_monthly.EPS_y3.*alpha;
eps_3yr_ahead = thomson_monthly.EPS_y3.*(1-alpha)+thomson_monthly.EPS_y4.*alpha;
eps_4yr_ahead = thomson_monthly.EPS_y4.*(1-alpha)+thomson_monthly.EPS_y5.*alpha;
thomson_clean = array2table([datenum(thomson_monthly.date) eps_1yr_ahead eps_2yr_ahead eps_3yr_ahead eps_4yr_ahead]);
thomson_clean.Properties.VariableNames = {'date', 'eps_1yr_ahead','eps_2yr_ahead','eps_3yr_ahead','eps_4yr_ahead'}; 
thomson_clean.date = datetime(thomson_clean.date, 'ConvertFrom', 'datenum');
thomson_clean = table2timetable(thomson_clean);

desc = {'Mean analyst forecast for Earnings per share for the S&P500 over the next 12 months'; ...
    'Mean analyst forecast for Earnings per share for the S&P500 over the next 24 months'; ...
    'Mean analyst forecast for Earnings per share for the S&P500 over the next 36 months'; ...
    'Mean analyst forecast for Earnings per share for the S&P500 over the next 48 months'};

temp_data = table2array(thomson_clean);
first_day = {}; last_day = {};
for j = 1:size(temp_data,2)
    first_nonnan = find(~isnan(temp_data(:,j)),1,'first');
    last_nonnan = find(~isnan(temp_data(:,j)),1,'last');
    first_day = [first_day; max(datenum(1900,1,1),datenum(thomson_clean.date(first_nonnan)))]; % first preiod
    last_day = [last_day; datenum(thomson_clean.date(last_nonnan))]; % last period
end
export_to_excel = [export_to_excel;[{'eps_1yr_ahead';'eps_2yr_ahead';'eps_3yr_ahead';'eps_4yr_ahead'},repmat({'ThomsonReuters I/B/E/S'},4,1),repmat({'monthly'},4,1),first_day,last_day,desc]];

% fred
fred_clean = retime(fred_ts,'Monthly', 'fillwithmissing'); % convert to monthly (take last day of month)

% compustat
compustat_clean = retime(compustat_ts,'Monthly', 'previous'); % convert to monthly (take last day of month)

% baker wurgler -- net issuance
bw_clean = retime(bw_ts,'Monthly','fillwithmissing'); % convert to monthly (take last day of month)

% lettau ludvingson -- cay
cay_clean = retime(cay_ts,'Monthly','fillwithmissing'); % convert to monthly (take last day of month)

% combine all
data_ts = synchronize(damodaran_clean,cfo_clean,ff_clean,mom_clean,ff_portfolios_clean,mom_portfolios_clean,yields_clean,tips_clear,nyfed_clear,shiller_clean,thomson_clean,fred_clean,compustat_clean,bw_clean,cay_clean);

%% restricting data to be post 1870

ind = data_ts.Time < datetime(1870,12,31);
data_ts = data_ts(ind == 0,:);

lagged = lag(data_ts, 1); %lagged data by one month

%% Create campbell-thompson variables

% dividend-price ratio
dp = 100*lagged.D./data_ts.P;
data_ts.dp = dp;

desc = {'Dividend price ratio'};
temp_data = data_ts.dp;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'dp'},{'Shiller'},{'monthly'},first_day,last_day,desc]];

% earnings-price ratio
ep = 100*lagged.E./data_ts.P;
data_ts.ep = ep;

desc = {'Earnings price ratio'};
temp_data = data_ts.ep;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'ep'},{'Shiller'},{'monthly'},first_day,last_day,desc]];

% inflation
lag_cpi = lagged.CPI;
inflation = 12*100*(data_ts.CPI./lagged.CPI -1);
data_ts.inflation = inflation;

desc = {'Log change in CPI'};
temp_data = data_ts.inflation;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'inflation'},{'Shiller'},{'monthly'},first_day,last_day,desc]];

% term spread
term_spread = data_ts.SVENY10-data_ts.SVENY01;
data_ts.term_spread = term_spread;

desc = {'10 year minus 1 year yield'};
temp_data = data_ts.term_spread;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'term_spread'},{'Fed Board'},{'monthly'},first_day,last_day,desc]];

% net issuance
share_equity_issuance = 100*data_ts.equity_issuance./(data_ts.equity_issuance+data_ts.debt_issuance);
data_ts.share_equity_issuance = share_equity_issuance;

desc = {'Share of equity issuance'};
temp_data = data_ts.share_equity_issuance;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'share_equity_issuance'},{'Jeffrey Wurgler'},{'monthly'},first_day,last_day,desc]];

% book to market

% assign to each month the value of book_per_share in december of the
% current year (as in Goyal and Welch)
data_ts.book_per_share = circshift(data_ts.book_per_share,-12);
desc = {'Assign to each month of the year the value of books_per_share in December of the same year (as in Goyal and Welch)'};
temp_data = data_ts.book_per_share;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'book_per_share'},{'Compustat'},{'monthly'},first_day,last_day,desc]];

% assign to each month the value of book_per_share in december of the previous year
lag12m = lag(data_ts, 12);
data_ts.lag_book_per_share = lag12m.book_per_share;

desc = {'Assign to each month of the year the value of books_per_share in December of the previous year'};
temp_data = data_ts.lag_book_per_share;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'lag_book_per_share'},{'Compustat'},{'monthly'},first_day,last_day,desc]];

% nominal book equity
data_ts.book_equity = data_ts.lag_book_per_share;

desc = {'Is the same as lag_book_per_share (created to follow Campbell Thompson names)'};
temp_data = data_ts.book_equity;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'book_equity'},{'Compustat'},{'monthly'},first_day,last_day,desc]];

% real book equity
data_ts.real_book_equity = data_ts.book_equity./data_ts.CPI;

desc = {'book_equity/CPI'};
temp_data = data_ts.real_book_equity;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max((datenum(1900,1,1)),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'real_book_equity'},{'Compustat and Shiller'},{'monthly'},first_day,last_day,desc]];

% construct book-to-market
data_ts.bm = 100*data_ts.lag_book_per_share./data_ts.sp500_price;

desc = {'lag_book_per_share divided by sp500_price'};
temp_data = data_ts.bm;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'bm'},{'Compustat'},{'monthly'},first_day,last_day,desc]];

% smooth earnings price ratio (1/CAPE)
data_ts.smooth_EP = 100./data_ts.CAPE;

desc = {'1/CAPE'};
temp_data = data_ts.smooth_EP;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'smooth_EP'},{'Shiller'},{'monthly'},first_day,last_day,desc]];

% return on equity
roe = data_ts.E./((1+data_ts.inflation/(100*12)).*(data_ts.real_book_equity));
roe = filter(ones(120,1)/120,1,roe);
data_ts.roe = roe;

desc = {'real return on equity'};
temp_data = data_ts.roe;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'roe'},{'Compustat and Shiller'},{'monthly'},first_day,last_day,desc]];

% historical average of the risk-free rate
Real_Rfree = data_ts.RateGS10-data_ts.inflation;
nan_temp = isnan(Real_Rfree);
historical_Real_Rfree = nancumsum(Real_Rfree)./(1:length(Real_Rfree))';
historical_Real_Rfree(nan_temp) = NaN;
data_ts.real_rf_mean = historical_Real_Rfree/10; 

desc = {'Historical average of the real risk free rate'};
temp_data = data_ts.real_rf_mean;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'real_rf_mean'},{'Fama-French and Shiller'},{'monthly'},first_day,last_day,desc]];

% historical average of the payout ratio
payout = data_ts.D./data_ts.E;
nan_temp = isnan(payout);
historical_payout = 100*cumsum(payout,'omitnan')./(1:length(payout))';
historical_payout(nan_temp) = NaN;
data_ts.payout_mean = historical_payout; 

desc = {'Historical average of the payout ratio D/E'};
temp_data = data_ts.payout_mean;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'payout_mean'},{'Shiller'},{'monthly'},first_day,last_day,desc]];

% historical average of real earnings growth
real_earnings = data_ts.E./data_ts.CPI;
earnings_growth = 100*12*(real_earnings(2:end)./real_earnings(1:end-1) - 1);
nan_temp = isnan(earnings_growth);
earnings_growth = cumsum(earnings_growth, 'omitnan')./(1:length(earnings_growth))';
earnings_growth(nan_temp) = NaN;
earnings_growth(end+1) = NaN;
data_ts.real_earnings_gr_mean = earnings_growth; 

desc = {'Historical average of real earnings growth'};
temp_data = data_ts.real_earnings_gr_mean;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'real_earnings_gr_mean'},{'Shiller'},{'monthly'},first_day,last_day,desc]];

% growth-adjusted dividend-price ratio
data_ts.dp_growth_adj = data_ts.dp+(1-data_ts.payout_mean/100).*data_ts.roe;

desc = {'Growth-adjusted dividend-price ratio'};
temp_data = data_ts.dp_growth_adj;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'dp_growth_adj'},{'Compustat and Shiller'},{'monthly'},first_day,last_day,desc]];

% growth-adjusted earnings-price ratio
data_ts.ep_growth_adj = (data_ts.payout_mean/100).*data_ts.ep+(1-data_ts.payout_mean/100).*data_ts.roe;

desc = {'Growth-adjusted earnings-price ratio'};
temp_data = data_ts.ep_growth_adj;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day = datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'ep_growth_adj'},{'Compustat and Shiller'},{'monthly'},first_day,last_day,desc]];

% growth-adjusted smooth earnings-price ratio
data_ts.smooth_ep_growth_adj = (data_ts.payout_mean/100).*data_ts.smooth_EP+(1-data_ts.payout_mean/100).*data_ts.roe;

desc = {'Growth-adjusted smooth earnings-price ratio'};
temp_data = data_ts.smooth_ep_growth_adj;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'smooth_ep_growth_adj'},{'Compustat and Shiller'},{'monthly'},first_day,last_day,desc]];

% growth-adjusted book-to-market ratio
data_ts.bm_growth_adj = data_ts.roe.*(1+data_ts.payout_mean/100).*(data_ts.bm/100)-1;

desc = {'Growth-adjusted book-to-market ratio'};
temp_data = data_ts.bm_growth_adj;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'bm_growth_adj'},{'Compustat and Shiller'},{'monthly'},first_day,last_day,desc]];

% growth-adjusted dividend-price ratio minus the real rate
data_ts.dp_growth_adj_rf = data_ts.dp_growth_adj-data_ts.real_rf_mean;

desc = {'Growth-adjusted dividend-price ratio minus the risk-free rate'};
temp_data = data_ts.dp_growth_adj_rf;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'dp_growth_adj_rf'},{'Compustat, Shiller and Fama-French'},{'monthly'},first_day,last_day,desc]];

% growth-adjusted earnings-price ratio minus the real rate
data_ts.ep_growth_adj_rf = data_ts.ep_growth_adj-data_ts.real_rf_mean;

desc = {'Growth-adjusted earnings-price ratio minus the risk-free rate'};
temp_data = data_ts.ep_growth_adj_rf;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'ep_growth_adj_rf'},{'Compustat, Shiller and Fama-French'},{'monthly'},first_day,last_day,desc]];

% growth-adjusted smooth earnings-price ratio minus the real rate
data_ts.smooth_ep_growth_adj_rf = data_ts.smooth_ep_growth_adj-data_ts.real_rf_mean;

desc = {'Growth-adjusted smooth earnings-price ratio minus the risk-free rate'};
temp_data = data_ts.smooth_ep_growth_adj_rf;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'smooth_ep_growth_adj_rf'},{'Compustat, Shiller and Fama-French'},{'monthly'},first_day,last_day,desc]];

% growth-adjusted book-to-market ratio minus the real rate
data_ts.bm_growth_adj_rf = data_ts.bm_growth_adj-data_ts.real_rf_mean;

desc = {'Growth-adjusted book-to-market ratio minus the risk-free rate'};
temp_data = data_ts.bm_growth_adj_rf;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'bm_growth_adj_rf'},{'Compustat, Shiller and Fama-French'},{'monthly'},first_day,last_day,desc]];

% real dividend growth
data_ts.real_D_growth = 100*12*(data_ts.RealDividend./lagged.RealDividend-1);

desc = {'Growth rate of real dividends'};
temp_data = data_ts.real_D_growth;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'real_D_growth'},{'Shiller'},{'monthly'},first_day,last_day,desc]];

% real earnings growth
data_ts.real_E_growth = 100*12*(data_ts.RealEarnings./lagged.RealEarnings-1);

desc = {'Growth rate of real earnings'};
temp_data = data_ts.real_E_growth;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'real_E_growth'},{'Shiller'},{'monthly'},first_day,last_day,desc]];

% real price growth (ex-dividend real return)
data_ts.real_P_growth = 100*12*(data_ts.RealPrice./lagged.RealPrice-1);

desc = {'Growth rate of the real price of the S&P500 (i.e. ex-dividend returns)'};
temp_data = data_ts.real_P_growth;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'real_P_growth'},{'Shiller'},{'monthly'},first_day,last_day,desc]];

% dividend yield
data_ts.dividend_yield = 100*data_ts.D./lagged.P;

desc = {'Dividend divided by lagged S&P500 price'};
temp_data = data_ts.dividend_yield;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first preiod
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'dividend_yield'},{'Shiller'},{'monthly'},first_day,last_day,desc]];

% realized equity risk premium
data_ts.realized_erp = 100*12*((data_ts.P+data_ts.D/12)./lagged.P-1)-data_ts.RF;

desc = {'Realized excess market returns (cum dividend)'};
temp_data = data_ts.realized_erp;
first_nonnan = find(~isnan(temp_data),1,'first');
last_nonnan = find(~isnan(temp_data),1,'last');
first_day = max(datenum(1900,1,1),datenum(data_ts.Time(first_nonnan))); % first period
last_day =  datenum(data_ts.Time(last_nonnan)); % last period
export_to_excel = [export_to_excel;[{'realized_erp'},{'Shiller'},{'monthly'},first_day,last_day,desc]];

%% save
save([output_dir filesep 'variables.mat'],'data_ts')
save([next_dir filesep 'variables.mat'],'data_ts')

%% write to csv
file_name = [excel_output filesep 'Constructed_variables.csv'];
export_to_excel = cell2table(export_to_excel);
export_to_excel.Properties.VariableNames = {'Names', 'Source', 'Frequency', ...
    'FirstPeriod', 'LastPeriod', 'Description'}; % headers of excel file
export_to_excel.FirstPeriod = datestr(export_to_excel.FirstPeriod);
export_to_excel.LastPeriod = datestr(export_to_excel.LastPeriod);
writetable(export_to_excel, file_name);

warning(stateWarning); % restore warning state

clearvars -except root_dir

disp('Finished create_variables.m')