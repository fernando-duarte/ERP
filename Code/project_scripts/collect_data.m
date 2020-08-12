%% collect_data

%%% This file reads all the original data and converts it into matlab
%%% financial time series objects.

stateWarning = warning('off','all'); % suppress warnings

%% set directories

% main directories
input_dir = [root_dir filesep 'Input'];
next_dir = [root_dir filesep 'Temp' filesep 'Create variables' filesep 'Input'];
excel_dir = [root_dir filesep 'Output'];

lastday = datetime('today');

%% Get Damodaran data

filename = [input_dir filesep 'histimpl.xls'];
histimpl = readtable(filename, 'Sheet', 7);

histimpl = histimpl(1:end-1,:);

dd_index = histimpl{:,14};
fcfe_index = histimpl{:,16};

date = datenum(histimpl{:,1}, 12 ,1);
datanames = {'damodaran_ddm','damodaran_fcfe'};
damodaran_ts = timetable(datetime(date,'ConvertFrom', 'datenum'), ...
                         dd_index,fcfe_index, 'VariableNames', datanames);
ind = damodaran_ts.Properties.RowTimes > lastday;
damodaran_ts = damodaran_ts(ind == 0,:);
%description
damodaran_ts.Properties.VariableDescriptions = ...
    {'ERP from a dividend discount model', ... 
    'ERP from a dividend discount model using free cash flow'};
n = width(damodaran_ts);
c(1:n) = {'annual'}; 
damodaran_ts.Properties.VariableUnits = c;

% save
save([next_dir filesep 'damodaran.mat'],'damodaran_ts')

clear histimpl year ddm_index fcfe_index date data data_init datanames 
clear description freq data_init header ending c

%% get board (bond yield) data

% nominal yield
filename = [input_dir filesep 'feds200628.csv'];
feds200628 = readtable(filename,'Headerlines', 9, 'TreatAsEmpty', {'NA'}); 

%creating values
yields_ts = table2timetable(feds200628);

% description of the data 
n = width(yields_ts);
c(1:n) = {'Bond Yields'};
yields_ts.Properties.VariableDescriptions = c;
c(1:n) = {'monthly'};
yields_ts.Properties.VariableUnits = c;
clear c

% tips 
filename = [input_dir filesep 'TIPS.csv'];
tips = readtable(filename); % get all the data
tips.DATE = datetime(tips.DATE);

tips_ts = table2timetable(tips);
ind = tips_ts.Properties.RowTimes > lastday;
tips_ts = tips_ts(ind == 0,:);
% description of the data
tips_ts.Properties.VariableDescriptions = ...
    {'yield on a 10 year TIPS';'yield on a 20 year TIPS'; ...
    'yield on a 30 year TIPS';'yield on a 5 year TIPS'; ... 
    'yield on a 7 year TIPS'};
n = length(tips_ts.Properties.VariableDescriptions);
c(1:n) = {'monthly'}; 
tips_ts.Properties.VariableUnits = c;

% save
save([next_dir filesep 'board.mat'],'yields_ts','tips_ts')

clear tips datevalue year month day dates data date ...
    feds200628 header headers header_location i c

%% get cfo data

% loading in data
filename = [input_dir filesep 'cfo_data_final.mat'];
data = load(filename);

data.data_out.premium1yr = data.data_out.ExpRet1YR - data.data_out.WGS1YR; 
data.data_out.premium10yr = data.data_out.ExpRet10YR - data.data_out.WGS10YR;
data.data_out = data.data_out(:,5:6);

cfo_ts = data.data_out;
ind = cfo_ts.Properties.RowTimes > lastday;
cfo_ts = cfo_ts(ind == 0,:);
% description of the data
cfo_ts.Properties.VariableDescriptions = ...
    {'one-year ahead ERP from CFO survey'; ...
    'ten-year ahead ERP from CFO survey'};
n = length(cfo_ts.Properties.VariableDescriptions); 
c(1:n) = {'quarterly'};
cfo_ts.Properties.VariableUnits = c;

% save
save([next_dir filesep 'cfo.mat'],'cfo_ts')

clear i month day monthnum month2num year date erp_index headers CFO 
clear filename datevalue DUKE_ERP data c

%% get fama-french data

filename = [input_dir filesep 'F-F_Research_Data_Factors.txt'];
fid = fopen(filename); % open files

% reads file until it finds headers
current_line = fgetl(fid);
while ~isempty(current_line)
    current_line = fgetl(fid);
end

% reads file
headers = textscan(fgetl(fid),'%s','Delimiter',' ', ...
    'MultipleDelimsAsOne',1);
headers = regexprep(headers{1},'[^a-zA-Z]','');
ncols = length(headers{1});
[data, ~] = textscan(fid,repmat('%f ',1,ncols),'Delimiter',...
    ' ','EmptyValue',NaN,'MultipleDelimsAsOne',1);
fclose(fid);

% create matlab dates
date = datenum(num2str(data{1}),'yyyymm');
data = cell2mat(data(2:end));
ff_ts = array2table([date data]);
ff_ts.Properties.VariableNames = {'date', 'MktRF', 'SMB', 'HML', 'RF'};
ff_ts.date = datetime(ff_ts.date, 'ConvertFrom', 'datenum');
ff_ts = table2timetable(ff_ts);
ind = ff_ts.Properties.RowTimes > lastday;
ff_ts = ff_ts(ind == 0,:);
% description of the data
ff_ts.Properties.VariableDescriptions = ...
    {'Realized excess returns for the market';'Size factor'; ...
    'Book-to-market factor';'risk free rate'};
n = length(ff_ts.Properties.VariableDescriptions); 
c(1:n) = {'monthly'};
ff_ts.Properties.VariableUnits = c;

clear c

% momentum factor
filename = [input_dir filesep 'F-F_Momentum_Factor.TXT'];
fid = fopen(filename); % open files

% reads file until it finds headers
current_line = fgetl(fid);
counter = 0;
while counter<2
    current_line = fgetl(fid);
    if isempty(current_line)
        counter = counter+1;
    end
end

% reads file
headers = textscan(fgetl(fid),'%s','Delimiter',' ', ...
    'MultipleDelimsAsOne',1);
ncols = length(headers{1});
[data, ~] = textscan(fid,repmat('%f ',1,ncols+1), ...
    'Delimiter',' ','EmptyValue',NaN,'MultipleDelimsAsOne',1);
fclose(fid);

% create matlab dates
date = datenum(num2str(data{1}),'yyyymm');
data = cell2mat(data(2));
mom_ts = array2table([date data]);
mom_ts.Properties.VariableNames = {'date', 'mom'};
mom_ts.date = datetime(mom_ts.date, 'ConvertFrom', 'datenum');
mom_ts = table2timetable(mom_ts); 
ind = mom_ts.Properties.RowTimes > lastday;
mom_ts = mom_ts(ind == 0,:);
% description of the data
mom_ts.Properties.VariableDescriptions = {'Momentum factor'};
mom_ts.Properties.VariableUnits = {'monthly'};

% size and book to market sorted portfolios
filename = [input_dir filesep '25_Portfolios_5x5.txt'];
fid = fopen(filename); % open files

% reads file until it finds headers
current_line = fgetl(fid);
while ~contains(current_line, 'Small')
    current_line = fgetl(fid);
end

% reads file
headers = textscan(fgetl(fid),'%s','Delimiter',' ', ...
    'MultipleDelimsAsOne',1);
ncols = length(headers{1});
[data, ~] = textscan(fid,repmat('%f ',1,ncols+1),'Delimiter',...
    ' ','EmptyValue',NaN,'MultipleDelimsAsOne',1);
fclose(fid);

% create matlab dates
date = datenum(num2str(data{1}),'yyyymm');
data = horzcat(data{2:26});
ff_portfolios_ts = array2table([date, data]);
headers = cell(26,1);
headers{1} = 'date';
headers(2:26,1) = cellfun(@(x) ['ff_p' num2str(x)],...
    num2cell((1:25)'),'UniformOutput',0);
ff_portfolios_ts.Properties.VariableNames = headers;
ff_portfolios_ts.date = datetime(ff_portfolios_ts.date, ...
    'ConvertFrom', 'datenum');
ff_portfolios_ts = table2timetable(ff_portfolios_ts); 


% description of the data
headers = cell(25,1);
headers(:,1) = cellfun(@(x) ['Portfolio ', num2str(x), ...
    ' of size and bm sorted portfolios'],...
    num2cell((1:25)'),'UniformOutput',0);
ff_portfolios_ts.Properties.VariableDescriptions = headers;
n = length(ff_portfolios_ts.Properties.VariableDescriptions); 
c(1:n) = {'monthly'};
ff_portfolios_ts.Properties.VariableUnits = c;
ind = ff_portfolios_ts.Properties.RowTimes > lastday;
ff_portfolios_ts = ff_portfolios_ts(ind == 0,:);
clear c

% save
save([next_dir filesep 'ff.mat'],'ff_ts','mom_ts','ff_portfolios_ts')

% momentum to market sorted portfolios
filename = [input_dir filesep '10_Portfolios_Prior_12_2.txt'];
fid = fopen(filename); % open files

% reads file until it finds headers
current_line = fgetl(fid);
while ~contains(current_line, 'Average Value') 
    current_line = fgetl(fid);
end

% reads file
headers = textscan(fgetl(fid),'%s','Delimiter',' ', ...
    'MultipleDelimsAsOne',1);
ncols = length(headers{1});
[data, ~] = textscan(fid,repmat('%f ',1,ncols+1),'Delimiter', ...
    ' ','EmptyValue',NaN,'MultipleDelimsAsOne',1);
fclose(fid);

% reads replication
filename = [input_dir filesep '10_Portfolios_Prior_12_2_Value_weighted.txt'];
file_rep = dlmread(filename,'\t',1,0);
ff_dates = data{1,1}; max_datadate = ff_dates(end);
if file_rep(end,1) > max_datadate
    for i = 1 : 11
       data{1,i} = [data{1,i}; file_rep(~ismember(file_rep(:,1), ...
           ff_dates),i)];
    end
end

% create matlab dates
date = datenum(num2str(data{1}),'yyyymm');
data = horzcat(data{2:11});
mom_portfolios_ts = array2table([date data]);
headers = cell(11,1);
headers{1} = 'date';
headers(2:11,1) = cellfun(@(x) ['mom_p' num2str(x)],...
    num2cell((1:10)'),'UniformOutput',0);
mom_portfolios_ts.Properties.VariableNames = headers;
mom_portfolios_ts.date = datetime(mom_portfolios_ts.date, ...
    'ConvertFrom', 'datenum');
mom_portfolios_ts = table2timetable(mom_portfolios_ts);
ind = mom_portfolios_ts.Properties.RowTimes > lastday;
mom_portfolios_ts = mom_portfolios_ts(ind == 0,:);
% description of the data
headers = cell(10,1);
headers(1:10,1) = cellfun(@(x) ['Portfolio ', num2str(x), ...
    ' of momentum portfolios'],...
    num2cell((1:10)'),'UniformOutput',0);
mom_portfolios_ts.Properties.VariableDescriptions = headers;
n = width(mom_portfolios_ts);
c(1:n) = {'monthly'};  
mom_portfolios_ts.Properties.VariableUnits = c;

% save
save([next_dir filesep 'ff.mat'],'ff_ts','mom_ts','ff_portfolios_ts',...
    'mom_portfolios_ts')

clear counter data date dates_duke ff_dates headers file_rep ncols 
clear end_pointer i max_datadate fid c

%% get ny fed capital markets data

filename = [input_dir filesep 'FRBNY_ERP_QS_crs.mat'];
FRBNY_ERP_QS = load(filename);
nyfed_ts = array2table(FRBNY_ERP_QS.final(:,1:2));
nyfed_ts.Properties.VariableNames = {'date', 'NYFed_ERP'};
nyfed_ts.date = datetime(nyfed_ts.date,'ConvertFrom', 'datenum'); 
nyfed_ts = table2timetable(nyfed_ts);
ind = nyfed_ts.Properties.RowTimes > lastday;
nyfed_ts = nyfed_ts(ind == 0,:);
% description of the data
nyfed_ts.Properties.VariableDescriptions = ...
    {'ERP as constructed in Adrian, Crump and Moench (2013)'};
nyfed_ts.Properties.VariableUnits = {'monthly'};

% save
save([next_dir filesep 'nyfed.mat'],'nyfed_ts')

clear FRBNY_ERP_QS data dates headers

%% get shiller data

filename = [input_dir filesep 'ie_data.xls'];
shiller = readtable(filename,'Sheet',4, 'TreatAsEmpty', {'NA', '.', ''});
shiller = shiller(:, 1:13);
shiller.Properties.VariableNames = {'date', 'P', 'D', 'E', 'CPI', ...
    'date2', 'RateGS10', 'RealPrice', 'RealDividend', ...
    'RealTotalReturnPrice', 'RealEarnings', 'RealEarningsScaled', 'CAPE'};

% deleting unecessary data
shiller.date2 = [];
shiller.RealTotalReturnPrice = [];
shiller.RealEarningsScaled = [];
% ind = strcmp(shiller.date, {''});
shiller = shiller(1:end-1,:);

%creating dates
dates = shiller{:,1};

year = ceil(dates);
month = ceil((dates-floor(dates))*100);
% ind = strcmp(month(:,2),{''});
% month = str2num(month);
% month(ind == 1,:) = 10;
dates = datenum(strcat(string(year), '-', string(month), '-', '01'), ...
    'yyyy-mm-dd');

data = zeros(size(shiller, 1), size(shiller, 2)); 
for i = 2:size(shiller,2)
    data(:,i-1) = shiller{:,i};
end
data(:,size(shiller, 2)) = dates;
data = array2table(data);
data.Properties.VariableNames = {'P', 'D', 'E', 'CPI', ...
    'RateGS10', 'RealPrice', 'RealDividend', ...
    'RealEarnings', 'CAPE', 'date'};
shiller_ts= array2table([data.date data.CAPE data.CPI, data.D, ...
    data.E, data.P, data.RateGS10, data.RealDividend, ...
    data.RealEarnings, data.RealPrice]);

shiller_ts.Properties.VariableNames = {'date', 'CAPE', 'CPI', 'D', 'E', ...
    'P', 'RateGS10', 'RealDividend', 'RealEarnings','RealPrice'};
shiller_ts.date = datetime(shiller_ts.date, 'ConvertFrom', 'datenum');
shiller_ts = table2timetable(shiller_ts); 
ind = shiller_ts.Properties.RowTimes > lastday;
shiller_ts = shiller_ts(ind == 0,:);

% description of the data
shiller_ts.Properties.VariableDescriptions = {'Cyclycally Adjusted Price-Earnings ratio'; ...
    'Consumer Price Index'; 'Nominal dividends for the S&P 500'; ...
    'Nominal earnings for the S&P 500'; 'Nominal price for the S&P 500';...
    '10 year nominal treasury yield';'Real dividends for the S&P 500';...
    'Real earnings for the S&P 500';'Real price for the S&P 500';};
c(1:width(shiller_ts)) = {'monthly'};  
shiller_ts.Properties.VariableUnits = c;

% save
save([next_dir filesep 'shiller.mat'],'shiller_ts')

clear i data data1 data2 dates headers month year days shiller
clear real_earnings real_price sp500 cpi earnings c

%% get thomson reuters data

filename = [input_dir filesep 'EPS_estimates.csv'];
fid = fopen(filename); % open files
headers = textscan(fgetl(fid),'%s','Delimiter',',');                        % get headers
headers = regexprep(headers{1},'[^a-zA-Z]','');                             % make headers suitable for variable names

% read data
ncols = length(headers);
[data, ~] = textscan(fid,repmat('%s ',1,ncols),'Delimiter',',', ...
    'EmptyValue',NaN);
fclose(fid);

% convert strings to numbers
data{strcmpi(headers,'MEANEST')} = str2double(data{strcmpi(headers,'MEANEST')});
data{strcmpi(headers,'ACTUAL')} = str2double(data{strcmpi(headers,'ACTUAL')});

% separate the time series for 1, 2, 3, 4 , 5 fiscal year forecasts
yr1_flag = strcmpi(data{strcmpi(headers,'FPI')},'1');
yr2_flag = strcmpi(data{strcmpi(headers,'FPI')},'2');
yr3_flag = strcmpi(data{strcmpi(headers,'FPI')},'3');
yr4_flag = strcmpi(data{strcmpi(headers,'FPI')},'4');
yr5_flag = strcmpi(data{strcmpi(headers,'FPI')},'5');

data1 = cellfun(@(x) x(yr1_flag),data,'UniformOutput',0);
data2 = cellfun(@(x) x(yr2_flag),data,'UniformOutput',0);
data3 = cellfun(@(x) x(yr3_flag),data,'UniformOutput',0);
data4 = cellfun(@(x) x(yr4_flag),data,'UniformOutput',0);
data5 = cellfun(@(x) x(yr5_flag),data,'UniformOutput',0);

% create time series for forecasts

date = datenum(data1{strcmpi(headers,'STATPERS')},'ddmmmyyyy');
datas = data1{strcmpi(headers,'MEANEST')};
output1 = array2table([date datas]); 
output1.Properties.VariableNames = {'date', 'EPS_y1'};
output1.date = datetime(output1.date, 'ConvertFrom', 'datenum');
output1 = table2timetable(output1);

date = datenum(data2{strcmpi(headers,'STATPERS')},'ddmmmyyyy');
datas = data2{strcmpi(headers,'MEANEST')};
output2 = array2table([date datas]); 
output2.Properties.VariableNames = {'date', 'EPS_y2'};
output2.date = datetime(output2.date, 'ConvertFrom', 'datenum');
output2 = table2timetable(output2);

date = datenum(data3{strcmpi(headers,'STATPERS')},'ddmmmyyyy');
datas = data3{strcmpi(headers,'MEANEST')};
output3 = array2table([date datas]); 
output3.Properties.VariableNames = {'date', 'EPS_y3'};
output3.date = datetime(output3.date, 'ConvertFrom', 'datenum');
output3 = table2timetable(output3);

date = datenum(data4{strcmpi(headers,'STATPERS')},'ddmmmyyyy');
datas = data4{strcmpi(headers,'MEANEST')};
output4 = array2table([date datas]); 
output4.Properties.VariableNames = {'date', 'EPS_y4'};
output4.date = datetime(output4.date, 'ConvertFrom', 'datenum');
output4 = table2timetable(output4);

date = datenum(data5{strcmpi(headers,'STATPERS')},'ddmmmyyyy');
datas = data5{strcmpi(headers,'MEANEST')};
output5 = array2table([date datas]); 
output5.Properties.VariableNames = {'date', 'EPS_y5'};
output5.date = datetime(output5.date, 'ConvertFrom', 'datenum');
output5 = table2timetable(output5);
    
% create time series for actual earnings per share
[unique_fpdates,ind1,~]= unique(data{strcmpi(headers,'FPEDATS')});
actual_eps = array2table([datenum(unique_fpdates,'ddmmmyyyy') data{strcmpi(headers,'ACTUAL')}(ind1)]);
actual_eps.Properties.VariableNames = {'date', 'EPS'};
actual_eps.date = datetime(actual_eps.date, 'ConvertFrom', 'datenum');
actual_eps = table2timetable(actual_eps); 

% merge all EPS time series and save
thomson_ts = synchronize(actual_eps, output1,output2,output3,output4,output5);
ind = thomson_ts.Properties.RowTimes > lastday;
thomson_ts = thomson_ts(ind == 0,:);
% description of the data
thomson_ts.Properties.VariableDescriptions = ...
    {'Realized Earnings per Share for the S&P 500'; ...
    'Mean analyst forecast for Earnings per share for the S&P500 for the current year'; ...
    'Mean analyst forecast for Earnings per share for the S&P500 for next year'; ...
    'Mean analyst forecast for Earnings per share for the S&P500 for 2 years ahead'; ...
    'Mean analyst forecast for Earnings per share for the S&P500 for 3 years ahead'; ...
    'Mean analyst forecast for Earnings per share for the S&P500 for 4 years ahead'};
thomson_ts.Properties.VariableUnits = {'annual';'monthly';'monthly';'monthly';'monthly';'monthly'};

% save
save([next_dir filesep 'thomson.mat'],'thomson_ts')

clear actual_eps dates date yr1 yr2 yr3 yr4 yr5 yr1_flag yr2_flag yr3_flag
clear yr4_flag yr5_flag headers data* ind1 ind2 fid ncols unique_fpdates
clear output1 output2 output3 output4 output5 
%% get fred data

% get default spread
filename = [input_dir filesep 'Default_spread.csv'];
default_spreads = readtable(filename);
default_spreads.Properties.VariableNames = {'date', 'BAA_AAA'};
default_spreads.date = datetime(default_spreads.date);
fred_ts1 = table2timetable(default_spreads);

% get recession indicator
filename = [input_dir filesep 'recession.csv'];
recession_dummy = readtable(filename);
recession_dummy.date = datenum(recession_dummy.date);
recession_dummy.date = datetime(recession_dummy.date, 'ConvertFrom', ...
    'datenum');
fred_ts2 = table2timetable(recession_dummy); 

% get short-term nominal bond yields
filename = [input_dir filesep 'Short_term_yields.csv'];
short_term_yields = readtable(filename);
short_term_yields.DATE = datetime(short_term_yields.DATE);
fred_ts3 = table2timetable(short_term_yields); 

% get 10 year yields
filename = [input_dir filesep 'tenyearyields.csv'];
tenyearyields = readtable(filename);
tenyearyields.Properties.VariableNames = {'date', 'RLONG'};

% tenyearyields.date = datetime(tenyearyields.date);
% tenyearyields.RLONG = str2double(tenyearyields.RLONG);
fred_ts4 = table2timetable(tenyearyields); 
fred_ts4 = retime(fred_ts4, 'monthly', 'lastvalue');

fred_ts = synchronize(fred_ts1, fred_ts2, fred_ts3, fred_ts4); 
fred_ts = retime(fred_ts, 'monthly','fillwithmissing');
ind = fred_ts.Properties.RowTimes > lastday;
fred_ts = fred_ts(ind == 0,:);
ind = fred_ts.Properties.RowTimes < datetime(1900,1,1);
fred_ts = fred_ts(ind == 0,:);

% description of the data
fred_ts.Properties.VariableDescriptions = {'Baa - Aaa', ...
    'NBER recession ind','yield on a 3 month bond', ...
    'yield on a 6 month bond', 'ten year yield'};
fred_ts.Properties.VariableUnits = {'monthly', 'monthly', 'monthly', ...
    'monthly', 'monthly'}; 
% save
save([next_dir filesep 'fred.mat'],'fred_ts')

clear monthnum headers data dates year month day default_spreads 
clear recession_dummy recession_ts fred_ts1 fred_ts2 fred_ts3 fred_ts4

%% get compustat data

filename = [input_dir filesep 'sp500_book_to_market.xlsx'];
book2market = readtable(filename);
book2market = book2market(:,[1 2 5]);
book2market.Properties.VariableNames = {'date', 'book_per_share', ...
    'sp500_price'};
book2market.date = datetime(book2market.date); 
compustat_ts = table2timetable(book2market);
ind = compustat_ts.Properties.RowTimes > lastday;
compustat_ts = compustat_ts(ind == 0,:);
% description of the data
compustat_ts.Properties.VariableDescriptions = ...
    {'Book value per share';'S&P 500 closing price'};
compustat_ts.Properties.VariableUnits = {'annual';'monthly'};

% save
save([next_dir filesep 'compustat.mat'],'compustat_ts')

clear monthnum month dates year day headers datevalue ...
    book2market monthstr current_line end_pointer

%% get sentiment and equity/debt issuance data

filename = [input_dir filesep 'bw_data.csv'];
investor_sentiment = readtable(filename, 'TreatAsEmpty', {'.'});
investor_sentiment = investor_sentiment(:,[1 2 8 9]);
investor_sentiment.Properties.VariableNames = {'date', 'sentiment', ...
    'equity_issuance', 'debt_issuance'};
dates = datenum(string(table2array(investor_sentiment(:,1))),'yyyymm');
dates = datetime(dates, 'ConvertFrom', 'datenum');
investor_sentiment.date = dates;

bw_ts = table2timetable(investor_sentiment);
ind = bw_ts.Properties.RowTimes > lastday;
bw_ts = bw_ts(ind == 0,:);
% description of the data
bw_ts.Properties.VariableDescriptions = {'Sentiment measure from Investor Sentiment in the Stock Market, Journal of Economic Perspectives, 2007','Equity issuance','Debt issuance'};
bw_ts.Properties.VariableUnits = {'monthly', 'monthly', 'monthly'};

% save
save([next_dir filesep 'bw.mat'],'bw_ts')

clear investor_sentiment headers data dates data_init header_location

%% get cay data

filename = [input_dir filesep 'cay_updated.csv'];
cay = readtable(filename, 'ReadVariableNames', true); % open files
cay = cay(:, [1 5]);
cay.Properties.VariableNames = {'date', 'cay'};
cay.date = datetime(cay.date);
cay_ts = table2timetable(cay);
ind = cay_ts.Properties.RowTimes > lastday;
cay_ts = cay_ts(ind == 0,:);

% description of the data
cay_ts.Properties.VariableDescriptions = {'CAY variable from Lettau and Ludvingson "Understanding Trend and Cycle in Asset Values: Reevaluating the Wealth Effect on Consumption," American Economic Review, 94 (1), 279-299, March 2004'};
cay_ts.Properties.VariableUnits = {'quarterly'};

% save
save([next_dir filesep 'cay.mat'],'cay_ts')

clear headers ncols data end_pointer year month dates current_line year
clear data datevalue s n ind tenyear*

%% save all in a single file

save([next_dir filesep 'data.mat'],'*_ts')

%% write to csv

% variables to write
variable_names = who('-regexp','_ts')';

% setting up export file 
output_for_excel = {};
for i = 1:length(variable_names)
    variable = eval(variable_names{i});
    column1 = variable.Properties.VariableNames';                           % variablenames
    column3 = variable.Properties.VariableUnits';                           % frequency
    column2 = variable.Properties.VariableDescriptions';                    % description of time series
    mat = table2array(variable);
    column4 = {};column5 = {};
    for j = 1:width(variable)
        first_nonnan = find(~isnan(mat(:,j)),1,'first');
        last_nonnan = find(~isnan(mat(:,j)),1,'last');
        column4 = [column4; max(datenum(1900,1,1), ... 
            datenum(variable.Properties.RowTimes(first_nonnan)))];
        column5 = [column5; ...
            datenum(variable.Properties.RowTimes(last_nonnan))];
    end
    output_for_excel = [output_for_excel; ...
        [column1 column2 column3 column4 column5]];
end

output_for_excel = cell2table(output_for_excel);

output_for_excel.Properties.VariableNames = {'Names', 'Source', ...         % headers of excel file
    'Frequency', 'FirstPeriod', 'LastPeriod'}; 
output_for_excel.FirstPeriod = datestr(output_for_excel.FirstPeriod);
output_for_excel.LastPeriod = datestr(output_for_excel.LastPeriod);

file_name = [excel_dir filesep 'Original_data.csv'];                        % output file name

writetable(output_for_excel, file_name);                                    % writing to csv

clearvars -except root_dir

disp('Finished collect_data.m')
