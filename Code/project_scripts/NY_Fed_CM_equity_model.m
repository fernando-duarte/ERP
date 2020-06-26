clear all

%% Setting Directories
root_dir = [filesep 'home' filesep 'rcecrs02' filesep 'ERP' filesep 'Data' filesep 'Input' ...
    filesep 'NY Fed Capital Markets' filesep];
cd(root_dir)

%% Downloading data from web

%%% downloading data

url = 'http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/25_Portfolios_5x5_CSV.zip';
filename = '/home/rcecrs02/ERP/Data/Input/NY Fed Capital Markets/ff25port_csv.zip';
outfile = websave(filename, url);

url = 'http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_Factors_TXT.zip';
filename = '/home/rcecrs02/ERP/Data/Input/NY Fed Capital Markets/F-F-Research_Data_Factors.zip';
outfile = websave(filename, url);

%%% unzipping data

! unzip "/home/rcecrs02/ERP/Data/Input/NY Fed Capital Markets/ff25port_csv.zip"

! unzip "/home/rcecrs02/ERP/Data/Input/NY Fed Capital Markets/F-F-Research_Data_Factors.zip"

%% Cleaning Data

%%% loading in data

ff25port = readtable('25_Portfolios_5x5.CSV','HeaderLines', 16);
Research_Factors = readtable('F-F_Research_Data_Factors.txt', 'ReadVariableNames', false);

%%% cleaning/organizing data

ff25port.Var1 = datenum(string(ff25port.Var1), 'yyyymm');
ff25port.Properties.VariableNames = {'date','sblh11','sblh12','sblh13', ...
    'sblh14', 'sblh15', 'sblh21', 'sblh22', 'sblh23', 'sblh24', 'sblh25', ...
    'sblh31', 'sblh32', 'sblh33', 'sblh34', 'sblh35', 'sblh41', 'sblh42', ...
    'sblh43', 'sblh44', 'sblh45', 'sblh51', 'sblh52', 'sblh53', 'sblh54', 'sblh55'};
Research_Factors_final = Research_Factors(3:1119,1:5);
Research_Factors_final.Properties.VariableNames =  {'date', 'MktRF','SMB','HML','RF'};
Research_Factors_final.date = datenum(string(Research_Factors_final.date),'yyyymm');
Research_Factors_final.MktRF = str2num(char(Research_Factors_final.MktRF));
Research_Factors_final.SMB = str2num(char(Research_Factors_final.SMB));
Research_Factors_final.HML = str2num(char(Research_Factors_final.HML));
Research_Factors_final.RF = str2num(char(Research_Factors_final.RF));

%%% downloading haver data

fbaa = readtable('fbaa.csv');
tb1m = readtable('tb1m.csv');
tb1y = readtable('tb1y.csv');
tb10y = readtable('tb10y.csv');
tb20y = readtable('tb20y.csv');
sdy5 = readtable('sdy5.csv');

fbaa.date = datenum(fbaa.date, 'ddmmmyyyy');
tb1m.date = datenum(tb1m.date, 'ddmmmyyyy');
tb1y.date = datenum(tb1y.date, 'ddmmmyyyy');
tb10y.date = datenum(tb10y.date, 'ddmmmyyyy');
tb20y.date = datenum(tb20y.date, 'ddmmmyyyy');
sdy5.date = datenum(sdy5.date, 'ddmmmyyyy');

% creating final dataset
haver_data = array2table([fbaa.date, fbaa.fbaa_usecon]);
haver_data.Properties.VariableNames = {'date', 'FBAA'};

%merging in tb1m
haver_data = outerjoin(haver_data, tb1m);
ind = (isnan(haver_data.date_haver_data));
haver_data.date_haver_data(ind==1) = haver_data.date_tb1m(ind==1);
haver_data = haver_data(:,1:3);
haver_data.Properties.VariableNames = {'date', 'FBAA', 'TB1M'};

%merging in tb1y
haver_data = outerjoin(haver_data, tb1y);
ind = (isnan(haver_data.date_haver_data));
haver_data.date_haver_data(ind==1) = haver_data.date_tb1y(ind==1);
haver_data = haver_data(:,1:4);
haver_data.Properties.VariableNames = {'date', 'FBAA', 'TB1M', 'TB1Y'};

%merging in tb10y
haver_data = outerjoin(haver_data, tb10y);
ind = (isnan(haver_data.date_haver_data));
haver_data.date_haver_data(ind==1) = haver_data.date_tb10y(ind==1);
haver_data = haver_data(:,1:5);
haver_data.Properties.VariableNames = {'date', 'FBAA', 'TB1M', 'TB1Y', 'TB10Y'};

%merging in tb20y
haver_data = outerjoin(haver_data, tb20y);
ind = (isnan(haver_data.date_haver_data));
haver_data.date_haver_data(ind==1) = haver_data.date_tb20y(ind==1);
haver_data = haver_data(:,1:6);
haver_data.Properties.VariableNames = {'date', 'FBAA', 'TB1M', 'TB1Y', 'TB10Y','TB20Y'};


%merging in sdy5
haver_data = outerjoin(haver_data, sdy5);
ind = (isnan(haver_data.date_haver_data));
haver_data.date_haver_data(ind==1) = haver_data.date_sdy5(ind==1);
haver_data = haver_data(:,1:7);
haver_data.Properties.VariableNames = {'date', 'FBAA', 'TB1M', 'TB1Y', 'TB10Y','TB20Y', 'SDY5'};

%%% creating final dataset

final_data = outerjoin(haver_data,ff25port);
ind = (isnan(final_data.date_haver_data));
final_data.date_haver_data(ind==1) = final_data.date_ff25port(ind==1);
final_data.date_ff25port = [];
final_data.Properties.VariableNames{'date_haver_data'} = 'date';

final_data = outerjoin(final_data,Research_Factors_final);
ind = (isnan(final_data.date_final_data));
final_data.date_final_data(ind==1) = final_data.date_Research_Factors_final(ind==1);
final_data.date_Research_Factors_final = [];
final_data.Properties.VariableNames{'date_final_data'} = 'date';

clear fbaa ff25port haver_data ind Research_Factors* sdy5 tb*

filename = 'Fed_data';
save(filename, 'final_data');

%% Processing Data to Compute Equity Premium

%%% adding path to functions
funcdir = [pwd filesep 'functions'];
addpath(funcdir)

%%% setting up dates
startdate = datenum('1/1/1955');
trainingdate = datenum('1/1/1960');
%enddate = datenum(year(today()), month(today())-1, 1);
enddate = datenum('1/1/2014');
freq = 12; %monthly

%%% seting options 

lags = 0;
EW = 0; 

%%% importing data

load(filename)

%%% finding start/enddates in data

datestart = find(final_data.date == startdate);
datetraining = find(final_data.date == trainingdate);
dateend = find(final_data.date == enddate);

%%% restricting data to those dates

training_data = final_data(datestart:datetraining,:);
full_data = final_data(datestart:dateend,:);

%%% matching datasets with equity premium

date = table2array(full_data(:,1));

FactorsRaw = struct();
FactorsRaw.in = [table2array(full_data(:,7)), table2array(full_data(:,2)) - table2array(full_data(:,6)), ...
    table2array(full_data(:,36)), table2array(full_data(:,33))];
FactorsRaw.dates = table2array(full_data(:,1));
FactorsRaw.names = {'DIV', 'DEF', 'RF', 'RM_RF'};

Factors = FactorsRaw;

Returns = struct();
Returns.in = [table2array(full_data(:,8:32))];
Returns.dates = table2array(full_data(:,1));
Returns.names = cell(1, size(Returns.in,2));
for i = 1:size(Returns.in,2)
    Returns.names{i} =  full_data.Properties.VariableNames{7+i};
end

FFFacs = struct();
FFFacs.in = [table2array(full_data(:,33:36))];
FFFacs.dates = table2array(full_data(:,1));
FFFacs.names = {'Mkt-RF', 'SMB', 'HML', 'RF'};

%%% Transforming Returns to Excess Returns; add RMkt

Returns.in = [FFFacs.in(:,1), ... 
    Returns.in - repmat(FFFacs.in(:,4),1,size(Returns.in,2))];
Returns.names = ['RMkt' Returns.names];
%%% Add Equal-Weighted Return Series to Returns 

if EW == 1
    Returns.in = [mean(Returns.in,2) Returns.in];
    Returns.names = ['EW Ret' Returns.names];
end

%%% Store Factors and Returns in Cells (for OOS forecasting)

Factors.rec = cell(size(Factors.in,1),1);
Returns.rec = cell(size(Returns.in,1),1);
for t=1:datetraining - datestart
    Returns.rec{t} = Returns.in(1:datetraining - datestart,:);
    Factors.rec{t} = Factors.in(1:datetraining - datestart,:);
end
for t = datetraining - datestart + 1: dateend - datestart + 1
    Returns.rec{t} = Returns.in(1:t,:);
    Factors.rec{t} = Factors.in(1:t,:);
end

%%% TSOLS and OLS

%%% matching functions date requirements

date1 = datestart;
date2 = dateend;
date3 = datetraining;

date_all = FactorsRaw.dates;
tforecast = [date(2:end); date(end)+ (date(end)-date(end-1))];

%Set TSOLS options
pricedcol = 4; % columns of the factors that are priced
forecastcol = [1 2 3]; % columns of the factors that are forecasting variables
returnscol = 1:size(Returns.in,2); % columns of returns that are priced

%Prepare IS and OOS Factors and Returns for TSOLS code
    PricedFacs.names = Factors.names(:,pricedcol); PricedFacs.in=Factors.in(:,pricedcol);
    PricedFacs.rec = cell(size(Factors.rec));
    Forecasters.names = Factors.names(:,forecastcol); Forecasters.in=Factors.in(:,forecastcol);
    Forecasters.rec = cell(size(Factors.rec));
    Rets.names = Returns.names(:,returnscol); Rets.in = Returns.in(:,returnscol);
    Rets.rec = cell(size(Returns.rec));
    for t = 1:size(Factors.rec)
        PricedFacs.rec{t} = Factors.rec{t}(:,pricedcol);
        Forecasters.rec{t} = Factors.rec{t}(:,forecastcol);
        Rets.rec{t} = Returns.rec{t}(:,returnscol);
    end
    
disp(['Forecasting Vars: ' Factors.names(forecastcol)]);
disp(['Priced Factors: ' Factors.names(pricedcol)]);
disp(['XS Returns: ' Returns.names(returnscol)]);

%%% Do Three Stage OLS

DoTSOLS 

%%% For dates w/o return data, generate Equity Risk Premium using last

%%% lambda estimate and current forecasters. 

lambda_all =  cell(size(FactorsRaw.in,1),1);

for t = 1:n
    lambda_all{t} = out.lambda.rec{t};
end
for t = n+1:size(lambda_all,1)
    lambda_all{t} = out.lambda.rec{end};
end

eq_premium = zeros(size(FactorsRaw.in,1),1);
for t = 1:size(lambda_all,1)
    eq_premium(t,1)=[1 FactorsRaw.in(t,forecastcol)]*lambda_all{t}';
end

% Correct date for output by adding a month to all values (because it's end
% of month data, eq premium is essentially for the first of the next month)
date_premium_vec = datevec(date_all); date_premium_vec(:,2) = date_premium_vec(:,2)+1;
date_premium = datenum(date_premium_vec); date_premium_str = datestr(date_premium_vec);

%%% figure to check equity premium

figure; plot(date_premium,eq_premium); datetick('x','mmm-yyyy'); title('Equity Risk Premium'); axis tight;

%%% writing output to excel





