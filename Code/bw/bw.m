%% Set Directory 

root_dir = [filesep 'home' filesep 'rcecrs02' filesep 'ERPv3'];
stata_dir = [filesep 'home' filesep 'rcecrs02' filesep 'ERPv3' filesep 'Data'];
cd(root_dir)

%% Add the Matlab APIs to path
wrds_api_path = 'Matlab subfunctions/WRDS API';
fred_api_path = 'Matlab subfunctions/FredFetch-master';
matlab_subfunctions = [root_dir filesep 'Matlab subfunctions'];
addpath(matlab_subfunctions)
addpath(genpath(wrds_api_path)); addpath(genpath(fred_api_path));

%% Download WRDS Data

% Connect to the WRDS server
wrds_username = 'smitch15';
wrds_password = 'D@nkey10';
a.SSH2conn = ssh2_config('wrds.wharton.upenn.edu',wrds_username,wrds_password,22);

% Transfer data to appropriate folders. 
%scp_get(a.SSH2conn,'firm_data.csv', 'Data/Input/Baker Wurgler/','/home/frb-ny/smitch15/wrds/DATA/value/');
%scp_get(a.SSH2conn,'firm_data2.csv', 'Data/Input/Baker Wurgler/','/home/frb-ny/smitch15/wrds/DATA/value/');
%scp_get(a.SSH2conn,'m2b.csv', 'Data/Input/Baker Wurgler/','/home/frb-ny/smitch15/wrds/DATA/value/');
%scp_get(a.SSH2conn,'m2b2.csv', 'Data/Input/Baker Wurgler/','/home/frb-ny/smitch15/wrds/DATA/value/');
%scp_get(a.SSH2conn,'firm_data3.csv', 'Data/Input/Baker Wurgler/','/home/frb-ny/smitch15/wrds/DATA/value/');
%scp_get(a.SSH2conn,'firm_data4.csv', 'Data/Input/Baker Wurgler/','/home/frb-ny/smitch15/wrds/DATA/value/');
scp_get(a.SSH2conn,'firm_data5.csv', 'Data/Input/Baker Wurgler/','/home/frb-ny/smitch15/wrds/DATA/value/');
ssh2_close(a.SSH2conn);

%% BW Investor Sentiment with updating data

%% Data Downloads

% first day ipo returns weighted by ipo volume

url = 'https://site.warrington.ufl.edu/ritter/files/2019/01/IPOALL_2018.xlsx';
filename = 'Data/Input/Baker Wurgler/ritter_ipo_data.xlsx';
outfile = websave(filename,url);

ritter_ipo = readtable(filename, 'ReadVariableNames', 0);
endindex = find(isnan(ritter_ipo.Var2));
ritter_ipo = ritter_ipo(1:endindex-1,1:6);

ind2000 = (ritter_ipo.Var2 == [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 ...
    18 19 20 21 22 23]);
ind2000 = sum(ind2000,2);
ritter_ipo.year(ind2000 == 0) = 19;
ritter_ipo.year(ind2000 == 1) = 20;
ritter_ipo.Var2 = num2str(ritter_ipo.Var2,'%02d');
ritter_ipo.year = [num2str(ritter_ipo.year) num2str(ritter_ipo.Var2)];
date = datenum(str2num(ritter_ipo.year), str2num(char(ritter_ipo.Var1))+1, 1);
ripo = str2double(string(ritter_ipo.Var3));
nipo = ritter_ipo.Var4;

ipo_data = struct();
ipo_data.date = date-1;
ipo_data.ripo = ripo;
ipo_data.nipo = nipo;

%%% measures are avaliable on FRED

fred = fred.latest({'INDPRO','PCEDG','PCEND','PCES','USREC','PAYEMS','CPIAUCSL'});

date = fred.date;
indpro = fred.value(:,1);
consdur = fred.value(:,2);
consnon = fred.value(:,3);
consserv = fred.value(:,4);
recess = fred.value(:,5);
employ = fred.value(:,6);
cpi = fred.value(:,7);

fred_data = struct();
fred_data.date = datenum(year(date),month(date)+1,day(date))-1;
fred_data.date = fred_data.date;
fred_data.indpro = indpro;
fred_data.consdur = consdur;
fred_data.consnon = consnon;
fred_data.consserv = consserv;
fred_data.recess = recess;
fred_data.employ = employ;
fred_data.cpi = cpi; 

%%% equity/debt issuance

% url = 'https://fweb.rsma.frb.gov/cgi-bin/fame/fame.pl?output=CSV&bdate=%27*%27&edate=%27*%27&deci=Auto&transform=level&freq=%27no%20conversion%27&selection=US%27EQFSNBCRSS_N.M&display_series=YES';
% % filename = 'Data/Input/Baker Wurgler/equity.csv';
% % websave(filename, url);
% % 
% url = 'https://fweb.rsma.frb.gov/cgi-bin/fame/fame.pl?output=CSV&bdate=%27*%27&edate=%27*%27&deci=Auto&transform=level&freq=%27no%20conversion%27&selection=US%27DTBCFBNMD_N.M&display_series=YES';
% filename = 'Data/Input/Baker Wurgler/debt.csv';
% o = weboptions('CertificateFilename', '');
% webread(url, o)
% outfile = websave(filename, url); 

equity = readtable('Data/Input/Baker Wurgler/equity_share.xlsx');
dateval = char(string(equity.yearmo));
yearval = str2num(dateval(:,1:4));
monthval = str2num(dateval(:,5:6)); 

date = datenum(yearval, monthval+1, 1);
date = date - 1;

equity_data = struct();
equity_data.date = date;
equity_data.se = equity.se;
equity_data.sd = equity.sd;
equity_data.s = equity.s;

%%% value weighted dividend premium
% 
% data = readtable('Data/Input/Baker Wurgler/firm_data.csv');
% data = readtable('Data/Input/Baker Wurgler/m2b.csv');
% data = readtable('Data/Input/Baker Wurgler/m2b2.csv');
% data = readtable('Data/Input/Baker Wurgler/firm_data3.csv');
% data = readtable('Data/Input/Baker Wurgler/firm_data2.csv');

data = readtable('Data/Input/Baker Wurgler/firm_data5.csv');

%%% processing data
data.date = datenum(char(string(data.DATE)), 'yyyymmdd');
ind = (data.sic > 4899).*(data.sic< 4950);
data = data(ind == 0,:);
ind = (data.sic>5999).*( data.sic < 7000);
data = data(ind == 0,:);

%%% creating data for output
data = array2table([data.date data.PERMCO abs(data.mktvl) abs(data.BE) ...
    abs(data.BA) abs(data.dvt) abs(data.lag_dvt)]);
data.Properties.VariableNames = {'date', 'permco', 'mktvl', 'BE', 'BA','dvt', 'lag_dvt'};

%%% market to book measure
data.BE = data.BE.*1000000; %converting to dollar values
data.BA = data.BA.*1000000; %converting to dollar values
data.mktvl = data.mktvl * 1000; %converting to dollar values
data.mb = data.mktvl./data.BE;

%%% restricting values such that in 2009 dollars these firms have greater
%%% than $250,000 in book equity and assets greater than $500,000

fred = fred.latest({'CPIAUCSL'});
ind = month(fred.date) == 1;
cpi_date = fred.date(ind == 1,:);
cpi = fred.value(ind == 1,:);

for i = 1:length(cpi_date)
    ind = (year(data.date) == year(cpi_date(i)));
    data.BE_cpi(ind == 1) = data.BE(ind == 1)*cpi(50)/cpi(i);
    data.BA_cpi(ind == 1) = data.BA(ind == 1)*cpi(50)/cpi(i);
end

ind = (data.BE_cpi > 250000).*(data.BA_cpi > 500000);
data = data(ind == 1,:);

%%% dropping data with nan dividend values
ind = isnan(data.dvt);
data = data(ind == 0,:);

%eliminating observations with no marketvalue, bookvalue, or dividend data
ind = (isnan(data.mb) ==0).*(isnan(data.BE) == 0).*(data.BE>0);
data = data(ind == 1,:);

% Determening Dividend payers
ind = data.dvt > 0;
div = data(ind == 1,:);
no_div = data(ind==0,:);

% determening market to book by month for div and no_div
dates = unique(data.date);
ewm2b = zeros(size(dates,1),2); %column 1 is div, column 2 is no div
vwm2b = zeros(size(dates,1),2);
for i = 1:length(dates)
    ind = ismember(div.date, dates(i)); 
    ewm2b(i,1) = mean(div.mb(ind == 1));
    vwm2b(i,1) = sum(div.mb(ind ==1).*div.BE(ind == 1));
    vwm2b(i,1) = vwm2b(i,1)./(sum(div.BE(ind == 1)));
    ind = ismember(no_div.date, dates(i)); 
    ewm2b(i,2) = mean(no_div.mb(ind == 1));
    vwm2b(i,2) = sum(no_div.mb(ind ==1).*no_div.BE(ind==1));
    vwm2b(i,2) = vwm2b(i,2)./(sum(no_div.BE(ind == 1)));
end

ewm2b = log(ewm2b);
vwm2b = log(vwm2b); 
ewm2b(:,3) = 100*(ewm2b(:,1) - ewm2b(:,2));
vwm2b(:,3) = 100*(vwm2b(:,1) - vwm2b(:,2));
div_prem = vwm2b(:,3); 

%%% plotting to compare
plot(dates, [div_prem(:,1) bw.pdnd((height(bw) - length(dates)+1):end)])
corr([div_prem(:,1) bw.pdnd((height(bw) - length(dates)+1):end)], 'rows', 'complete')
datetick('x')
title('Dividend Premium Comparison')
ylabel('Dividend Premium')
legend('My Measure', 'BW Measure')

% saving data
pdnd = struct();
pdnd.date = dates;
pdnd.date = datenum(year(dates), month(dates)+1,1);
if month(dates) == 12
    pdnd.date = datenum(year(dates), 12, 31);
end
pdnd.date = pdnd.date- 1; 
pdnd.value = div_prem; 

%%% closed end fund discount

cefd = struct();
cefd.date = pdnd.date;
cefd.value = zeros(size(pdnd.date,1),1);

%% loading in historical data

bw = readtable('Data/Input/Baker Wurgler/bw_data.csv','TreatAsEmpty', {'.'});
bw.ripo = str2double(bw.ripo);
bw = bw(1:732,:);
datest = char(string(bw.x_yearmo));
bw_dates = datenum(str2num(datest(:,1:4)),str2num(datest(:,5:6))+1,1)-1;

clear con* debt* emp* end* fred gap header* ind* monthly* no* div_pay* ...
    equity_m* equity_y* recess nipo ripo ritter_ipo table_* years yearval 
clear x i equity_share date cpi ans a dividend* firm_data
clear div* dat* equity ew* vw* filename* month* o* start_* vw* y*  


%% using measures to update the BW spreadsheet

%%% setting up data to export
dates = union(fred_data.date, union(pdnd.date, union(equity_data.date, union(cefd.date, ipo_data.date)))); 
dates = dates(year(dates) > 1950);
data_out = struct();
variables = {'yearmo' 'pdnd' 'ripo' 'nipo' 'cefd' 'se' 'sd' 's' 'blank' 'indpro' 'consdur' 'consnon' 'consserv' 'recess' 'employ' 'cpi'};
for i = 1:length(variables)
    data_out.(variables{i}) = nan(size(dates,1),size(dates,2));
end

%%% filling in prexisting bw data
ind1 = ismember(bw_dates,dates);
ind2 = ismember(dates, bw_dates);
for i = 2:length(variables) 
    data_out.(variables{i})(ind2 == 1) = table2array(bw(ind1 == 1,i+2));
end

%%% replacing with updated data
ind1 = ismember(equity_data.date,dates);
ind2 = ismember(dates, equity_data.date);
data_out.se(ind2 == 1) = equity_data.se(ind1==1);
data_out.sd(ind2 == 1) = equity_data.sd(ind1==1);
data_out.s(ind2 == 1) = equity_data.s(ind1==1);

ind1 = ismember(fred_data.date,dates);
ind2 = ismember(dates, fred_data.date);
for i = 10:length(variables)
    data_out.(variables{i})(ind2 == 1) = fred_data.(variables{i})(ind1 == 1);
end

ind1 = ismember(ipo_data.date,dates);
ind2 = ismember(dates, ipo_data.date);
data_out.ripo(ind2 == 1) = ipo_data.ripo(ind1==1);
data_out.nipo(ind2 == 1) = ipo_data.nipo(ind1==1);
 
% ind1 = ismember(cefd.date,dates);
% ind2 = ismember(dates, cefd.date);
% data_out.cefd(ind2 == 1) = cefd.value(ind1==1);

ind1 = ismember(pdnd.date,dates);
ind2 = ismember(dates, pdnd.date);
data_out.pdnd(ind2 == 1) = pdnd.value(ind1==1);

yearmo = [year(dates) month(dates)];
yearmo = num2str(yearmo,'%02d');
yearmo = [yearmo(:,1:4) yearmo(:,7:8)];
yearmo = str2num(yearmo);
data_out.yearmo = yearmo;

export_to_csv = struct2table(data_out);
filename = 'Data/Input/Baker Wurgler/bw_data_input.csv';
writetable(export_to_csv, filename);


%% calling stata code on excel sheet to create sentiment measures

cd(stata_dir)

!/apps/stata16/stata-se -b "replicate_bw_sentiment.do"

cd(root_dir)





