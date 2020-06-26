%% NY_Fed_CM_DAPM

%%% This file recreates the NY FED ERP measure as produced in Adrian,
%%% Crump and Moench.

%% Set directories
cd(root_dir)

clearvars -except root_dir

dir = [root_dir filesep 'Input'];
cd(dir)

%% Downloading Data from Web

%%% Fama French
url = 'http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_Factors_TXT.zip';
filename = [dir filesep 'F-F--Research_Data_Factors.zip'];
outfile = websave(filename, url);
unzip(filename);

url = 'http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/Portfolios_Formed_on_ME_CSV.zip';
filename = [dir filesep 'Portfolios_Formed_on_ME.zip'];
outfile = websave(filename, url);
unzip(filename);

cd(root_dir)

%%% Shiller
filename = ['ie_data.xls'];
shiller = readtable(filename,'Sheet',4, 'TreatAsEmpty', {'NA', '.', ''});
shiller = shiller(7:end-1, 1:13);
shiller.Properties.VariableNames = {'date', 'P', 'D', 'E', 'CPI', 'date2', ...
    'RateGS10', 'RealPrice', 'RealDividend', 'RealTotalReturnPrice', ...
    'RealEarnings', 'RealEarningsScaled', 'CAPE'};

% deleting unecessary data
shiller.date2 = [];
shiller.RealTotalReturnPrice = [];
shiller.RealEarningsScaled = [];
ind = strcmp(shiller.date, {''});
shiller = shiller(ind ==0,:);

%creating dates
dates = char(table2array(shiller(:,1)));
year = dates(:,1:4);
month = dates(:,6:7);
ind = strcmp(month(:,2),{''});
month = str2num(month);
month(ind == 1,:) = 10;
shiller_dates = datenum(strcat(year, '-', string(month), '-', '01'), 'yyyy-mm-dd');

for i = 2:size(shiller,2)
    final = zeros(size(shiller,1),1);
    temp = cellfun(@str2num,shiller{:,i}, 'UniformOutput', 0);
    ind = strcmp(shiller{:,i}, {''});
    ind2 = strcmp(shiller{:,i}, {'NA'});
    ind3 = (ind == 0).*(ind2==0);
    final(ind3 == 1,:) = cell2mat(temp);
    final(ind == 1,:) = NaN;
    final(ind2 == 1,:) = NaN;
    data(:,i-1) = array2table(final);
end
data.date = dates;
data.Properties.VariableNames = {'P', 'D', 'E', 'CPI', ...
    'RateGS10', 'RealPrice', 'RealDividend', ...
    'RealEarnings', 'CAPE', 'date'};


clear year month

%%% gs30
url = 'https://fred.stlouisfed.org/graph/fredgraph.csv?bgcolor=%23e1e9f0&chart_type=line&drp=0&fo=open%20sans&graph_bgcolor=%23ffffff&height=450&mode=fred&recession_bars=on&txtcolor=%23444444&ts=12&tts=12&width=1168&nt=0&thu=0&trc=0&show_legend=yes&show_axis_titles=yes&show_tooltip=yes&id=GS30&scale=left&cosd=1977-02-01&coed=2019-09-01&line_color=%234572a7&link_values=false&line_style=solid&mark_type=none&mw=3&lw=2&ost=-99999&oet=99999&mma=0&fml=a&fq=Monthly&fam=avg&fgst=lin&fgsnd=2009-06-01&line_index=1&transformation=lin&vintage_date=2019-10-28&revision_date=2019-10-28&nd=1977-02-01';
filename = [dir filesep 'gs30.csv'];
outfile = websave(filename,url);

gs30 = readtable(filename);
gs30_dates = datenum(gs30.DATE);

%%% gs20
url = 'https://fred.stlouisfed.org/graph/fredgraph.csv?bgcolor=%23e1e9f0&chart_type=line&drp=0&fo=open%20sans&graph_bgcolor=%23ffffff&height=450&mode=fred&recession_bars=on&txtcolor=%23444444&ts=12&tts=12&width=1168&nt=0&thu=0&trc=0&show_legend=yes&show_axis_titles=yes&show_tooltip=yes&id=GS20&scale=left&c';
filename = [dir filesep 'gs20.csv'];
outfile = websave(filename,url);

gs20 = readtable(filename);
gs20_dates= datenum(gs20.DATE);

%%% tb1m
tb1m = readtable('tb1m.csv');
tb1m_date = datenum(tb1m.date);
tb1m_value = table2array(tb1m(:,1));

%%% FRED 
fred_data = fred.latest({'GS1','GS3','GS5','GS7', 'GS10' ... 
    'BAA', 'CPIAUCSL', 'CPILFESL', 'TB3MS', 'PAYEMS'});

%% Merging Downloaded Data

fred_dates = fred_data.date;
dates = union(fred_dates, union(gs20_dates, ...
    union(gs30_dates, union(tb1m_date,shiller_dates))));

ind = ismember(dates,fred_dates);
fred_values(ind == 0,1:10) = NaN;
fred_values(ind == 1,:) = fred_data.value;

x = 100*diff(fred_values(:,8))./fred_values(1:end-1,8);
fred_values(2:end,8) = x;

ind = ismember(dates,gs30_dates);
gs30_data(ind == 0,:) = NaN;
gs30_data(ind == 1,:) = str2double(gs30.GS30);

ind = ismember(dates,gs20_dates);
gs20_data(ind == 0,:) = NaN;
gs20_data(ind == 1,:) = str2double(gs20.GS20);

ind = ismember(dates,tb1m_date);
tb1m_data(ind == 0,:) = NaN;
tb1m_data(ind == 1,:) = tb1m_value;

ind = ismember(dates,shiller_dates);
SDY5(ind == 0,:) = NaN;
SDY5(ind == 1,:) = data.D./data.P;
SPECAPE(ind == 0,:) = NaN;
SPECAPE(ind == 1,:) = data.CAPE;
SP500(ind == 0,:) = NaN;
SP500(ind == 1,:) = data.P;

data = array2table([dates tb1m_data fred_values gs20_data gs30_data ...
    SDY5 SPECAPE SP500]);
data.Properties.VariableNames = {'date', 'GS1TR', 'GS1' 'GS3' 'GS5', 'GS7', ...
    'GS10', 'BAA', 'CPIAUCSL', 'CPILFESL', 'TB3MS', 'PAYEMS', 'GS20', ...
    'GS30', 'SDY5', 'SPECAPE','SP500'};

%%% reorganizing variables

data = data(:,[1:7 13:14 8:12 15 17 16]);

clear dates filename fred* gs* ind s* S* url x i 

%% Cleaning Data

%%% Fama French

startrow = 13;
lastupdate = datevec(datenum(2019,9,19)); 
today_date = datenum(datetime('today'));
date_vec = datevec(today_date);
endrow = 1129 + 12*(date_vec(1) - lastupdate(1)) + date_vec(2) - lastupdate(2);
size = dlmread('Portfolios_Formed_on_ME.CSV', ',',[startrow, 0, endrow,19]);
ind = isnan(size(:,1));
size(ind == 1,:) = [];
size(:,1) = datenum(num2str(size(:,1)),'yyyymm');
size = array2table(size);

clear ind

Research_Factors = readtable('F-F_Research_Data_Factors.txt', 'ReadVariableNames', false);
startrow = 3;
endrow = 1119 + 12*(date_vec(1) - lastupdate(1)) + date_vec(2) - lastupdate(2);
Research_Factors_final = Research_Factors(startrow:endrow,1:5);
Research_Factors_final.Properties.VariableNames =  {'date', 'MktRF','SMB','HML','RF'};
Research_Factors_final.date = datenum(string(Research_Factors_final.date),'yyyymm');
Research_Factors_final.MktRF = str2num(char(Research_Factors_final.MktRF));
Research_Factors_final.SMB = str2num(char(Research_Factors_final.SMB));
Research_Factors_final.HML = str2num(char(Research_Factors_final.HML));
Research_Factors_final.RF = str2num(char(Research_Factors_final.RF));
Research_Factors_final.MKT = Research_Factors_final.MktRF + Research_Factors_final.RF;

% wrds cmt data

wrds = readtable('cmt_data.csv');
cmt = unstack(wrds, {'value'}, {'cmtcode'});
cmt.Properties.VariableNames = {'date', 'cmt2', 'cmt1', 'cmt5', 'cmt7', ...
    'cmt10', 'cmt20', 'cmt30'};
date = datetime(datenum(cmt.date),'ConvertFrom', 'datenum');
date = datenum(date.Year, date.Month,1);
cmt.date = date;  

% saving the data

final_cmt = cmt;
final_size = [size(:,1),size(:,11:20)];
final_size.Properties.VariableNames{'size1'} = 'date';
final_FF = Research_Factors_final;
final_haver = data;

clear ans dates fbaa FCM*  ftbs3e lanagra merge mpcuslfe pcun cmt* lastupdate
clear Research*  sdy5 size SP500E SPECAE tb1m wrds a SPECAPE endrow url data 
clear final_data tb1m* startrow 

%% Merging the Different Datasets to match dates

final_data = outerjoin(final_FF,final_cmt, 'MergeKeys', 1);
final_data = outerjoin(final_data,final_size, 'MergeKeys', 1);
final_data = outerjoin(final_data,final_haver, 'MergeKeys', 1);

clear final_cmt final_FF final_haver final_size

final_data = table2array(final_data);

final_FF = [final_data(:,1),final_data(:,2:6)];
final_size = [final_data(:,1),final_data(:,14:23)];
final_cmt = [final_data(:,1),final_data(:,7:13)];
final_haver = [final_data(:,1),final_data(:,24:end)];

%% Converting to Structs

FF = struct();
FF.data = final_FF(:,2:end);
FF.date = final_FF(:,1);
FF.names = {'MKTRF', 'SMB', 'HML', 'RF', 'MKT'};

RF = struct();
RF.data = final_FF(:,5);
RF.date = FF.date;
RF.names = {'RF'};
clear Research_Factors* 

size10 = struct();
size10.data = final_size(:,2:end);
size10.date = final_size(:,1);
size10.names = {'size1','size2','size3','size4','size5','size6','size7', ...
    'size8','size','size10'};

clear size

haver = struct();
haver.data = final_haver(:,2:end);
haver.date = final_haver(:,1);
haver.names = {'FTBTR1M', 'FCM1E','FCM3E','FCM5E','FCM7E','FCM10E','FCM20E', ...
    'FCM30E','FBAAE','SDY5','PCUN', 'MPCUSLFE','FTBS3E','LANAGR','SP500E','SPECAPE',};

cmts = struct();
cmts.data = final_cmt(:,2:end);
cmts.date = final_cmt(:,1);
cmts.names = {'cmt1','cmt2','cmt5','cmt7','cmt10','cmt20','cmt30'};

clear final*

%% Creating Data for ERP Calculation

mkt = FF.data(:,end);
term = haver.data(:,6) - haver.data(:,12);
term_3alt = haver.data(:,6) - haver.data(:,3);
term_alt = haver.data(:,4) - haver.data(:,12);
TSY10 = haver.data(:,6);
TSY5 = haver.data(:,4);
TSY3M = haver.data(:,12);
CINF = haver.data(:,11);
HML = FF.data(:,3); 
SP500_HAV = nan(size(haver.data,1),1);
for i = 2:size(haver.date,1)
    SP500_HAV(i,1) = 100*((haver.data(i,15) ./ (haver.data(i-1,15))) - 1);
end
SMB = FF.data(:,2); 
PAYROLL = haver.data(:,13);
NFP_A = nan(size(haver.data,1),1);
for i = 11:length(PAYROLL)
    NFP_A(i) = 100*log(PAYROLL(i,1)/PAYROLL(i-10,1));
end 
CP = nan(size(haver.data,1),1);
dy2 = log(haver.data(:,14));
cape = log(haver.data(:,end)/100);
DivEarn = cape - dy2;

%% Converting Data to Struct

MonthlyFactors = struct();
MonthlyFactors.data = [mkt, term, term_3alt, term_alt, TSY10, TSY5, TSY3M ...
    CINF, HML, SP500_HAV, SMB, PAYROLL, NFP_A, CP, dy2, cape, DivEarn];
MonthlyFactors.date = haver.date;
MonthlyFactors.names = {'MKT','TERM','TERM_3alt','TERM_alt','TSY10','TSY5', ...
    'TSY3M','CINF','HML','SP500_HAV','SMB','PAYROLL','NFP_A','CP','dy2', 'cape', 'DivEarn'};

clear CINF HML TERM* CP dy2 cape Div* TSY* SMB SP5* term* size* mkt i FF 
clear haver NFP* PAYROLL RF cmts

%% Computing Equity Risk Premium

%dates
Tstart_sample = datenum('196401', 'yyyymm');
Tstart_is = datenum('196401', 'yyyymm');
Tend = MonthlyFactors.date(end); 
currentdate = '20140826'; 

%loading in data

facs = MonthlyFactors.data;
dates = MonthlyFactors.date;
facstr = MonthlyFactors.names';

T0 = find(dates == Tstart_sample);
T1_is = find(dates == Tstart_is);
T2 = find(dates == Tend);

%whether to use Ken French data
IndUseKenFrench = false;

%prepare data 

X1list = char('MKT','SMB');
X2list = char('TSY10');
X3list = char('TERM','dy2');

X1full = facs(:, cell2mat(strindx_vec(facstr,X1list,'exact')));
X2full = facs(:, cell2mat(strindx_vec(facstr,X2list,'exact')));
X3full = facs(:, cell2mat(strindx_vec(facstr,X3list,'exact')));

faclist = char(X1list,X2list,X3list);

%state returns in basis point per month

X1full = X1full/100;
X2full = X2full/100;
X3full(:,1) = X3full(:,1)/100;    

%%% in sample analysis

if ~isempty(X1full), X1 = X1full(T1_is:T2,:); else X1 = []; end     
if ~isempty(X2full), X2 = X2full(T1_is:T2,:); else X2 = []; end 
if ~isempty(X3full), X3 = X3full(T1_is:T2,:); else X3 = []; end

%%% deleting NaNs
ind = sum(isnan([X1 X2 X3]),2);

X1 = X1(ind == 0,:);
X2 = X2(ind == 0,:);
X3 = X3(ind == 0,:);
final_dates = dates(T1_is:T2);
final_dates = final_dates(ind == 0,:);

%%% running analysis

[T,k1] = size(X1);
k2 = size(X2,2);
k3 = size(X3,2);
kF = k2+k3;
k = k1+k2+k3;


load(['DAPM_params_', currentdate, '.mat']);

horsel = 1:120;
maxh = max(horsel);
F_tilde_fcsts = zeros(T,1+kF,maxh);

if IndUseKenFrench
    mu = outputs.Psi_CB(1,:)'; Phi = outputs.Psi_CB(2:end,:)';
    mu_fcsts = zeros(k,maxh);
    mu_fcsts(:,1) = mu;
    Phi_fcsts = zeros(k,k,maxh);
    X_fcsts = zeros(T,k,maxh);
    X = [X1 X2 X3];
    Z = [ones(T,1) X];
    for h = 1:maxh
        if h > 2
            mu_fcsts(:,h) = [eye(k) + sum(Phi_fcsts(:,:,1:h-1),3)]*mu;
        end
        Phi_fcsts(:,:,h) = Phi^h;
        X_fcsts(:,:,h) = Z*[mu_fcsts(:,h) Phi_fcsts(:,:,h)]';
        F_tilde_fcsts(:,:,h) = [ones(T,1) X_fcsts(:,k1+1:end,h)];    
    end
else
    X_QS = [X2(2:end,:) X3(2:end,:)];
    Z_QS = [ones(size(X2,1)-1,1) X2(1:end-1,:) X3(1:end-1,:)];
    Psi_QS = ((Z_QS'*Z_QS)\Z_QS'*X_QS)';
    mu = Psi_QS(:,1);
    Phi = Psi_QS(:,2:end);
    %
    mu_fcsts = zeros(kF,maxh);
    mu_fcsts(:,1) = mu;
    Phi_fcsts = zeros(kF,kF,maxh);
    X_fcsts = zeros(T,kF,maxh);
    X = [X2 X3];
    Z = [ones(T,1) X];
    for h = 1:maxh
        if h > 2
            mu_fcsts(:,h) = [eye(kF) + sum(Phi_fcsts(:,:,1:h-1),3)]*mu;
        end
        Phi_fcsts(:,:,h) = Phi^h;
        X_fcsts(:,:,h) = Z*[mu_fcsts(:,h) Phi_fcsts(:,:,h)]';
        F_tilde_fcsts(:,:,h) = [ones(T,1) X_fcsts(:,:,h)];    
    end
end

for hor = 1:length(horsel)
    h = horsel(hor);
    Premth_MKT{hor} = ones(T,1);    
    if (hor==1)
        Premth_MKT_cont = ones(T,1+kF);
        for i = 1:(1+kF)
            Lamt_h_i = F_tilde_fcsts(:,i,h)*outputs.Lam_CB(1,i);
            Premth_MKT_cont(:,i)= Premth_MKT_cont(:,i).*(ones(T,1) + Lamt_h_i);
        end
    end
    for hh = 1:h
        Lamt_h = F_tilde_fcsts(:,:,hh)*outputs.Lam_CB';
        Premth_MKT{h} = Premth_MKT{h}.*(ones(T,1) + Lamt_h(:,1));
    end
    Premth_MKT_ann{hor} = 100*(Premth_MKT{hor}.^(12/h)-1); % produce annual rates of return
end

horse1_QS = [1 12 60 120];
horse1str = {'1 month ERP', '1 year ERP', '5 year ERP', '10 year ERP'};
Premth_MKT_QS = cell2mat(Premth_MKT_ann);

Fed_CM_ERP = Premth_MKT_QS(:,horse1_QS);

%%% saving output
final = [final_dates, Fed_CM_ERP];
filename = [root_dir filesep 'Input' filesep 'FRBNY_ERP_QS_crs.mat'];
save(filename, 'final') 

cd(root_dir)

clearvars -except root_dir

disp('Finished NY_Fed_CM_DAPM.m')