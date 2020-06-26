%%% This file computes CAY measure as produced by Lettau and Ludvigson in
%%% 2001. 

%% Setting Directory
dir = [root_dir filesep 'Input'];

%% Downloading Data

%%% fred data for C and A
fred = fred.latest({'A794RX0Q048SBEA', 'A229RX0Q048SBEA'});
fred_label = {'PCE', 'RDI'};

str = datestr(datenum(datetime('today')),'yyyy-mm-dd');
url = ['https://fred.stlouisfed.org/graph/fredgraph.csv?bgcolor=%23e1e9f0&' ...
    'chart_type=line&drp=0&fo=open%20sans&graph_bgcolor=%23ffffff&height=4' ...
    '50&mode=fred&recession_bars=on&txtcolor=%23444444&ts=12&tts=12&width=1' ...
    '168&nt=0&thu=0&trc=0&show_legend=yes&show_axis_titles=yes&show_tooltip' ...
    '=yes&id=HNONWPDPI&scale=left&cosd=1946-10-01&coed=' str  ...
    '&line_color=%234572a7&link_values=false&line_style=solid&mark_type=' ...
    'none&mw=3&lw=2&ost=-99999&oet=99999&mma=0&fml=a&fq=Quarterly%2C%20End' ...
    '%20of%20Period&fam=avg&fgst=lin&fgsnd=2009-06-01&line_index=1&trans' ...
    'formation=lin&vintage_date=' str '&revision_date=' str '&nd=1946-10-01'];
filename = [dir filesep 'hhnw_data.csv'];
outfile = websave(filename,url);

hhnw = readtable(filename);
hhnw_date = datenum(hhnw.DATE);
hhnw_data = str2double(table2array(hhnw(:,2)));
ind = hhnw_date > datenum(1951,12,31);
hhnw_date = hhnw_date(ind ==1);
hhnw_data = hhnw_data(ind ==1);

%%% bea data to compute Y 
filename =  [dir filesep 'rdic_data.xlsx'];
url = 'https://apps.bea.gov/histdata/Releases/GDP_and_PI/2019/Q4/Advance_January-31-2020/Section2all_xls.xlsx';
outfile = websave(filename,url, 'Timeout', 60);

rdic = table2array(readtable(filename,'Sheet',3));
rdic = rdic(7:end,24:end);
rdic_date = char(string(rdic(1,:)'));
rdic_year = str2num(rdic_date(:,1:4));
rdic_quarters = str2num(rdic_date(:,6));
rdic_date = datenum(rdic_year, 3*rdic_quarters-2,1);
rdic_data = 1000*(sum([str2num(char(rdic(3,:))),str2num(char(rdic(17,:)))],2) - sum([str2num(char(rdic(26,:))), str2num(char(rdic(27,:))).*str2num(char(rdic(3,:)))./(str2num(char(rdic(3,:)))+str2num(char(rdic(10,:)))+str2num(char(rdic(13,:)))+str2num(char(rdic(14,:))))],2))./str2num(char(rdic(44,:)));
rdic_data = rdic_data.*str2num(char(rdic(43,:)))./str2num(char(rdic(42,:)));

%% Matching Dates

date = intersect(fred.date, intersect(hhnw_date, rdic_date));
ind = date > datenum(1951,12,31); 
date = date(ind == 1); 

ind = ismember(fred.date, date);
fred.value = fred.value(ind == 1,:);

ind = ismember(hhnw_date, date);
hhnw_data = hhnw_data(ind == 1,:);

ind = ismember(rdic_date, date);
rdic_data = rdic_data(ind == 1,:);

%% Determening C, A, and Y

c = fred.value(:,1);
a = hhnw_data.*fred.value(:,2)/100;
y = rdic_data;

%%% taking logs and converting to timetable
c = log(c);
a = log(a);
y = log(y);

dateindex = datetime(date, 'ConvertFrom','datenum');
A = timetable(dateindex,a);
C = timetable(dateindex,c);
Y = timetable(dateindex,y);

%% Running regressions to compute CAY

k = 8; %number of lags

Xbeta = [table2array(A,2) table2array(Y,2)];
fixedA = zeros(size(Xbeta,1),k*2);
for i = 1:k*2
    fixedA(:,i) = table2array(lag(A,i-k),2) - table2array(lag(A,i-1-k),2);
end
fixedY = zeros(size(Xbeta,1),k*2);
for i = 1:k*2
    fixedY(:,i) = table2array(lag(Y,i-k),2) - table2array(lag(Y,i-1-k),2);
end
Xfixed = [fixedA fixedY];

X = [ones(size(Xbeta,1),1) Xbeta Xfixed];

beta = regress(c(k+1:end,1),X(k+1:end,:));

cay_value = c - beta(2,1) * a - beta(3,1) * y - beta(1,1);

%% Exporting data to csv

export_to_excel = array2table([date c a y cay_value]);
export_to_excel.Properties.VariableNames = {'dates', 'c', 'a', 'y', 'cay'};
export_to_excel.dates = datestr(export_to_excel.dates);
filename = [dir filesep 'cay_updated.csv']; 
writetable(export_to_excel, filename); 

clearvars -except root_dir

cd(root_dir)

disp('Finished cay.m')
