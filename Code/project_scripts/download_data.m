%% download_data

%%% This file downloads data required to compute different ERP measures. It
%%% uses WRDS/Compustat data, FRED data, data from Shiller, the FED, and
%%% other websites. 

%% Download WRDS Data

%%% pathing in Sas Studio WRDS
wrds_code_dir = '/home/frb-ny/fedraj/WRDS/Code';
wrds_out_dir = '/home/frb-ny/fedraj/WRDS/Output';

%%% Connect to the WRDS server
wrds_username = 'fedraj';
wrds_password = 'Feddata2020$%';
a.SSH2conn = ssh2_config('wrds.wharton.upenn.edu',wrds_username,wrds_password,22);

% Update the SAS code on the server
scp_put(a.SSH2conn,'wrdspull.sas', wrds_code_dir,'Code/wrds_scripts/','wrdspull.sas');
scp_put(a.SSH2conn,'get_crspm.sas',wrds_code_dir,'Code/wrds_scripts/','get_crspm.sas');
scp_put(a.SSH2conn,'get_cmt.sas',wrds_code_dir,'Code/wrds_scripts/','get_cmt.sas');

% Remotely run the SAS code and download the data files to local folders
ssh2_command(a.SSH2conn, ['qsas ' wrds_code_dir '/wrdspull.sas'])
ssh2_command(a.SSH2conn, ['qsas ' wrds_code_dir '/get_crspm.sas'])
ssh2_command(a.SSH2conn, ['qsas ' wrds_code_dir '/get_cmt.sas'])

% Transfer data to appropriate folders. 
scp_get(a.SSH2conn,'sp500_book_to_market.xlsx', 'Input/', wrds_out_dir);
scp_get(a.SSH2conn,'EPS_estimates.csv', 'Input/', wrds_out_dir);
scp_get(a.SSH2conn,'cmt_data.csv', 'Input/',wrds_out_dir);
try
    scp_get(a.SSH2conn,'crsp_m.txt', 'Input/',wrds_out_dir); %download works but causes an error sometimes
catch
end

% close connection
ssh2_close(a.SSH2conn);

%% Download assorted data from the internet
mkdir([root_dir filesep 'Code' filesep 'python_scripts' filesep 'pyData']) 

% downloading data from many different websites through python
py_path = [root_dir filesep 'Code' filesep 'python_scripts'];
cd(py_path)
!/apps/Anaconda3-2019.03/bin/python publicDataGrabber.py

% copying files from one folder to another
path = [py_path filesep 'pyData'];
cd(path)
movefile('*', [root_dir filesep 'Input']);
cd(root_dir)

% creating certificate for web acess 
o = weboptions('CertificateFilename',"");

% downloading other data
url = 'https://us.spindices.com/documents/additional-material/sp-500-eps-est.xlsx';
filename = [root_dir filesep 'Input' filesep 'sp-500-eps-est_download.xlsx'];
outfile = websave(filename,url, o);

url = 'http://www.federalreserve.gov/data/yield-curve-tables/feds200628.csv';
filename = [root_dir filesep 'Input' filesep 'feds200628.csv'];
outfile = websave(filename,url, o);

url = 'http://people.stern.nyu.edu/jwurgler/data/Investor_Sentiment_Data_20190327_POST.xlsx';
filename = [root_dir filesep 'Input' filesep 'Investor_Sentiment_data_20190327.xlsx'];
outfile = websave(filename,url, 'Timeout', 60);

dates = datestr(datenum(datetime('today')),'yyyy-mm-dd');
url = ['https://fred.stlouisfed.org/graph/fredgraph.csv?' ... 
    'bgcolor=%23e1e9f0&chart_type=line&drp=0&fo=open%20sans&' ... 
    'graph_bgcolor=%23ffffff&height=450&mode=fred&recession_bars=on' ... 
    '&txtcolor=%23444444&ts=12&tts=12&width=1168&nt=0&thu=0&trc=0&' ... 
    'show_legend=yes&show_axis_titles=yes&show_tooltip=yes&id=DGS10&'...
    'scale=left&cosd=1962-01-02&coed=' dates '&line_color=%234572a7' ... 
    '&link_values=false&line_style=solid&mark_type=none&mw=3&lw=2&ost=' ...
    '-99999&oet=99999&mma=0&fml=a&fq=Daily&fam=avg&fgst=lin&fgsnd=' ...
    '2009-06-01&line_index=1&transformation=lin&vintage_date=' ...
    dates '&revision_date=' dates '&nd=1962-01-02'];
filename = [root_dir filesep 'Input' filesep 'tenyearyields.csv'];
outfile = websave(filename,url, o);

%% Replicate Fama French momentum portfolio 
stata_path = [root_dir filesep 'Code' filesep 'stata_scripts'];
cd(stata_path)
!/apps/stata16/stata-se -b "replicate_mom_portfolios.do"
cd(root_dir)

%% Download FRED Data

% Download TIPS data from FRED
fred = fred.latest({'FII10','FII20','FII30','FII5','FII7'});
DATE = datetime(fred.date,'ConvertFrom','datenum');
DATE = datestr(DATE,'yyyy-mm-dd');
data = fred.value;
FII10 = data(:,1); FII20 = data(:,2); FII30 = data(:,3); FII5 = data(:,4);
FII7 = data(:,5);
T = table(DATE,FII10,FII20,FII30,FII5,FII7);
filename = [root_dir filesep 'Input' filesep 'TIPS.csv'];
writetable(T,filename); 
clear fred T DATE data 

% Download BAA and AAA bond yields from FRED
fred = fred.latest({'BAA','AAA'});
observation_date = datetime(fred.date,'ConvertFrom','datenum');
observation_date = datestr(observation_date,'yyyy-mm-dd');
data = fred.value;
BAA_AAA = data(:,1)-data(:,2);
T = table(observation_date,BAA_AAA);
filename = [root_dir filesep 'Input' filesep 'Default_spread.csv'];
writetable(T,filename); 
clear fred T 

% Download 3-month and 6-month T-Bill yields
fred = fred.latest({'TB3MS','TB6MS'});
DATE = datetime(fred.date,'ConvertFrom','datenum');
DATE = datestr(DATE,'yyyy-mm-dd');
data = fred.value;
TB3MS = data(:,1); TB6MS = data(:,2);
T = table(DATE,TB3MS,TB6MS);
filename = [root_dir filesep 'Input' filesep 'Short_term_yields.csv'];
writetable(T,filename); 
clear fred T

% recessions
fred = fred.latest({'USREC'});
data = array2table([fred.date fred.value]);
data.Properties.VariableNames = {'date', 'recession'};
data.date = datestr(data.date);
filename = [root_dir filesep 'Input' filesep 'recession.csv'];
writetable(data,filename);

%% End

clearvars -except root_dir
disp('Finished download_data.m')