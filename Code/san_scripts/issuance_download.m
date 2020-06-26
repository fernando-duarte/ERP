% This file downloads equity and debt issuance

clear,clc,close all;
%% Declare Path Directories

if isunix
    proj_dir = '/san/RDS/Work/cmf/b1crs02/';
else
    proj_dir = '\\rb\b1\NYRESAN\RDS\Work\cmf\b1crs02';
end

cd(proj_dir);

%% Downloading Data

url = 'https://fweb.rsma.frb.gov/cgi-bin/fame/fame.pl?output=CSV&bdate=%27*%27&edate=%27*%27&deci=Auto&transform=level&freq=%27no%20conversion%27&selection=US%27DTBCFBNMD_N.M&display_series=YES';
outfile = [proj_dir filesep 'Projects' filesep 'etp_personal' filesep ...
    'issuance_download' filesep 'debt.csv'];
websave(outfile, url);

url = 'https://fweb.rsma.frb.gov/cgi-bin/fame/fame.pl?output=CSV&bdate=%27*%27&edate=%27*%27&deci=Auto&transform=level&freq=%27no%20conversion%27&selection=US%27EQFSNBCRSS_N.M&display_series=YES';
outfile = [proj_dir filesep 'Projects' filesep 'etp_personal' filesep ...
    'issuance_download' filesep 'equity.csv'];
websave(outfile, url);