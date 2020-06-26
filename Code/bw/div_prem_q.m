clear 
clc

%%% loading in data
CRSP = readtable("/Volumes/rcecrs02/ERPv3/CRSP.csv");
COMPU = readtable("/Volumes/rcecrs02/ERPv3/CRSP_COMPUQ.csv");

%%% cleaning
CRSP.date = datenum(char(string(CRSP.date)), 'yyyymmdd');
COMPU.date = datenum(char(string(COMPU.datadate)), 'yyyymmdd');

CRSP.qtr = quarter(CRSP.date);
CRSP.year = year(CRSP.date);
COMPU.qtr = COMPU.fqtr;
COMPU.year = year(COMPU.date); 
ind = COMPU.qtr == 1;
COMPU.year(ind) = COMPU.year(ind) - 1;
COMPU.qtr(ind) = 4;
COMPU.qtr(ind == 0) = COMPU.qtr(ind == 0) - 1;
COMPU.PERMNO = COMPU.LPERMNO;
COMPU.PERMCO = COMPU.LPERMCO;

%%% merging by PERMNO and quarter
data = innerjoin(CRSP,COMPU, 'keys', {'year', 'qtr', 'PERMNO'});
data = sortrows(data, [1:2]);

%%% creating variables
data.div_COMPU = data.dvpspq;
data.div_CRSP = (data.vwretd - data.vwretx);
data.txditcq(isnan(data.txditcq)) = 0;

data.SHE = data.seqq;
ind = isnan(data.SHE);
data.SHE(ind) = data.ceqq(ind) + data.pstkq(ind);
ind = isnan(data.SHE);
data.SHE(ind) = data.atq(ind) - data.ltq(ind) - data.mibq(ind);
data.SHE = data.SHE + data.txditcq;
data.PS = data.pstkrq;
ind = isnan(data.PS);
data.PS(ind) = data.pstkq(ind);

data.BE = (data.SHE - data.PS);
data.BE = data.BE * 1000000;
data.BA = data.atq*1000000;
data.COMPU_mktvl = data.mkvaltq;
data.CRSP_mktvl = abs(data.PRC).*data.SHROUT*1000; 
data.mb = data.CRSP_mktvl./data.BE;

%%% deleting duplicates
data = unique(data, 'stable'); 

%%% dropping data in wrong industry
ind = (data.sic > 4899).*(data.sic< 4950);
data = data(ind == 0,:);
ind = (data.sic>5999).*( data.sic < 7000);
data = data(ind == 0,:);

%%% dropping data if too small book value or book assets
ind = (data.BE > 250000).*(data.BA > 500000);
data = data(ind == 1,:);

%%% dropping data if missing key variables
ind = (isnan(data.mb) ==0).*(isnan(data.BE) == 0).*(data.BE>0).*(isnan(data.div_COMPU)==0);
data = data(ind == 1,:);

% Determening Dividend payers
ind = data.div_COMPU > 0;
div = data(ind == 1,:);
ind = data.div_COMPU == 0;
no_div = data(ind,:);

% determening market to book by month for div and no_div
dates = unique(data.date_CRSP);
ewm2b = zeros(size(dates,1),2); %column 1 is div, column 2 is no div
vwm2b = zeros(size(dates,1),2);
for i = 1:length(dates)
    ind = ismember(div.date_CRSP, dates(i)); 
    ewm2b(i,1) = mean(div.mb(ind == 1));
    vwm2b(i,1) = sum(div.mb(ind ==1).*div.BE(ind == 1));
    vwm2b(i,1) = vwm2b(i,1)./(sum(div.BE(ind == 1)));
    ind = ismember(no_div.date_CRSP, dates(i)); 
    ewm2b(i,2) = mean(no_div.mb(ind == 1));
    vwm2b(i,2) = sum(no_div.mb(ind ==1).*no_div.BE(ind==1));
    vwm2b(i,2) = vwm2b(i,2)./(sum(no_div.BE(ind == 1)));
end

ewm2b = log(ewm2b);
vwm2b = log(vwm2b); 
ewm2b(:,3) = 100*(ewm2b(:,1) - ewm2b(:,2));
vwm2b(:,3) = 100*(vwm2b(:,1) - vwm2b(:,2));
prem = vwm2b(:,3); 

%% Loading in BW data

root_dir = [filesep 'Volumes' filesep 'rcecrs02' filesep 'ERPv3'];
stata_dir = [filesep 'Volumes' filesep 'rcecrs02' filesep 'ERPv3' filesep 'Data'];
cd(root_dir)

bw = readtable('Data/Input/Baker Wurgler/bw_data.csv','TreatAsEmpty', {'.'});
bw.ripo = str2double(bw.ripo);
bw = bw(154:end-2,:);

%%% plotting to compare
plot(dates(1:579), [prem(1:579,1) bw.pdnd])
corr([prem(1:579,1) bw.pdnd], 'rows', 'complete')
datetick('x')
title('Dividend Premium Comparison')
ylabel('Dividend Premium')
legend('My Measure', 'BW Measure')
saveas(gcf, 'div_prem_comp_quarterly.png')