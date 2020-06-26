clear 
clc

%%% loading in data
CRSP = readtable("/Volumes/rcecrs02/ERPv3/CRSP.csv");
COMPU = readtable("/Volumes/rcecrs02/ERPv3/CRSP_COMPUA.csv");

%%% cleaning
CRSP.date = datenum(char(string(CRSP.date)), 'yyyymmdd');
COMPU.date = datenum(char(string(COMPU.datadate)), 'yyyymmdd');
CRSP.RET = str2double(CRSP.RET);
CRSP.RETX = str2double(CRSP.RETX);
CRSP.year = year(CRSP.date);
COMPU.year = year(COMPU.date); 
COMPU.PERMNO = COMPU.LPERMNO;
COMPU.PERMCO = COMPU.LPERMCO;
COMPU = sortrows(COMPU, [1 29]);


%% Various Types of CRSP/COMPU Merges

%%% matching by year (naive)
data = innerjoin(CRSP, COMPU, 'keys', {'year', 'PERMNO'});
data = sortrows(data, [1:2]);

%%% last avaliable book value
COMPU.start = nan(height(COMPU),1);
COMPU.dates = datenum(year(COMPU.date), month(COMPU.date), 1);
COMPU.final = nan(height(COMPU),1);

for i = 2:height(COMPU)-1
    if COMPU.GVKEY(i) == COMPU.GVKEY(i + 1)
        COMPU.final(i) = datenum(year(COMPU.date(i+1)), month(COMPU.date(i+1)), 1);
        COMPU.lagdvt 
    end
    if COMPU.GVKEY(i) == COMPU.GVKEY(i-1)
        COMPU.start(i) = datenum(year(COMPU.date(i-1)), month(COMPU.date(i-1)), 1);
    end
end

data = outerjoin(CRSP,COMPU, 'keys', {'PERMNO'});
ind = (data.date_CRSP >= data.start) .* (data.date_CRSP <= data.dates);
data = data(ind == 1, :);

%%% book value in last december
COMPU.dec_date = year(COMPU.date); 

ind = month(CRSP.date) == 12;
CRSP.dec_date(ind == 1) = year(CRSP.date(ind == 1)); 
CRSP.dec_date(ind == 0) = year(CRSP.date(ind == 0)) - 1;

data = innerjoin(CRSP,COMPU, 'keys', {'dec_date', 'PERMNO'});
data = sortrows(data, [1:2]);

%%% book value in last july
ind = month(COMPU.date) <= 7; 
COMPU.july_date(ind) = datenum(year(COMPU.date(ind)), 7.*ones(sum(ind),1),ones(sum(ind),1)); 
COMPU.july_date(ind == 0) = datenum(year(COMPU.date(ind==0))+1,7.*ones(sum(ind ==0),1),ones(sum(ind ==0),1)); 

ind = month(CRSP.date) < 7; 
CRSP.july_date(ind) = datenum(year(CRSP.date(ind))-1,7.*ones(sum(ind),1),ones(sum(ind),1));
CRSP.july_date(ind ==0) =  datenum(year(CRSP.date(ind == 0)), 7.*ones(sum(ind ==0),1),ones(sum(ind ==0),1));
data = innerjoin(CRSP,COMPU, 'keys', {'july_date', 'PERMNO'});
data = sortrows(data, [1:2]);

%% GIVEN MATCH, now compute output

%%% deleting duplicates
data = unique(data, 'stable'); 

%%% dropping data if wrong SHRCD code or Exchange
ind = (data.SHRCD < 12).*(data.SHRCD > 9);
data = data(ind == 1, :);

ind = (data.EXCHCD > 0) .* (data.EXCHCD < 4);
data = data(ind == 1,:);

%%% dropping data in wrong industry
ind = (data.sic > 4899).*(data.sic< 4950);
data = data(ind == 0,:);
ind = (data.sic>5999).*( data.sic < 7000);
data = data(ind == 0,:);

%%% dropping data if missing certain variables
ind = (isnan(data.PRC) == 0).*(isnan(data.SHROUT) == 0)...
    .*(isnan(data.dvpsx_f) == 0); %asset total, prc, shares outstanding, dividend
data = data(ind == 1, :); 

ind = (isnan(data.pstkr)).* (isnan(data.pstkl)).*(isnan(data.upstk));
data = data(ind ==0,:);

ind = (isnan(data.seq)).*(isnan(data.lt)).*(isnan(data.ceq));
data = data(ind == 0,:);

ind = (data.PRC == 0);
data = data(ind == 0,:);

%%% creating variables
data.div_COMPU = data.dvpsx_f;
data.div_CRSP = (data.RET - data.RETX);
data.txditc(isnan(data.txditc)) = 0;

data.SHE = data.seq;
ind = isnan(data.SHE);
data.SHE(ind) = data.ceq(ind) + data.pstk(ind);
ind = isnan(data.SHE);
data.SHE(ind) = data.at(ind) - data.lt(ind) - data.mib(ind);
data.SHE = data.SHE + data.txditc;
data.PS = data.pstkl;
ind = isnan(data.PS);
data.PS(ind) = data.pstkr(ind);
ind = isnan(data.PS);
data.PS(ind) = data.upstk(ind);
ind = isnan(data.PS);
data.PS(ind) = 0;
 
data.BE = (data.SHE-data.PS)* 1000000;
data.BA = (data.at)* 1000000;
data.mktvl = abs(data.PRC).*data.SHROUT*1000;
data.mb= data.mktvl./data.BE; 

%%% dropping data if too small book value or book assets
ind = (data.BE > 250000).*(data.BA > 500000);
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
bw = bw(48:end,:);

%%% plotting to compare
plot(dates(1:687), [prem(1:687,1) bw.pdnd])
corr([prem(1:687,1) bw.pdnd], 'rows', 'complete')
datetick('x')
title('Dividend Premium Comparison')
ylabel('Dividend Premium')
legend('My Measure', 'BW Measure')
saveas(gcf, 'div_prem_comp_naive_match.png')