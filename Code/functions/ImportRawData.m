
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import data
[Returns.in, Returns.names] = xlsread(DataFile,RetSheet);
[FFFacs.in, FFFacs.names] = xlsread(DataFile,FFSheet);
[Factors.in, Factors.names] = xlsread(DataFile,FacSheet);

Returns.names = Returns.names(1,2:end);
FFFacs.names = FFFacs.names(1,2:end);
Factors.names = Factors.names(1,2:end);

if freq==12||freq==4 %For mly or qly data
    dater = num2str(Returns.in(:,1)); Returns.in(:,1) = [];
    dater = [str2num(dater(:,1:4)) str2num(dater(:,5:6))];
    
    datef = num2str(FFFacs.in(:,1)); FFFacs.in(:,1) = [];
    datef = [str2num(datef(:,1:4)) str2num(datef(:,5:6))];

    date = num2str(Factors.in(:,1)); Factors.in(:,1) = [];
    date = [str2num(date(:,1:4)) str2num(date(:,5:6))];
end

if freq==52 %For wly data
    dater = num2str(Returns.in(:,1)); Returns.in(:,1) = [];
    dater = [str2num(dater(:,1:4)) str2num(dater(:,5:6)) str2num(dater(:,7:8))];
    
    datef = num2str(FFFacs.in(:,1)); FFFacs.in(:,1) = [];
    datef = [str2num(datef(:,1:4)) str2num(datef(:,5:6)) str2num(datef(:,7:8))];

    date = num2str(Factors.in(:,1)); Factors.in(:,1) = [];
    date = [str2num(date(:,1:4)) str2num(date(:,5:6)) str2num(date(:,7:8))];
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find start and end dates

Startdate = datevec(startdate); Enddate = datevec(enddate); Trainingdate = datevec(trainingdate);
Y1 = Startdate(1); M1 = Startdate(2); Y2 = Enddate(1); M2 = Enddate(2); Y3 = Trainingdate(1); M3 = Trainingdate(2);

date1r = find(dater(:,1)==Y1 & dater(:,2)==M1,1,'first');
if isempty(date1r); display('Return start date not found'); return; end
date1f = find(datef(:,1)==Y1 & datef(:,2)==M1,1,'first');
if isempty(date1f); display('FF Factor start date not found'); return; end
date1 = find(date(:,1)==Y1 & date(:,2)==M1,1,'first');
if isempty(date1); display('Factor start date not found'); return; end

date2r = find(dater(:,1)==Y2 & dater(:,2)==M2,1,'last');
if isempty(date2r); display('Return end date not found'); return; end
date2f = find(datef(:,1)==Y2 & datef(:,2)==M2,1,'last');
if isempty(date2f); display('FF Factor end date not found'); return; end
date2 = find(date(:,1)==Y2 & date(:,2)==M2,1,'last');
if isempty(date2); display('Factor end date not found'); return; end

date3r = find(dater(:,1)==Y3 & dater(:,2)==M3,1,'last');
if isempty(date3r); display('Return training date not found'); return; end
date3f = find(datef(:,1)==Y3 & datef(:,2)==M3,1,'last');
if isempty(date3f); display('FF Factor training date not found'); return; end
date3 = find(date(:,1)==Y3 & date(:,2)==M3,1,'last');
if isempty(date3); display('Factor training date not found'); return; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parse data
FactorsRaw.in = Factors.in(date1:end,:); %store factors from start date on for forecasting
FactorsRaw.dates = date(date1:end,:);
FactorsRaw.names = Factors.names;

Factors.in = Factors.in(date1:date2,:); %store data for 
FFFacs.in = FFFacs.in(date1f:date2f,:);
Returns.in = Returns.in(date1r:date2r,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Transform Returns to Excess Returns; add RMkt
Returns.in = [FFFacs.in(:,strcmpi('Mkt-RF',FFFacs.names)), ...
    Returns.in - repmat(FFFacs.in(:,strcmpi('RF',FFFacs.names)),1,size(Returns.in,2))];
Retuns.names = ['RMkt' Returns.names];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add Equal-Weighted Return Series to Returns (optional)
if EW == 1
    Returns.in = [mean(Returns.in,2) Returns.in];
    Returns.names = ['EW Ret' Returns.names];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Store Factors and Returns in Cells (for OOS forecasting)
%Factors

Factors.rec = cell(size(Factors.in,1),1);
Returns.rec = cell(size(Returns.in,1),1);
for t=1:date3-date1
    Returns.rec{t} = Returns.in(1:date3-date1,:);
    Factors.rec{t} = Factors.in(1:date3-date1,:);
end
for t = date3-date1+1:date2-date1+1
    Returns.rec{t} = Returns.in(1:t,:);
    Factors.rec{t} = Factors.in(1:t,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create a date vector, date_all vector and a tforecast vector for charts
date = datenum([date(date1:date2,:) ones(date2-date1+1,1)]);
date_all = datenum([FactorsRaw.dates ones(size(FactorsRaw.dates,1),1)]);
tforecast = [date(2:end); date(end)+ (date(end)-date(end-1))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display Factor Names and Return Names
display(['Raw Factors: ' Factors.names]);
display(['Raw Returns: ' Returns.names]);


