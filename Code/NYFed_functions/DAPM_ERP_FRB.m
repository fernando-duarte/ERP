
clear all; close all;
pth = pwd;
ind = strfind(pth,'Code');
datapath = [pth(1,1:ind-1),'Data'];
resultspath = [pth(1,1:ind-1),'Results'];

datafile = [datapath '\eprem_series.xlsx'];

% Sample dates
Tstart_sample = 198501;
Tstart_is = 198501;
Tend = 201511;  %!% Change this to latest month there is data for

[rawfacs,facstr] = xlsread(datafile,'facs');

facstr = facstr(1,2:end)';
dates = rawfacs(:,1);
facs = rawfacs(:,2:end);

T0 = find(dates == Tstart_sample);
T1_is = find(dates == Tstart_is);
T2 = find(dates == Tend);
plotdates = floor(Tstart_is/100)+(Tstart_is-100*floor(Tstart_is/100))/12:1/12:floor(Tend/100)+(Tend-100*floor(Tend/100))/12';
datesnum1 = datenum([floor(Tstart_sample/100) (Tstart_sample-100*floor(Tstart_sample/100)) 1 0 0 0]);
datesnum2 = datenum([floor(Tend/100) (Tend-100*floor(Tend/100)) 1 0 0 0]);

%% Prepare data
retlist = char('10size','zcbret');
Rfull = []; allrets = [];
for i = 1:size(retlist,1),
    [rets{i},retstr{i}] = xlsread(datafile,deblank(retlist(i,:))); %
    dates = rets{i}(:,1);
    dates_num = datenum([floor(dates/100) dates-100*floor(dates/100) ones(length(dates),1) zeros(length(dates),3)]);
    [~,t1] = intersect(dates_num,datesnum1);
    [~,t2] = intersect(dates_num,datesnum2);
    Rfull = [Rfull rets{i}(t1:t2,2:end)];
    allrets = char(allrets,char(retstr{i}(1,2:end)'));
end;
%RFfull = xlsread(datafile,'Rf'); 
%RFfull = RFfull(:,2); 
allrets = allrets(2:end,:);

X1list = char('MKT','SMB'); % Use MKT and SMB as cross-sectional pricing factors   
X2list = char('TSY10');      % Use TSY10 as both cross-sectional pricing and predictive factors
X3list = char('TERM','LN_DY','LN_E/D'); % Use TERM and dy2 as predictive factors
%X3list = char('TERM','LN_DY','LN_P/E','LN_E/D'); % Use TERM and dy2 as predictive factors

X1full = facs(:,cell2mat(strindx_vec(facstr,X1list,'exact'))); 
X2full = facs(:,cell2mat(strindx_vec(facstr,X2list,'exact'))); 
X3full = facs(:,cell2mat(strindx_vec(facstr,X3list,'exact'))); 

faclist = char(X1list,X2list,X3list);

% State all returns and factors in basis points per month
X3full(:,1) = X3full(:,1)/100;
X1full = X1full/100;
X2full = X2full/100;

%% In-sample analysis
%Rx = Rfull/100 - repmat(RFfull(T1_is:T2,:)/100,1,size(Rfull,2)); %  - repmat(RFfull(T1_is:T2,:),1,size(Rfull,2)) - repmat(X4full(T1_is:T2,:),1,size(Rfull,2)
Rx = Rfull/100;
nanind = []; 
for n = 1:size(Rx,2),
    if ~isempty(find(isnan(Rx(:,n)),1)),
        nanind = [nanind n];
    end;
end;
Rx(:,nanind) = [];    
    
if ~isempty(X1full), X1 = X1full(T1_is:T2,:); else X1 = []; end;        
if ~isempty(X2full), X2 = X2full(T1_is:T2,:); else X2 = []; end;        
if ~isempty(X3full), X3 = X3full(T1_is:T2,:); else X3 = []; end;
if ~isempty(find(isnan([X1 X2 X3]),1)), disp('NaN in factors!!!'); return; end;

[T,k1] = size(X1);
k2 = size(X2,2);
k3 = size(X3,2);
kF = k2+k3;
k = k1+k2+k3;

N = size(Rx,2);
Re = Rx(2:end,:);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outputs = DAPM_noTVB_QS(Rx,X1,X2,X3,[],0); % run DAPM with constant betas and VAR parameters

%currentdate = datestr(now, 'YYYYmmDD');
%save(['DAPM_params_', currentdate, '.mat'], 'outputs');


% varlist = char(fieldnames(outputs));
% for i = 1:size(varlist,1),
%     if length(deblank(varlist(i,:)))>2,
%         eval([deblank(varlist(i,:)),' = outputs.',deblank(varlist(i,:)),';']);
%     end;
% end;

%% Calculate Risk Premia

horsel = 1:120;
maxh = max(horsel);
F_tilde_fcsts = zeros(T,1+kF,maxh);
mu = outputs.Psi_CB(1,:)'; Phi = outputs.Psi_CB(2:end,:)';
mu_fcsts = zeros(k,maxh);
mu_fcsts(:,1) = mu;
Phi_fcsts = zeros(k,k,maxh);
X_fcsts = zeros(T,k,maxh);
X = [X1 X2 X3];
Z = [ones(T,1) X];
for h = 1:maxh,
    if h > 2,
        mu_fcsts(:,h) = [eye(k) + sum(Phi_fcsts(:,:,1:h-1),3)]*mu;
    end;
    Phi_fcsts(:,:,h) = Phi^h;
    X_fcsts(:,:,h) = Z*[mu_fcsts(:,h) Phi_fcsts(:,:,h)]';
    F_tilde_fcsts(:,:,h) = [ones(T,1) X_fcsts(:,k1+1:end,h)];    
end;    

for hor = 1:length(horsel),
    h = horsel(hor);
    Premth_MKT{hor} = ones(T,1);    
    if (hor==1)
        Premth_MKT_cont = ones(T,(1+kF));
        for i = 1:(1+kF)
            Lamt_h_i = F_tilde_fcsts(:,i,h)*outputs.Lam_CB(1,i);
            Premth_MKT_cont(:,i)= Premth_MKT_cont(:,i).*(ones(T,1) + Lamt_h_i);
        end
    end
    for hh = 1:h,
        Lamt_h = F_tilde_fcsts(:,:,hh)*outputs.Lam_CB';
        Premth_MKT{h} = Premth_MKT{h}.*(ones(T,1) + Lamt_h(:,1));
    end
    Premth_MKT_ann{hor} = 100*(Premth_MKT{hor}.^(12/h)-1); % produce annual rates of return
end

% Plot ERP
figure(1)
1200*[Premth_MKT{1}-1 Premth_MKT_cont-1]; 
plot(plotdates,1200*[Premth_MKT{1}-1])

% look at contributions of individual series to ERP_t,1
figure(2)
plot(plotdates,Premth_MKT_cont(:,2:end));
legend('TSY10','TERM','DY','ED');