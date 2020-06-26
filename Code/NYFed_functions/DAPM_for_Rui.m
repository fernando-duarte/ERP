%% Instructions
% Step one: Update Tend to most recent date for which there is data
% Step three: look for %!% - this marks things that need to be updated

clear all; close all;
pth = pwd;
ind = strfind(pth,'Code');
datapath = [pth(1,1:ind-1),'Data'];
resultspath = [pth(1,1:ind-1),'Results'];
addpath([pwd,'/Code'])

datafile = [datapath '\AllData_150508.xlsx'];

% Sample dates
Tstart_sample = 196401;
Tstart_is = 196401;
Tend = 201504;  %!% Change this to latest month there is factor data for
currentdate = '20140826'; %!% date of last run of DAPM_for_QS; need to load DAPM_params_"currentdate"
IndUseKenFrench = false;

[rawfacs,facstr] = xlsread(datafile,'MonthlyFactors');

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
X1list = char('MKT','SMB'); % Use MKT and SMB as cross-sectional pricing factors   
X2list = char('TSY10')      % Use TSY10 as both cross-sectional pricing and predictive factors
X3list = char('TERM','dy2') % Use TERM and dy2 as predictive factors

X1full = facs(:,cell2mat(strindx_vec(facstr,X1list,'exact'))); 
X2full = facs(:,cell2mat(strindx_vec(facstr,X2list,'exact'))); 
X3full = facs(:,cell2mat(strindx_vec(facstr,X3list,'exact'))); 

faclist = char(X1list,X2list,X3list);

% State all returns and factors in basis points per month
X3full(:,1) = X3full(:,1)/100;
X1full = X1full/100;
X2full = X2full/100;

%% In-sample analysis
    
if ~isempty(X1full), X1 = X1full(T1_is:T2,:); else X1 = []; end;        
if ~isempty(X2full), X2 = X2full(T1_is:T2,:); else X2 = []; end;        
if ~isempty(X3full), X3 = X3full(T1_is:T2,:); else X3 = []; end;
if ~isempty(find(isnan([X2 X3]),1)), disp('NaN in factors!!!'); return; end;

[T,k1] = size(X1);
k2 = size(X2,2);
k3 = size(X3,2);
kF = k2+k3;
k = k1+k2+k3;

load(['DAPM_params_', currentdate, '.mat']);

horsel = 1:240;
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
    for h = 1:maxh,
        if h > 2,
            mu_fcsts(:,h) = [eye(k) + sum(Phi_fcsts(:,:,1:h-1),3)]*mu;
        end;
        Phi_fcsts(:,:,h) = Phi^h;
        X_fcsts(:,:,h) = Z*[mu_fcsts(:,h) Phi_fcsts(:,:,h)]';
        F_tilde_fcsts(:,:,h) = [ones(T,1) X_fcsts(:,k1+1:end,h)];    
    end;
else
    X_QS = [X2(2:end,:) X3(2:end,:)];
    Z_QS = [ones(size(X2,1)-1,1) X2(1:end-1,:) X3(1:end-1,:)];
    Psi_QS = ((Z_QS'*Z_QS)\Z_QS'*X_QS)'
    mu = Psi_QS(:,1);
    Phi = Psi_QS(:,2:end);
    %
    mu_fcsts = zeros(kF,maxh);
    mu_fcsts(:,1) = mu;
    Phi_fcsts = zeros(kF,kF,maxh);
    X_fcsts = zeros(T,kF,maxh);
    X = [X2 X3];
    Z = [ones(T,1) X];
    for h = 1:maxh,
        if h > 2,
            mu_fcsts(:,h) = [eye(kF) + sum(Phi_fcsts(:,:,1:h-1),3)]*mu;
        end;
        Phi_fcsts(:,:,h) = Phi^h;
        X_fcsts(:,:,h) = Z*[mu_fcsts(:,h) Phi_fcsts(:,:,h)]';
        F_tilde_fcsts(:,:,h) = [ones(T,1) X_fcsts(:,:,h)];    
    end;   
end

for hor = 1:length(horsel),
    h = horsel(hor);
    Premth_MKT{hor} = ones(T,1);    
    if (hor==1)
        Premth_MKT_cont = ones(T,4);
        for i = 1:4
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Produce Excel Spreadsheet for QS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
horsel_QS = horsel;
horselstr = cell(size(horsel'));
for i = horsel
    horselstr{i} = [num2str(i) '-month ERP'];
end
Premth_MKT_QS = cell2mat(Premth_MKT_ann);

xlswrite('FRBNY_ERP_Rui.xls',Premth_MKT_QS(:,horsel_QS),'FRBNY_ERP','B2');
xlswrite('FRBNY_ERP_Rui.xls',dates(T1_is:T2),'FRBNY_ERP','A2');
xlswrite('FRBNY_ERP_Rui.xls',{'Dates/Horizon'},'FRBNY_ERP','A1:A1');
xlswrite('FRBNY_ERP_Rui.xls',horselstr','FRBNY_ERP','B1');

