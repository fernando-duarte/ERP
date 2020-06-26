%% Instructions
% Step one: Update Tend to most recent date for which there is data
% Step three: look for %!% - this marks things that need to be updated

clear all; close all;
pth = pwd;
ind = strfind(pth,'Code');
datapath = [pth(1,1:ind-1),'Data'];
resultspath = [pth(1,1:ind-1),'Results'];
addpath([pwd,'/Code'])

datafile = [datapath '\AllData_180319.xlsx'];

% Sample dates
Tstart_sample = 196401;
Tstart_is = 196401;
Tend = 201801;  %!% Change this to latest month there is forecasting variables data for (i.e., X2 and X3 variables)
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
X1list = char('MKT','SMB');  % Use MKT and SMB as cross-sectional pricing factors   
X2list = char('TSY10');      % Use TSY10 as both cross-sectional pricing and predictive factors
X3list = char('TERM','dy2'); % Use TERM and dy2 as predictive factors

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
        Premth_MKT_cont = ones(T,1+kF);
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Produce Excel Spreadsheet for QS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
horsel_QS = [1 12 60 120];
horselstr = {'1 month ERP'; '1 year ERP'; '5 year ERP'; '10 year ERP'};
Premth_MKT_QS = cell2mat(Premth_MKT_ann);

xlswrite('FRBNY_ERP_QS.xls',Premth_MKT_QS(:,horsel_QS),'FRBNY_ERP','B2');
xlswrite('FRBNY_ERP_QS.xls',dates(T1_is:T2),'FRBNY_ERP','A2');
xlswrite('FRBNY_ERP_QS.xls',{'Dates/Horizon'},'FRBNY_ERP','A1:A1');
xlswrite('FRBNY_ERP_QS.xls',horselstr','FRBNY_ERP','B1');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Produce ERP Surface Plot for QS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
i = 1;
figure(i)
surf(plotdates(1:end),(1:120)./12,Premth_MKT_QS(1:end,:)','EdgeColor','none');
set(gca,'YTick',[2:2:10],'xdir','reverse');
set(gca,'YTick',[2:2:100],'XTick',floor(plotdates(12:120:end)));
colormap jet
ttt = colorbar('peer',gca,'location','North');
%colorbar;
shading 'interp';
brighten(.5);
grid on
set(gca,'GridLineStyle',':','Xcolor',[.0 .0 .0],'Ycolor',[.0 .0 .0],'Zcolor',[.0 .0 .0]);
axis tight
ylabel('','Fontsize',9);
zlabel('Expected Returns','Fontsize',9);
ylabel('Horizon (years)','Fontsize',9); 
xlabel('Date','Fontsize',9);
title('ACM Equity Risk Premium Surface','Fontsize',10,'Fontname','Arial'); ...
file='ERP_surface.png';
print(gcf,'-dpng','-r300',file);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Look at contributions of individual series to ERP_t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
decomp = 1;

if decomp == 1
    for n =  [1 12 60 120]
        i = i+1;
        figure(i)
        y = Premth_MKT_ann{n};
        [b_hat] = Z\y;
        
        decomp1 = b_hat(1+1)*X(:,1);
        decomp2 = b_hat(1+2)*X(:,2);
        decomp3 = b_hat(1+3)*X(:,3)+b_hat(1);
        decomp  = [decomp1,decomp2,decomp3];
        
       
        plot(plotdates,y,'LineWidth',2);
        hold on
        plot(plotdates,decomp,'LineWidth',1);
        hold off
        legend({'ERP',X2list,'TERM','DY + Constant'},'Location','Southwest',...
            'FontSize',12);
        title([num2str(n),'-Month ERP Decomposition'],'FontSize',16)
        set(figure(i),'PaperOrientation','portrait','PaperType','usletter',...
            'InvertHardCopy','off','PaperPositionMode','auto','Position',...
            [217  258 1200 950],'Color',[1 1 1]);
        print(figure(i),'-dpng',['ERP_decomp_',num2str(n),'_month'])
    end
end
close all;

return

    for n =  [1 12 60 120]
        i = i+1;
        figure(i)
        y = Premth_MKT_ann{n};
        [b_hat] = Z\y;
        
        decomp1 = b_hat(1+1)*X(:,1);
        decomp2 = b_hat(1+2)*X(:,2);
        decomp3 = b_hat(1+3)*X(:,3)+b_hat(1);
        decomp  = [decomp1,decomp2,decomp3];
        
        start = find(plotdates==1990);
       
        plot(plotdates(start:end),y(start:end),'LineWidth',2);
        hold on
        plot(plotdates(start:end),decomp(start:end,:),'LineWidth',1);
        hold off
        legend({'ERP',X2list,'TERM','DY + Constant'},'Location','Southeast',...
            'FontSize',12);
        title([num2str(n),'-Month ERP Decomposition'],'FontSize',16)
        set(figure(i),'PaperOrientation','portrait','PaperType','usletter',...
            'InvertHardCopy','off','PaperPositionMode','auto','Position',...
            [217  258 1200 950],'Color',[1 1 1]);
        print(figure(i),'-dpng',['ERP_decomp_',num2str(n),'_month_short'])
    end
    close all;
    
    for n =  [1 12 60 120]
        i = i+1;
        figure(i)
        y = Premth_MKT_ann{n};
        [b_hat] = Z\y;
        
        decomp1 = b_hat(1+1)*X(:,1);
        decomp2 = b_hat(1+2)*X(:,2);
        decomp3 = b_hat(1+3)*X(:,3)+b_hat(1);
        decomp  = [decomp1,decomp2,decomp3];
        
        start = find(plotdates==1990);
       
        plot(plotdates(start:end),y(start:end),'LineWidth',2);
        hold on
        plot(plotdates(start:end),decomp(start:end,:),'LineWidth',1);
        hold off
        legend({'ERP',X2list,'TERM','DY + Constant'},'Location','Southeast',...
            'FontSize',12);
        title([num2str(n),'-Month ERP Decomposition'],'FontSize',16)
        set(figure(i),'PaperOrientation','portrait','PaperType','usletter',...
            'InvertHardCopy','off','PaperPositionMode','auto','Position',...
            [217  258 1200 950],'Color',[1 1 1]);
        print(figure(i),'-dpng',['ERP_decomp_',num2str(n),'_month_short'])
    end
    close all;    
    
    for n =  [1 12 60 120]
        i = i+1;
        figure(i)
        y = Premth_MKT_ann{n};
        [b_hat] = Z\y;
        
        decomp1 = b_hat(1+1)*X(:,1);
        decomp2 = b_hat(1+2)*X(:,2);
        decomp3 = b_hat(1+3)*X(:,3)+b_hat(1);
        decomp  = [decomp1,decomp2,decomp3];
        
        start = find(plotdates==2000);
       
        plot(plotdates(start:end),y(start:end),'LineWidth',2);
        hold on
        plot(plotdates(start:end),decomp(start:end,:),'LineWidth',1);
        hold off
        legend({'ERP',X2list,'TERM','DY + Constant'},'Location','Southeast',...
            'FontSize',12);
        title([num2str(n),'-Month ERP Decomposition'],'FontSize',16)
        set(figure(i),'PaperOrientation','portrait','PaperType','usletter',...
            'InvertHardCopy','off','PaperPositionMode','auto','Position',...
            [217  258 1200 950],'Color',[1 1 1]);
        print(figure(i),'-dpng',['ERP_decomp_',num2str(n),'_month_shorter'])
    end
    close all;
    
    for n =  [1 12 60 120]
        i = i+1;
        figure(i)
        y = Premth_MKT_ann{n};
        [b_hat] = Z\y;
        
        decomp1 = b_hat(1+1)*X(:,1);
        decomp2 = b_hat(1+2)*X(:,2);
        decomp3 = b_hat(1+3)*X(:,3)+b_hat(1);
        decomp  = [decomp1,decomp2,decomp3];
        
        start = find(plotdates==2000);
       
        plot(plotdates(start:end),y(start:end),'LineWidth',2);
        hold on
        plot(plotdates(start:end),y(start:end)+100*Z(start:end,2),'LineWidth',2);
        hold off
        legend({'ERP', 'Expected Stock Returns'},'Location','Southeast',...
            'FontSize',12);
        title([num2str(n),'-Month ERP and Expected Stock Returns'],'FontSize',16)
        set(figure(i),'PaperOrientation','portrait','PaperType','usletter',...
            'InvertHardCopy','off','PaperPositionMode','auto','Position',...
            [217  258 1200 950],'Color',[1 1 1]);
        print(figure(i),'-dpng',['ERP_ExpStockRet_',num2str(n),'_month_shorter'])
    end
    close all;    
    
    
    start = find(plotdates==2010);
    Z_recent = Z(start:end,2:end);
    plot(100*Z_recent(:,1))
    plot(100*Z_recent(:,2))
    plot(Z_recent(:,3))
    
    
    
    
    
