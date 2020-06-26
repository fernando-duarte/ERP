%% Three-stage OLS where priced factors and forecasting variables can be
%% restricted.
clear out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X: Forecasting variables
% Xt: Forecasting variables at time t only
% PFacs: Priced Factors
% U: Priced innovations

% CONSTANT WINDOW

% Define Variables for Y (Rets(t+1)), Xt (Forecasters(t)
X = Forecasters.in; PFacs = PricedFacs.in;
Xt = X(1:end-1, :); Zt = [ones(size(Xt,1),1) Xt];
Rx = Rets.in; Y = Rx(2:end,:);
l = size(PFacs,2); [n,k] = size(X);

% 1. OLS: Regress priced factors on lagged forecasting variables and a
% constant; take residuals
[bpriced, U, ~, R2b] = ols(PFacs(2:end,:),Xt);

% 2. Regress returns on Factors and U
[gamma, E, ~, R2g, yhat] = ols(Y, [Xt U]);

gamma1 = gamma(1:size(Zt,2),:)';
beta = gamma(size(Zt,2)+1:end,:)';

% 3. Regress betas on gammas to estimate lambdas
[lambda, ~, ~, lR2] = ols(gamma1,beta,false);

%Store Parameter Estimates
out.U.in = U; out.beta.in = beta; out.lambda.in = lambda;

%Estimate X*lambdas
out.XlamU.in = Zt*lambda' + out.U.in;
out.Xlam.in = [ones(size(X,1),1) X] * lambda';

%HPR
out.HPRest.in = [ones(size(X,1),1) X]*lambda'*beta';

%RMSE, R^2
out.eps.in = out.HPRest.in(1:end-1,:)-Rx(2:end,:);
out.RMSE.ts.in = sqrt(mean(out.eps.in.^2));
out.RMSE.xs.in = sqrt(mean(out.eps.in.^2,2));


TSS = sum((Rx(2:end,:) - repmat(mean(Rx(2:end,:),1),size(Rx(2:end,:),1),1)).^2);
RSS = sum((Rx(2:end,:) - out.HPRest.in(1:end-1,:)).^2);
out.Rsqr.ts.in = 1-RSS./TSS;

MeanRe = mean(Rx(2:end,:));
TSS = sum((MeanRe - mean(MeanRe)).^2);
RSS = sum((MeanRe - mean(out.HPRest.in(1:end-1,:),1)).^2);
out.Rsqr.xs.in = 1-RSS./TSS;

%RECURSIVE (EXPANDING WINDOW)
%Build Structure for results
out.U.rec=cell(size(Forecasters.in,1),1);
out.beta.rec=cell(size(Forecasters.in,1),1);
out.lambda.rec=cell(size(Forecasters.in,1),1);
out.HPRest.rec=zeros(size(Rets.in));

% TRAINING PERIOD

t_initial = date3-date1; %find end of training period

% Define Variables for Y (Rets(t+1)), Xt (Forecasters.in(t)
X = Forecasters.rec{t_initial}; PFacs = PricedFacs.rec{t_initial};
Xt = X(1:end-1, :); Zt = [ones(size(Xt,1),1) Xt];
Rx = Rets.rec{t_initial}; Y = Rx(2:end,:);

% 1. OLS: Regress priced factors on lagged forecasting variables and a
% constant; take residuals
bpriced = olsgmm(PFacs(2:end,:),Zt,lags,-1);
U = PFacs(2:end,:)-Zt*bpriced;

% 2. Regress returns on Factors and U
[gamma] = olsgmm(Y,[Zt U],lags,-1);

gamma1 = gamma(1:size(Zt,2),:)'; a = gamma(1,:)'; c = gamma(2:size(Zt,2),:)';
beta = gamma(size(Zt,2)+1:end,:)'; e_hat = Y - [Zt U]*gamma;

% 3. Regress betas on gammas to estimate lambdas
[lambda] = olsgmm(gamma1,beta,lags,-1);

    %Store results
for T=1:t_initial
    out.beta.rec{T} = beta;
    out.lambda.rec{T} = lambda;
    out.HPRest.rec(T,:) = [1 X(T,:)]*lambda'*beta';
    out.U.rec{T}=U(1:T-1,:); 
    out.Xlam.rec(T,1:l)=[1 X(T,:)]*out.lambda.rec{T}'; 
end
    
%EXPANDING WINDOW
for T = t_initial+1:n
    % Define Variables for Y (Rets(t+1)), Xt (Forecasters(t)
    X = Forecasters.rec{T}; PFacs = PricedFacs.rec{T};
    Xt = X(1:end-1, :); Zt = [ones(size(Xt,1),1) Xt];
    Rx = Rets.rec{T}; Y = Rx(2:end,:);

    % 1. OLS: Regress priced factors on lagged forecasting variables and a
    % constant; take residuals
    bpriced = olsgmm(PFacs(2:end,:),Zt,lags,-1);
    U = PFacs(2:end,:)-Zt*bpriced;
    
    % 2. Regress returns on Factors and U
    [gamma] = olsgmm(Y,[Zt U],lags,-1);

    gamma1 = gamma(1:size(Zt,2),:)'; a = gamma(1,:)'; c = gamma(2:size(Zt,2),:)';
    beta = gamma(size(Zt,2)+1:end,:)'; e_hat = Y - [Zt U]*gamma;

    % 3. Regress betas on gammas to estimate lambdas
    [lambda] = olsgmm(gamma1,beta,lags,-1);


    % Store results
    out.beta.rec{T} = beta;
    out.lambda.rec{T} = lambda;
    out.HPRest.rec(T,:) = [1 X(end,:)]*lambda'*beta';
    out.U.rec{T} = U;
    out.Xlam.rec(T,1:l)=[1 X(end,:)]*out.lambda.rec{T}';
end

out.eps.rec = out.HPRest.rec(1:end-1,:)-Rets.in(2:end,:);
out.RMSE.ts.rec = sqrt(mean(out.eps.rec.^2));
out.RMSE.xs.rec = sqrt(mean(out.eps.rec.^2,2));