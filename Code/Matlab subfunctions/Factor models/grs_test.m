function [grs alpha beta stats] = grs_test(returns,factors,method,lags)
%GRS_TEST Gibbons-Ross-Shanken test of zero pricing errors in a factor model
%   GRS = GRS_TEST(RETURNS,FACTORS) gives the grs statistic (a number) for
%   the null hypothesis that that all intercepts (pricing errors) of
%   time-series regression of RETURNS on FACTORS are jointly zero. FACTORS is a
%   T-by-K matrix of T observations of K pricing factors. RETURNS is a T-by-N
%   matrix of T observations of the excess returns of N traded assets. Note: the matrix of
%   factors should not include a column of ones.
%
%   [GRS ALPHA] = GRS_TEST(RETURNS,FACTORS) gives the 1-by-N vector of intercepts
%   (or pricing errors).
%
%   [GRS ALPHA BETA] = GRS_TEST(RETURNS,FACTORS) gives the K-by-N
%   matrix of regression coefficients.
%
%   [GRS ALPHA BETA STATS] = GRS_TEST(RETURNS,FACTORS) gives a structure
%   with the following statistics: R^2 of the regression (stats.R2),
%   adjusted $R^2$ of the regression (stats.adjR2), average absolute
%   intercept (stats.absDev), standard errors of ALPHA (stats.stdErr),
%   t-statisitics for ALPHA (stats.tstat), Sharpe Ratio of the intercepts
%   (stats.SR), p-value of the grs test (stats.p), variance-covariance of
%   FACTORS (stats.sigmaF), variance-covariance the error term
%   (stats.sigmaEpsilon), Sharpe ratio of the ex-post tangency portfolio constructed
%   using FACTORS only (stats.sharpeF), Sharpe ratio of the ex-post tangency portfolio
%   constructed using FACTORS and RETURNS (stats.sharpeAll), 1, 5 and 10%
%   critical values for the grs test stats.criticalValue.
%
%   [...] = GRS_TEST(RETURNS,FACTORS,METHOD) uses the string METHOD to select alternative
%   ways to estimate the model and perform the test. Options are 'GRS'
%   (default), 'GMM' and 'Asymptotic GRS'.
%
%   For more detailed information, see grs_help.m and grs_help.html
%
%   References:
%      [1] Goyal A. (2011) Empirical cross-sectional asset pricing: a
%      survey, Financ Mark Portf Manag.
%
%   Copyright 2011 Fernando M. Duarte
%   $Revision: 0.1 $  $Date: 01/26/2012$

%% Check integrity of inputs
if nargin < 2
    error('grs_test:inputs','Not enough input arguments')
elseif nargin == 2
    T = size(factors);
    if size(returns,1)~=T
        error('grs_test:size','The returns and factors have a different number of observations')
    end
    method = 'GRS'; % pick the default method
elseif nargin == 3
    if ~any(strcmpi(method,{'GRS','GMM','Asymptotic GRS'}))
        error('grs_test:method',['The method ' method ' is not supported'])
    end
    T = size(factors,1);
    lags = floor(T^(1/4)); % default lags for Newey-West
elseif nargin == 4 && strcmpi(method,'GMM')
    if rem(lags,1) ~= 0
        error('grs_test:lags','The number of lags must be an integer')
    end
else
    error('grs_test:inputs','Too many input arguments')
end

%% Eliminate missing data
% check for missing data
missing = any(isnan(returns),2) | any(isnan(factors),2);
infinite = any(isinf(returns),2) | any(isinf(factors),2);
% eliminate rows with missing data
returns(missing | infinite) = [];
factors(missing | infinite) = [];
% total number of observations may have changed
[T K] = size(factors);
N = size(returns,2);

%% Test rank condition
if T-K <= 0
    grs = nan;
    alpha = nan(1,N);
    beta = nan(K,N);
    
    stats.sigmaF = nan(K,K); % variance covariance matrix of factors, normalized by T-K
    stats.sigmaEpsilon = nan(N,N); % variance covariance matrix of residuals, normalized by T
    stats.R2 = nan;
    stats.adjR2 = nan;
    stats.absDev = nan;
    stats.stdErr = nan(1,N); % standard errors of alpha
    stats.tstat = nan(1,N); % tstats for alpha
    stats.SR = nan; % Sharpe ratio of the intercepts
    stats.sharpeF = nan; % Sharpe ratio of the ex-post tangency portfolio using the K factors only
    stats.SharpeAll = nan; % Sharpe ratio of the ex-post tangency portfolio using the K factors and the N assets
    stats.p = nan;
    stats.criticalValue.tenPercent = nan;
    stats.criticalValue.fivePercent = nan;
    stats.criticalValue.onePercent = nan;
    
    return
end

%% Main program

if any(strcmpi(method,{'GRS','Asymptotic GRS'}))
    % collect the variables to run all regressions at once
    y = struct([]);
    x = struct([]);
    for j=1:N
        y(j).eq = returns(:,j);
    end
    for j=1:N
        x(j).eq = [ones(T,1) factors];
    end
    % run sur regression
    %addpath('regress','distrib')
    resultsOls = sur(N,y,x);
    %rmpath('regress','distrib')
    % find regression coefficients
    estimatesOls = [resultsOls.beta];
    alpha = estimatesOls(1,:);
    beta = estimatesOls(2:end,:);
    residuals = [resultsOls.resid];
    % find statistics
    % muHat = mean(returns);
    muHatF = mean(factors);
    stats.sigmaF = cov(factors)*(T-1)/(T-K); % variance covariance matrix of factors, normalized by T-K
    stats.sigmaEpsilon = cov(residuals,1); % variance covariance matrix of residuals, normalized by T
    if  strcmpi(method,'GRS')
        if T-N-K <= 0 || rcond(stats.sigmaEpsilon)<eps
            grs = nan;
            alpha = nan(1,N);
            beta = nan(K,N);
            
            stats.sigmaF = nan(K,K); % variance covariance matrix of factors, normalized by T-K
            stats.sigmaEpsilon = nan(N,N); % variance covariance matrix of residuals, normalized by T
            stats.R2 = nan;
            stats.adjR2 = nan;
            stats.absDev = nan;
            stats.stdErr = nan(1,N); % standard errors of alpha
            stats.tstat = nan(1,N); % tstats for alpha
            stats.SR = nan; % Sharpe ratio of the intercepts
            stats.sharpeF = nan; % Sharpe ratio of the ex-post tangency portfolio using the K factors only
            stats.SharpeAll = nan; % Sharpe ratio of the ex-post tangency portfolio using the K factors and the N assets
            stats.p = nan;
            stats.criticalValue.tenPercent = nan;
            stats.criticalValue.fivePercent = nan;
            stats.criticalValue.onePercent = nan;
            
            return
        end

        grs = ((T-N-K)/N)*(1+muHatF*(stats.sigmaF\muHatF'))^(-1)*(alpha*(stats.sigmaEpsilon\alpha')); % grs statistic
        stats.p = 1-fcdf(grs,N,T-N-K);
        stats.criticalValue.tenPercent = fsolve(@(cc) 1-fcdf(cc,N,T-N-K)-0.1,1,optimset('Display','none'));
        stats.criticalValue.fivePercent = fsolve(@(cc) 1-fcdf(cc,N,T-N-K)-0.05,1,optimset('Display','none'));
        stats.criticalValue.onePercent = fsolve(@(cc) 1-fcdf(cc,N,T-N-K)-0.01,1,optimset('Display','none'));
    else
        % method is 'Asymptotic GRS'
        
        if rcond(stats.sigmaEpsilon)<eps
            grs = nan;
            alpha = nan(1,N);
            beta = nan(K,N);
            
            stats.sigmaF = nan(K,K); % variance covariance matrix of factors, normalized by T-K
            stats.sigmaEpsilon = nan(N,N); % variance covariance matrix of residuals, normalized by T
            stats.R2 = nan;
            stats.adjR2 = nan;
            stats.absDev = nan;
            stats.stdErr = nan(1,N); % standard errors of alpha
            stats.tstat = nan(1,N); % tstats for alpha
            stats.SR = nan; % Sharpe ratio of the intercepts
            stats.sharpeF = nan; % Sharpe ratio of the ex-post tangency portfolio using the K factors only
            stats.SharpeAll = nan; % Sharpe ratio of the ex-post tangency portfolio using the K factors and the N assets
            stats.p = nan;
            stats.criticalValue.tenPercent = nan;
            stats.criticalValue.fivePercent = nan;
            stats.criticalValue.onePercent = nan;
            
            return
        end
        grs = T*(1+muHatF*(stats.sigmaF\muHatF'))^(-1)*(alpha*(stats.sigmaEpsilon\alpha')); % grs statistic
        stats.p = 1-chi2cdf(grs,N);
        stats.criticalValue.tenPercent = fsolve(@(cc) 1-chi2cdf(cc,N)-0.1,N,optimset('Display','none'));
        stats.criticalValue.fivePercent = fsolve(@(cc) 1-chi2cdf(cc,N)-0.05,N,optimset('Display','none'));
        stats.criticalValue.onePercent = fsolve(@(cc) 1-chi2cdf(cc,N)-0.01,N,optimset('Display','none'));
    end
    stats.R2 = mean([resultsOls.rsqr]);
    stats.adjR2 = mean([resultsOls.rbar]);
    stats.absDev = mean(abs(alpha));
    stats.stdErr = ((1+muHatF*(stats.sigmaF\muHatF'))*diag(stats.sigmaEpsilon)/T)'; % standard errors of alpha
    tstat = [resultsOls.tstat];
    stats.tstat = tstat(1,:); % tstats for alpha
    stats.sharpeF = (muHatF*(stats.sigmaF\muHatF'))^(1/2); % Sharpe ratio of the ex-post tangency portfolio using the K factors only
    
    if rcond(stats.sigmaEpsilon) > eps
        stats.SR = (alpha*(stats.sigmaEpsilon\alpha'))^(1/2); % Sharpe ratio of the intercepts
        stats.SharpeAll = (stats.SR^2+stats.sharpeF^2)^(1/2); % Sharpe ratio of the ex-post tangency portfolio using the K factors and the N assets
    else
        stats.SR = nan; % Sharpe ratio of the intercepts
        stats.SharpeAll = nan; % Sharpe ratio of the ex-post tangency portfolio using the K factors and the N assets
        
    end
    
else
    %% use GMM
    % collect the variables to run all regressions at once
    y = reshape(returns,T*N,1);
    x = kron(eye(N),[ones(T,1) factors]);
    % run newey-west regression
    %addpath('regress','distrib')
    resultsOls = nwest(y,x,lags);
    %rmpath('regress','distrib')
    % find regression coefficients
    estimatesOls = reshape(resultsOls.beta,K+1,N);
    alpha = estimatesOls(1,:);
    beta = estimatesOls(2:end,:);
    residuals = reshape(resultsOls.resid,T,N);
    % find statistics
    % muHat = mean(returns);
    muHatF = mean(factors);
    stats.sigmaF = cov(factors)*(T-1)/(T-K); % variance covariance matrix of factors, normalized by T-K
    stats.sharpeF = (muHatF*(stats.sigmaF\muHatF'))^(1/2); % Sharpe ratio of the ex-post tangency portfolio using the K factors only
    stats.sigmaEpsilon = cov(residuals,1); % variance covariance matrix of residuals, normalized by T
    varianceOfCoefficients = resultsOls.V; % variance covariance matrix of coefficients
    stats.stdErr = sqrt(diag(varianceOfCoefficients(1:K+1:end,1:K+1:end)))'; % standard errors of alpha
    if rcond(varianceOfCoefficients(1:K+1:end,1:K+1:end))<eps
        grs = nan;
        alpha = nan(1,N);
        beta = nan(K,N);
        
        stats.sigmaF = nan(K,K); % variance covariance matrix of factors, normalized by T-K
        stats.sigmaEpsilon = nan(N,N); % variance covariance matrix of residuals, normalized by T
        stats.R2 = nan;
        stats.adjR2 = nan;
        stats.absDev = nan;
        stats.stdErr = nan(1,N); % standard errors of alpha
        stats.tstat = nan(1,N); % tstats for alpha
        stats.SR = nan; % Sharpe ratio of the intercepts
        stats.sharpeF = nan; % Sharpe ratio of the ex-post tangency portfolio using the K factors only
        stats.SharpeAll = nan; % Sharpe ratio of the ex-post tangency portfolio using the K factors and the N assets
        stats.p = nan;
        stats.criticalValue.tenPercent = nan;
        stats.criticalValue.fivePercent = nan;
        stats.criticalValue.onePercent = nan;
    else
        grs = T*alpha*(varianceOfCoefficients(1:K+1:end,1:K+1:end)\alpha'); % grs statistic
        stats.p = 1-chi2cdf(grs,N);
        stats.criticalValue.tenPercent = fsolve(@(cc) 1-chi2cdf(cc,N)-0.1,N,optimset('Display','none'));
        stats.criticalValue.fivePercent = fsolve(@(cc) 1-chi2cdf(cc,N)-0.05,N,optimset('Display','none'));
        stats.criticalValue.onePercent = fsolve(@(cc) 1-chi2cdf(cc,N)-0.01,N,optimset('Display','none'));
        stats.R2 = resultsOls.rsqr;
        stats.adjR2 = resultsOls.rbar;
        stats.absDev = mean(abs(alpha));
        stats.tstat = resultsOls.tstat(1:K+1:end)'; % tstats for alpha
        if rcond(stats.sigmaEpsilon) > eps
            stats.SR = (alpha*(stats.sigmaEpsilon\alpha'))^(1/2); % Sharpe ratio of the intercepts
            stats.SharpeAll = (stats.SR^2+stats.sharpeF^2)^(1/2); % Sharpe ratio of the ex-post tangency portfolio using the K factors and the N assets
        else
            stats.SR = nan; % Sharpe ratio of the intercepts
            stats.SharpeAll = nan; % Sharpe ratio of the ex-post tangency portfolio using the K factors and the N assets
        end
    end
end


