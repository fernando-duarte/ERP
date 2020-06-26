function [lambda errors lambda_t errors_t stats] = fama_macbeth(returns,betas,method,shanken,var_method,constant)
%FAMA_MACBETH Fama-MacBeth procedure to estimate risk premia of pricing factors
%   LAMBDA = FAMA_MACBETH(RETURNS,BETAS) gives the 1-by-K vector of factor
%   risk premia. RETURNS is a T-by-N matrix of T observations of the excess
%   returns of N traded assets. BETAS is an N-by-1 cell, where each cell
%   contains a T-by-K matrix of factor exposures.
%
%   [LAMBDA ERRORS]= FAMA_MACBETH(RETURNS,BETAS) gives the 1-by-N vector of
%   pricing errors.
%
%   [LAMBDA ERRORS LAMBDA_T]= FAMA_MACBETH(RETURNS,BETAS) gives
%   a T-by-K vector that contains the time-series of estimates of the risk
%   premia.
%
%   [LAMBDA ERRORS LAMBDA_T ERRORS_T]= FAMA_MACBETH(RETURNS,BETAS) gives
%   a T-by-N vector that contains the time-series of estimates of the
%   pricing errors.
%
%   [LAMBDA ERRORS LAMBDA_T ERRORS_T STATS]= FAMA_MACBETH(RETURNS,BETAS) gives a structure
%   with the following statistics: the K-by-K variance-covariance matrix of the estimates of lambda
%   (stats.varLambda); the N-by-N variance-covariance matrix of the pricing errors
%   (stats.varErrors); the 1-by-K p-value of the null hypothesis that each
%   LAMBDA is equal to zero (stats.pLambda); the 1-by-K vector of t-statistics of LAMBDA
%   (stats.tLambda); the 1-by-N p-value of the null hypothesis that each
%   pricing error is equal to zero (stats.pError); the 1-by-N vector of t-statistics of
%   the pricing errors (stats.tErrors); the p-value of the null hypothesis that
%   all pricing errors are jointly zero (stats.pAllErrors).
%
%   [...] = FAMA_MACBETH(RETURNS,BETAS,METHOD) uses the string METHOD to select alternative
%   ways to estimate the model and statistics. Options are 'Standard' (the
%   default) and 'Constant betas'.
%
%   [...] = FAMA_MACBETH(RETURNS,BETAS,METHOD,SHANKEN) uses the the K-by-K
%   variance-covariance matrix of the factors to correct standard errors
%   for the error-in-variables incurred by estimating betas. Setting SHANKEN to
%   [] gives no adjustment (the default).
%
%   [...] = FAMA_MACBETH(RETURNS,BETAS,METHOD,SHANKEN,VAR_METHOD) computes a
%   generalized least-squares regression. Options are 'gls' (feasible generalized ls),
%   'wls' (feasible weighted ls) and 'ols' (the default). Set SHANKEN to [] to specify
%   CONSTANT without using the Shanken adjustment.
%
%   [...] = FAMA_MACBETH(RETURNS,BETAS,METHOD,SHANKEN,CONSTANT) If CONSTANT
%   is equal to 1, the cross-sectional regressions include a constant. If
%   CONSTANT is 0 (the default), the regressions don't include a constant.
%   Set SHANKEN to [] and GLS to eye(N) to specify CONSTANT without using the Shanken
%   adjustment or generalized least-squares.

%   References:
%      [1] Goyal A. (2011) Empirical cross-sectional asset pricing: a
%      survey, Financ Mark Portf Manag.
%
%   Copyright 2011 Fernando M. Duarte
%   $Revision: 0.1 $  $Date: 02/16/2012$

%% Check integrity of inputs

% enough input arguments
if nargin < 2
    error('fama_macbeth:inputs','Not enough input arguments')
end

[T,N] = size(returns);
if ~iscell(betas)
    error('fama_macbeth:inputs','The variable BETAS must be of type cell')
elseif length(betas)~=N
    error('fama_macbeth:inputs','BETAS and RETURNS have a different number of observations')
end

% right size for betas
for n=1:N
    nBetas = size(betas{n},1);
    if nBetas~=T
        error('fama_macbeth:inputs',['BETAS{ ' num2str(n) '} and RETURNS have a different number of observations'])
    end
end

% other checks for inputs
if nargin == 2
    method = 'Standard'; % pick the default method
    constant = 0; % pick the default of not including a constant in the regression
    shanken = []; % no shanken adjustment
    gls = eye(size(returns,2)); % ols
elseif nargin == 3
    if ~any(strcmpi(method,{'Standard','Constant betas'}))
        error('fama_macbeth:method',['The method ' method ' is not supported'])
    end
    constant = 0; % pick the default of not including a constant in the regression
    shanken = []; % no shanken adjustment
    gls = eye(size(returns,2)); % ols
elseif nargin == 4
    if (size(shanken,1) ~= size(betas{1},2)) || (size(shanken,2) ~= size(betas{1},2))
        error('fama_macbeth:shanken','The variable SHANKEN must be a square matrix with dimension equal to the number of factors')
    end
    if ~isposdef(shanken)
        error('fama_macbeth:shanken','The variable SHANKEN must be a positive definite matrix')
    end
    constant = 0; % pick the default of not including a constant in the regression
    gls = eye(size(returns,2)); % ols
elseif nargin == 5
    if ~strcmpi(var_method,'ols') && ~strcmpi(var_method,'wls') && ~strcmpi(var_method,'gls')
        error('fama_macbeth:var_method','The VAR_METHOD entered is not supported')
    end
    constant = 0; % pick the default of not including a constant in the regression
elseif nargin == 6
    if constant~=0 && constant~=1
        error('fama_macbeth:constant','The variable CONSTANT must be equal to 0 or 1')
    end
else
    error('fama_macbeth:inputs','Too many input arguments')
end

%% Main program
[T,N] = size(returns);
K = size(betas{1},2);

% transform data for gls. this gives ols if the weighting matrix is the identity (the
% default)

vecOne = ones(N,1);

if ~strcmpi(method,'Constant betas')
    
    
    % initialize variables
    lambda_t = nan(T,K);
    errors_t = nan(T,N);
    coeff_t = nan(T,K+1);
    for t=1:T
        % find betas for time t
        beta_t = cell2mat(cellfun(@(x) x(t,:)',betas,'UniformOutput',false))';
        
        if constant == 0
            % find weighting matrix
            wlsCoeff = regress(returns(t,:)',beta_t);
            wlsError = returns(t,:)'-beta_t*wlsCoeff;
            if strcmpi(var_method,'ols')
                gls = eye(N);
            elseif strcmpi(var_method,'wls')
                gls = T*diag(diag(wlsError*wlsError'));
            elseif strcmpi(var_method,'gls')
                gls = T*(wlsError*wlsError');
            end
            if isposdef(gls)
                sqrtGls = chol(gls);
            else
                sqrtGls = chol(diag(diag(gls)));
                warning('fama_macbeth:gls','Estimate of weighting matrix is not positive definite, using WLS instead');
            end
            beta_t = sqrtGls\beta_t;
            weightedReturns = sqrtGls\returns(t,:)';
            % run one cross-sectional regression for each time period
            % (without a constant)
            lambda_t(t,:) = regress(weightedReturns,beta_t);
            errors_t(t,:) = sqrtGls*(weightedReturns-beta_t*lambda_t(t,:)');
        elseif constant == 1
            % run one cross-sectional regression for each time period
            % (with a constant)
            
            % find weighting matrix
            wlsCoeff = regress(returns(t,:)',[vecOne beta_t]);
            wlsError = returns(t,:)'-[vecOne beta_t]*wlsCoeff;
            if strcmpi(var_method,'ols')
                gls = eye(N);
            elseif strcmpi(var_method,'wls')
                gls = T*diag(diag(wlsError*wlsError'));
            elseif strcmpi(var_method,'gls')
                gls = T*(wlsError*wlsError');
            end
            if isposdef(gls)
                sqrtGls = chol(gls);
            else
                sqrtGls = chol(diag(diag(gls)));
                warning('fama_macbeth:gls','Estimate of weighting matrix is not positive definite, using WLS instead');
            end
            
            beta_t = sqrtGls\[vecOne beta_t];
            weightedReturns = sqrtGls\returns(t,:)';
            
            coeff_t(t,:) = regress(weightedReturns,beta_t);
            %lambda0_t(t) = coeff_t(t,1);
            lambda_t(t,:) = coeff_t(t,2:end);
            errors_t(t,:) = sqrtGls*(weightedReturns-beta_t*coeff_t(t,:)');
        end
    end % time loop
    
    % compute time-average of lambda_t and errors_t
    lambda = mean(lambda_t,1);
    errors = mean(errors_t,1);
    
    % compute stats
    if strcmpi(method,'Standard')
        % compute variance-covariance
        stats.varLambda = bsxfun(@minus,lambda_t,lambda)'*bsxfun(@minus,lambda_t,lambda)/(T^2);
        stats.varErrors = bsxfun(@minus,errors_t,errors)'*bsxfun(@minus,errors_t,errors)/(T^2);
        % compute Chi^2 statistic and pvalue for the test that all errors are jointly
        % zero
        allErrorsChi2 = T*errors*(stats.varErrors\errors'); % Chi^2 with N degrees of freedom
        stats.pAllErrors = 1 - chi2cdf(allErrorsChi2,N);
        if ~isempty(shanken)
            % Shanken adjustment for measurement error of betas
            srFactors = lambda*(shanken\lambda'); % square of the Sharpe ratio of factors
            stats.varLambda = (1/T)*((1+srFactors)*(T*stats.varLambda-shanken)+shanken);
            stats.varErrors = stats.varErrors*(1+srFactors);
            % compute Chi^2 statistic and pvalue for the test that all errors are jointly
            % zero
            allErrorsChi2 = T*(1+srFactors)*errors*(stats.varErrors\errors'); % Chi^2 with N-K degrees of freedom
            stats.pAllErrors = 1 - chi2cdf(allErrorsChi2,N-K);
        end
        % compute tstats and p-values
        stats.tLambda = lambda./sqrt(diag(stats.varLambda/T))';
        [~,stats.pLambda] = ztest(stats.tLambda,0,1,[],[],1);
        stats.tErrors = errors./sqrt(diag(stats.varErrors/T))';
        [~,stats.pErrors] = ztest(stats.tErrors,0,1,[],[],1);
    end
    
elseif strcmpi(method,'Constant betas')
    % compute mean returns and mean betas
    meanRet = (mean(returns,1)');
    meanBetas = cell2mat(cellfun(@mean,betas,'UniformOutput',false)');
    
    % run regression
    if constant == 0
        % find the right weighting matrix
        wlsCoeff = regress(meanRet,meanBetas);
        wlsError = meanRet-meanBetas*wlsCoeff;
        if strcmpi(var_method,'ols')
            gls = eye(N);
        elseif strcmpi(var_method,'wls')
            gls = T*diag(diag(wlsError*wlsError'));
        elseif strcmpi(var_method,'gls')
            gls = T*(wlsError*wlsError');
        end
        if isposdef(gls)
            invGls = inv_posdef(gls);
            sqrtGls = chol(gls);
        else
            invGls = inv_posdef(diag(diag(gls)));
            sqrtGls = chol(diag(diag(gls)));
            warning('fama_macbeth:gls','Estimate of weighting matrix is not positive definite, using WLS instead');
        end
        sandwich = inv_posdef(meanBetas'*invGls*meanBetas); % frequently used matrix
        sigErrors = T*(wlsError*wlsError'); % estimate of T*E[error*error']
        if ~isposdef(sigErrors)
            warning('fama_macbeth:errors','Estimate of the error variance-covariance matrix is not positive definite');
        end
        lambda = sandwich*meanBetas'*invGls*meanRet;
        errors = meanRet - meanBetas*lambda;
        stats.varLambda = (1/T)*sandwich*(meanBetas'*sigErrors*meanBetas)*sandwich;
        proyect = eye(N)-(sqrtGls\meanBetas)*inv_posdef(meanBetas'*invGls*meanBetas)*(sqrtGls\meanBetas)'; % projection matrix
    elseif constant == 1
        meanBetasConst = [ones(N,1) meanBetas]; % include the constant
        wlsCoeff = regress(meanRet,meanBetasConst);
        wlsError = meanRet-meanBetasConst*wlsCoeff;
        if strcmpi(var_method,'ols')
            gls = eye(N);
        elseif strcmpi(var_method,'wls')
            gls = T*diag(diag(wlsError*wlsError'));
        elseif strcmpi(var_method,'gls')
            gls = T*(wlsError*wlsError');
        end
        if isposdef(gls)
            invGls = inv_posdef(gls);
            sqrtGls = chol(gls);
        else
            invGls = inv_posdef(diag(diag(gls)));
            sqrtGls = chol(diag(diag(gls)));
            warning('fama_macbeth:gls','Estimate of weighting matrix is not positive definite, using WLS instead');
        end
        sandwich = inv_posdef(meanBetasConst'*invGls*meanBetasConst); % frequently used matrix
        sigErrors = T*(wlsError*wlsError'); % estimate of T*E[error*error']
        if ~isposdef(sigErrors)
            warning('fama_macbeth:errors','Estimate of the error variance-covariance matrix is not positive definite');
        end
        coeff = sandwich*meanBetasConst'*invGls*meanRet;
        lambda = coeff(2:end);
        errors = meanRet - meanBetasConst*coeff;
        stats.varLambda = (1/T)*sandwich*(meanBetasConst'*sigErrors*meanBetasConst)*sandwich;
        stats.varLambda = stats.varLambda(2:end,2:end); % don't include cov of the constant term
        proyect = eye(N)-(sqrtGls\meanBetasConst)*inv_posdef(meanBetasConst'*invGls*meanBetasConst)*(sqrtGls\meanBetasConst)'; % projection matrix
    end
    % since returns and betas are constant, there's no time series
    % for lambda and errors
    lambda_t = nan;
    errors_t = nan;
    
    % compute variance of errors
    stats.varErrors = (1/T)*proyect*sigErrors*proyect';
    % test whether all errors are jointly zero
    allErrorsChi2 = T*errors'*(stats.varErrors\errors); % Chi^2 with N degrees of freedom
    stats.pAllErrors = 1 - chi2cdf(allErrorsChi2,N);
    
    if ~isempty(shanken)
        % Shanken adjustment for measurement error of betas
        srFactors = lambda'*(shanken\lambda); % square of the Sharpe ratio of factors
        stats.varLambda = (1/T)*((1+srFactors)*(T*stats.varLambda-shanken)+shanken);
        stats.varErrors = stats.varErrors*(1+srFactors);
        % compute Chi^2 statistic and pvalue for the test that all errors are jointly
        % zero
        allErrorsChi2 = (1+srFactors)*allErrorsChi2; % Chi^2 with N-K degrees of freedom
        stats.pAllErrors = 1 - chi2cdf(allErrorsChi2,N-K);
    end
    % compute tstats and p-values
    stats.tLambda =  lambda./sqrt(diag(stats.varLambda/T));
    [~,stats.pLambda] = ztest(stats.tLambda,0,1,[],[],2);
    stats.tErrors =  errors./sqrt(diag(stats.varErrors/T));
    [~,stats.pErrors] = ztest(stats.tErrors,0,1,[],[],2);
    
end % end method == 'Constant betas'

end % function

