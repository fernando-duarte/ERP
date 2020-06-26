function [hAlpha hBeta hLongRunAlpha hLongRunBeta] = optimal_bandwidth(Y,X,kernel,polyDegree,span)
% OPTIMAL_BANDWIDTH estimates the optimal bandwidth for nonparametric rolling least-squares
% estimators
%
%   HALPHA = OPTIMAL_BANDWIDTH(Y,X) gives the optimal bandwidth for the nonparametric
%	rolling least-squares estimators of regression intercepts.
%   X is a T-by-J matrix of T observations of
%   J independent variables. Y is a T-by-M matrix of T observations
%   of M dependent variables. HALPHA is a 1-by-M matrix of
%   bandwidths. Note: X should not include a column of ones.
%
%   [HALPHA HBETA] = OPTIMAL_BANDWIDTH(Y,X) also returns the a 1-by-M
%   matrix of optimal bandwidth for the slope coefficients.
%
%   [HALPHA HBETA HLONGRUNALPHA] = OPTIMAL_BANDWIDTH(Y,X) also returns the a 1-by-M
%   matrix of optimal bandwidth for the long run (time-integrated) intercepts
%
%   [HALPHA HBETA HLONGRUNALPHA HLONGRUNBETA] = OPTIMAL_BANDWIDTH(Y,X) also returns the a 1-by-M
%   matrix of optimal bandwidth for the long run (time-integrated) slope
%   coefficients
%
%   [...] = OPTIMAL_BANDWIDTH(Y,X,kernel) specifies the kernel. If not
%   specified, the default is 'gaussian'. Other options are
%   'one-sided-uniform','two-sided-uniform',
%   'one-sided exponential' or a used-supplied function handle (the kernel must integrate to one).
%
%   [...] = OPTIMAL_BANDWIDTH(Y,X,kernel,polyDegree) specifies the degree
%   of the polynomial used for the first-stage parametric estimation of the
%   least-squares estimators. If not specified, the default is 6.
%
%   [...] = OPTIMAL_BANDWIDTH(Y,X,kernel,polyDegree,span) specifies the
%   required span of the kernel. If span is non-empty, OPTIMAL_BANDWIDTH
%   just converts span into a bandwidth and returns the same value h for
%   HALPHA, HBETA, HLONGRUNALPHA and HLONGRUNBETA based on:
%       exponential kernel: h = span/(T*log(2)) -- so that span is a halflife
%       uniform kernel: h = 2*span/T -- so that span is half the size of the
%                                   window
%       gaussian kernel: h = span/size(Y,1) -- so that span is a std dev
%       otherwise: h = 1
%
%   References:
%      [1] Ang, A. and Kristensen, D. (2011) Testing Conditional Factor
%      Models
%
%   Copyright 2011 Fernando M. Duarte
%   Revision: 0.1   Date: 02/01/2012

%% Preliminary checks

if nargin < 2
    error('rolling_regress:TooFewInputs','Not enough input arguments');
elseif nargin == 2
    kernel = 'gaussian';
end
if nargin < 4 || isempty(polyDegree)
    polyDegree = 6;
end
if nargin < 5
    span = [];
end

% set up kernel
[~,kappa2]= find_kernel(kernel);

% Check that matrix (X) and left hand side (Y) have compatible dimensions
[T,J] = size(X);
if size(Y,1) ~= T
    error('rolling_regress:InvalidData','The returns and factors have a different number of observations');
end

% Remove missing values, if any
wasnan = (any(isnan(X),2) | any(isnan(X),2));
havenans = any(wasnan);
if havenans
    Y(wasnan) = [];
    X(wasnan,:) = [];
    T = length(Y);
end
M = size(Y,2);
timeIndex = (1:T)';

%% First-stage: find a parametric estimator of alpha(t) and beta(t), and use them to compute the first-stage bandwidth

if isempty(span)
    
    timeX = [ones(T,1) exp(kron(1:polyDegree,log(timeIndex/T)))]; % time polynomial
    firstStageX = [timeX repmat(timeX,1,J).*reshape(repmat(X,polyDegree+1,1),T,J*(polyDegree+1))]; % alpha(t)+beta(t)*X
    [firstStageT,firstStageJ] = size(firstStageX); % size of alpha(t)+beta(t)*X
    [firstStageQ,firstStageR,firstStagePerm] = qr(firstStageX,0); % qr decomposition to later run least squares
    firstStageP = sum(abs(diag(firstStageR)) > max(firstStageT,firstStageJ)*eps(firstStageR(1))); % check for linear dependence of regressors
    if firstStageP < firstStageJ % if there is linear dependence
        warning(message('stats:regress:RankDefDesignMat'));
        firstStageR = firstStageR(1:firstStageP,1:firstStageP);
        firstStageQ = firstStageQ(:,1:firstStageP);
        firstStagePerm = firstStagePerm(1:firstStageP); % pick which coefficients are zero
    end
    
    firstStageB = zeros(firstStageJ,M); % fill in zeros in elements corresponding to rows of firstStageX that were thrown out.
    firstStageB(firstStagePerm',:) = firstStageR \ (firstStageQ'*Y); % fit polynomial by least squares
    firstStagePolyB = mat2cell(reshape(firstStageB(polyDegree+2:end,:),(1+polyDegree),J*M),polyDegree+1,ones(1,J*M)); % put in desired form
    firstStagePolyA = mat2cell(firstStageB(1:polyDegree+1,:),polyDegree+1,ones(1,M)); % put in desired form
    firstStageBetaTilde = cellfun(@(x) polyder(polyder(flipud(x))), firstStagePolyB,'UniformOutput',false); % take second derivative of estimate of beta(t)
    firstStageAlphaTilde = cellfun(@(x) polyder(polyder(flipud(x))), firstStagePolyA,'UniformOutput',false); % take second derivative of estimate of beta(t)
    % compute bias of estimates
    firstStageBiasB = cell2mat(cellfun(@(pp) polyval(pp,(1:T)'/T),firstStageBetaTilde,'UniformOutput',false));
    firstStageBiasB = sum(reshape(firstStageBiasB.^2,T*J,M))/T;
    firstStageBiasA = cell2mat(cellfun(@(pp) polyval(pp,(1:T)'/T),firstStageAlphaTilde,'UniformOutput',false));
    firstStageBiasA = sum(firstStageBiasA.^2)/T;
    % find bandwidth
    firstStageError = Y-firstStageX*firstStageB; % error of the polynomial fit
    sigmaKK = firstStageError'*firstStageError/T; % variance-covariance matrix of errors
    firstStageVariance = (kappa2/T)*norm(inv_posdef(cov(X,1)),'fro')*diag(sigmaKK)'; % estimate of the norm of the variance of beta(t)
    firstStageBandwidthB = T^(-1/5)*(firstStageVariance./firstStageBiasB).^(1/5); % first stage optimal bandwidth
    firstStageBandwidthA = T^(-1/5)*(firstStageVariance./firstStageBiasA).^(1/5); % first stage optimal bandwidth
    
    %% Second Stage: find a nonparametric estimator of alpha(t) and beta(t) using the bandwidth found in the first stage, and use the estimators of alpha(t) and beta(t) to compute the second-stage bandwidth
    
    % Find the optimal bandwidth for alpha(t)
    [secondStageA, ~, secondStageStats] = rolling_regress(Y,X,kernel,firstStageBandwidthA); % second stage estimate of beta(t)
    drop_alpha = any(isnan(secondStageA),2);
    if all(drop_alpha)
        hAlpha=nan;
        hBeta=nan;
        hLongRunAlpha=nan;
        hLongRunBeta=nan;
        return
    end
    secondStageA(drop_alpha,:) = [];
    secondStageStats.varErrors(drop_alpha) = [];
    for qq=1:M
        secondStageStats.resid{qq}(drop_alpha) = [];
        secondStageStats.varX{qq}(drop_alpha) = [];
    end
    [T,~] = size(X);
    T = T-sum(drop_alpha);
    timeIndex = (1:T)';
    
    secondStageA = mat2cell(secondStageA,T,ones(1,M)); % put in cell
    secondStagePoly = cellfun(@(bb) polyfit(timeIndex/T,bb,polyDegree),secondStageA,'UniformOutput',false); % fit polynomial for each beta(t)
    secondStageAlphaTilde = cellfun(@(x) polyder(polyder(x)), secondStagePoly,'UniformOutput',false); % take second derivative of estimate of beta(t)
    % compute bias of estimates
    secondStageBiasA = cell2mat(cellfun(@(pp) polyval(pp,(1:T)'/T),secondStageAlphaTilde,'UniformOutput',false));
    secondStageBiasA = sum(reshape(secondStageBiasA.^2,T,M))/T;
    % compute variance of residuals
    sigmaKK = cellfun(@(xx) diag(xx),secondStageStats.varErrors,'UniformOutput',false);
    sigmaKK = mat2cell(num2cell(cell2mat(sigmaKK)'),T,ones(1,M));
    % compute the (time integrated) inverse of the variance of X
    invVarX = cellfun(@(outer) cellfun(@(inner) inv_posdef(inner),outer,'UniformOutput',false),secondStageStats.varX,'UniformOutput',false);
    integratedVar = cellfun(@(outerVarX,outerSigmaKK)  cellfun(@(vv,ss) vv*ss,outerVarX,outerSigmaKK,'UniformOutput',false),invVarX,sigmaKK,'UniformOutput',false);
    integratedVar = cellfun(@(xx) norm((kappa2/T)*sum(cat(3,xx{:}),3),'fro'),integratedVar,'UniformOutput',false);
    % second stage optimal bandwidth
    hAlpha = T^(-1/5)*(cell2mat(integratedVar)./secondStageBiasA).^(1/5);
    adjustment = 2/15; % must be smaller than 1/3
    hLongRunAlpha = hAlpha*T^(-adjustment);
    
    if nargout >=2
        % Find the optimal bandwidth for beta(t)
        [~, secondStageB, secondStageStats] = rolling_regress(Y,X,kernel,firstStageBandwidthB); % second stage estimate of beta(t)
        drop_beta = any(isnan(cell2mat(secondStageB)),2);
        if all(drop_beta)
            hAlpha=nan;
            hBeta=nan;
            hLongRunAlpha=nan;
            hLongRunBeta=nan;
            return
        end
        secondStageStats.varErrors(drop_beta) = [];
        for qq=1:M
            secondStageB{qq}(drop_beta,:) = []; 
            secondStageStats.resid{qq}(drop_beta) = [];
            secondStageStats.varX{qq}(drop_beta) = [];
        end
        [T,~] = size(X);
        T = T-sum(drop_beta);
        timeIndex = (1:T)';
        
        secondStageB = mat2cell(cell2mat(secondStageB),T,ones(1,J*M)); % separate betas
        secondStagePoly = cellfun(@(bb) polyfit(timeIndex/T,bb,polyDegree),secondStageB,'UniformOutput',false); % fit polynomial for each beta(t)
        secondStageBetaTilde = cellfun(@(x) polyder(polyder(x)), secondStagePoly,'UniformOutput',false); % take second derivative of estimate of beta(t)
        % compute bias of estimates
        secondStageBiasB = cell2mat(cellfun(@(pp) polyval(pp,(1:T)'/T),secondStageBetaTilde,'UniformOutput',false));
        %secondStageBiasB = cell2mat(cellfun(@(pp) polyval(pp,(1:T)'),secondStageBetaTilde,'UniformOutput',false));
        secondStageBiasB = sum(reshape(secondStageBiasB.^2,T*J,M))/T;
        % compute variance of residuals
        sigmaKK = cellfun(@(xx) diag(xx),secondStageStats.varErrors,'UniformOutput',false);
        sigmaKK = mat2cell(num2cell(cell2mat(sigmaKK)'),T,ones(1,M));
        % compute the (time integrated) inverse of the variance of X
        invVarX = cellfun(@(outer) cellfun(@(inner) inv_posdef(inner),outer,'UniformOutput',false),secondStageStats.varX,'UniformOutput',false);
        integratedVar = cellfun(@(outerVarX,outerSigmaKK)  cellfun(@(vv,ss) vv*ss,outerVarX,outerSigmaKK,'UniformOutput',false),invVarX,sigmaKK,'UniformOutput',false);
        integratedVar = cellfun(@(xx) norm((kappa2/T)*sum(cat(3,xx{:}),3),'fro'),integratedVar,'UniformOutput',false);
        
        % second stage optimal bandwidth
        hBeta = T^(-1/5)*(cell2mat(integratedVar)./secondStageBiasB).^(1/5);
        adjustment = 2/15; % must be smaller than 1/3
        hLongRunBeta = hBeta*T^(-adjustment);
    end
    
else
    if strcmpi(kernel,'lag one-sided exponential') || strcmpi(kernel,'one-sided exponential')
        h = span/(T*log(2));
    elseif strcmpi(kernel,'gaussian')
        h = span/T;
    elseif strcmpi(kernel,'lag one-sided-uniform') ||  strcmpi(kernel,'one-sided-uniform')
        h = 2*span/size(Y,1);
    else
        h = 1;
    end
    hAlpha= h;
    hBeta= h;
    hLongRunAlpha= h;
    hLongRunBeta= h;
end

end