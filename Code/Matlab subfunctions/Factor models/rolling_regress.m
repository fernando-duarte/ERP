function [alpha beta stats] = rolling_regress(Y,X,kernel,bandwidth)
% ROLLING_REGRESS estimates coefficients of rolling least-squares regressions.
%
%   ALPHA = ROLLING_REGRESS(Y,X) gives the intercepts of rolling weighted least-squares
%   regressions of the columns of Y on X using weights given by KERNEL.
%   X is a T-by-J matrix of T observations of
%   J independent variables. Y is a T-by-M matrix of T observations
%   of M dependent variables. ALPHA is a T-by-M matrix of
%   intercepts. Note: X should not include a column of ones.
%
%   [ALPHA BETA] = ROLLING_REGRESS(Y,X) also returns the a 1-by-M cell, each
%   containing a T-by-J matrix of regression coefficients.
%
%   [ALPHA BETA STATS] = ROLLING_REGRESS(Y,X) returns a structure with the
%   following fields:
%
%       STATS.KERNEL    : a string or function handle of the kernel used for
%                           the regression (it copies the variable KERNEL.)
%       STATS.BANDWIDTH : a 1-by-M vector of bandwidths used for the
%                           regression.
%       STATS.RESID     : a 1-by-M cell, where each cell contains a T-by-1
%                           vector of regression residuals.
%       STATS.VARX      : a 1-by-M cell, where each cell has T subcells.
%                           Each of the subcells has a J-by-J matrix with
%                           the estimated variance-covariance matrix of the factors
%                           X.
%       STATS.COVX      : a T-by-1 cell, where each cell has a JM-by-JM matrix with
%                           the estimated variance-covariance matrix of the
%                           factors X. The cell COVX(i,j) computes
%                           covariances using the bandwidths of assets
%                           i and j.
%       STATS.COVBETA   : a T-by-1 cell, where each cell has a JM-by-JM matrix with
%                           the estimated variance-covariance matrix of the
%                           regression coefficients BETA.
%       STATS.VARERRORS : a 1-by-T cell, where each cell has an M-by-M
%                           matrix that is the estimated variance-covariance matrix of the
%                           regression residuals.
%
%   [...] = ROLLING_REGRESS(Y,X,KERNEL) specifies the KERNEL. If empty or not specified,
%   the default is 'gaussian'. Other options are 'one_sided_uniform','two_sided_uniform'
%   'one_sided_exponential' or a used-supplied function handle (the KERNEL must integrate to one).
%
%   [...] = ROLLING_REGRESS(Y,X,KERNEL,BANDWIDTH) specifies the
%   KERNEL's BANDWIDTH. If BANDWIDTH is a number, the same BANDWIDTH is
%   used for all columns of Y. If BANDWIDTH is a 1-by-M vector, then the
%   kth element of the vector is used for the kth column of Y. If BANDWIDTH
%   is missing or not defined, it is assumed to be equal to a 1-by-M vector
%   of ones.
%
%   References:
%      [1] Ang, A. and Kristensen, D. (2011) Testing Conditional Factor
%      Models
%
%   Copyright 2011 Fernando M. Duarte
%   Revision: 0.1   Date: 02/06/2012

if  nargin < 2
    error('rolling_regress:TooFewInputs','Not enough input arguments');
elseif nargin == 2
    kernel = [];
    bandwidth = 1;
elseif nargin == 3
    bandwidth = 1;
end

% Check that matrix (X) and left hand side (Y) have compatible dimensions
[T,J] = size(X);
if size(Y,1) ~= T
    error('rolling_regress:InvalidData','The returns and factors have a different number of observations');
end

% find dimension of Y
M = size(Y,2);

% if any of the bandwidths is nan, return nan's for the function output
if any(isnan(bandwidth)) || any(isinf(bandwidth)) || any(bandwidth<=0)
    alpha = nan(T,M);
    beta = cell(1,M);
    beta = cellfun(@(x) nan(T,J),beta,'UniformOutput',false);
    stats.kernel = kernel;
    stats.bandwidth = repmat(bandwidth,1,M);
    stats.resid = cell(1,M);
    stats.resid = cellfun(@(x) nan(T,1),stats.resid,'UniformOutput',false);
    varX = cell(T,1);
    varX = cellfun(@(x) nan(J,J),varX,'UniformOutput',false);
    stats.varX = mat2cell(repmat(varX,1,M),T,ones(1,M));
    stats.varErrors = cell(1,T);
    stats.varErrors = cellfun(@(x) nan(M,M),stats.varErrors,'UniformOutput',false);
    return
end

% set up kernel
[kernel_fun,kappa2,~,one_sided_flag] = find_kernel(kernel);

% Remove missing values, if any
wasnan = (any(isnan(X),2) | any(isnan(Y),2));
havenans = any(wasnan);
if havenans
    Y(wasnan,:) = [];
    X(wasnan,:) = [];
    T = length(Y);
end

% if looking at a one-sided kernel, use start_time observations to
% initialize kernel
if one_sided_flag == false
    start_time = 1;
else
    start_time = max(floor(0.05*T),6);
end

% run regressions
timeIndex = (1:T)';

if all(abs(bandwidth-bandwidth(1))<eps) %same bandwidth for all regressions, use faster code
    
    % all regressions use the same bandwidth, reduce vector to a number
    bandwidth = bandwidth(1);
    
    % initialize variables
    alpha = nan(T,M);
    beta_temp = nan(T,J*M);
    varX = cellfun(@(x) nan(J,J),cell(T,1),'UniformOutput',false);
    stats.varErrors = cellfun(@(y) nan(M,M), cell(1,T),'UniformOutput',false);
    
    for t = start_time:T % loop through time to run rolling regressions
        kernelWeight = kernel_fun((timeIndex-t)/(bandwidth*T))/(bandwidth*T);

        sqrtKernel = sqrt(kernelWeight);
        weightedX = repmat(sqrtKernel,1,J+1).*[ones(T,1) X];
        weightedY = repmat(sqrtKernel,1,M).*Y;
        
        % Use the rank-revealing QR to remove dependent columns of weightedX.
        [Q,R,perm] = qr(weightedX,0);
        p = sum(abs(diag(R)) > max(T,J+1)*eps(R(1)));
        if (p < J+1) && (rank(X) > rank(weightedX))
            warning('rolling_regress:LinearDependent','A small bandwith made kernel regression rank deficient. Increase the bandwidth.');
        else
            if (p < J+1) && (rank(X) <= rank(weightedX))
                warning(message('stats:regress:RankDefDesignMat'));
                R = R(1:p,1:p);
                Q = Q(:,1:p);
                perm = perm(1:p);
            end
            
            % Compute the LS coefficients, filling in zeros in elements corresponding
            % to rows of weightedX that were thrown out.
            b = zeros(J+1,M);
            b(perm',:) = R \ (Q'*weightedY);
            
            % Separate intercept and slope
            alpha(t,:) = b(1,:);
            beta_temp(t,:) = reshape(b(2:end,:),1,J*M);
            
            if nargout == 3
                % find variance of X
                normalization = sum(kernelWeight);
                driftX = kernelWeight'*X/normalization;
                demeanedX = bsxfun(@minus,X,driftX);
                weightedDemeanedX = repmat(sqrtKernel,1,J).*demeanedX;
                varX{t} = weightedDemeanedX'*weightedDemeanedX/normalization;
                
                % find residuals' variance-covariance matrix
                weightedResid = weightedY - weightedX*b;
                stats.varErrors{1,t} = (weightedResid'*weightedResid)./(sqrtKernel'*sqrtKernel);
            end
        end
    end %end time iteration
    
    % put in desired form, 1-by-M cell for coefficients
    beta = mat2cell(beta_temp,T,repmat(J,1,M));
    
    if nargout == 3 % finish computing statistics
        
        %put varX in desired size format
        stats.varX = mat2cell(repmat(varX,1,M),T,ones(1,M));
        
        % find overall residual
        stats.resid = cellfun(@(yy,aa,bb) yy-sum([ones(T,1) X].*[aa bb],2),mat2cell(Y,T,ones(1,M)),mat2cell(alpha,T,ones(1,M)),beta,'UniformOutput',false);
        
        % find covariance of X
        stats.covX = cellfun(@(x) kron(x,ones(M,M)),varX,'UniformOutput',false);
        
        % find variance-covariance matrix of beta
        varX_inverse = cellfun(@inv_posdef,varX,'UniformOutput',false);
        stats.covBeta = cellfun(@(lambda,sigma) kappa2*kron(lambda,sigma)/(T*bandwidth),varX_inverse,stats.varErrors','UniformOutput',false);
        
        % put kernel and bandwidth into stats
        stats.kernel = kernel;
        stats.bandwidth = repmat(bandwidth,1,M);   
    end
    
else %different bandwidths for different regressions
    
    % initialize variables
    b = cell(T,1);driftX = cell(T,1);
    b = cellfun(@(x) cell(1,M),b,'UniformOutput',false);
    b = cellfun(@(x) cellfun(@(y) nan(J+1,1),x,'UniformOutput',false),b,'UniformOutput',false);
    varX = cellfun(@(x) cell(1,M),cell(T,1),'UniformOutput',false);
    varX = cellfun(@(x) cellfun(@(y) nan(J,J),x,'UniformOutput',false),varX,'UniformOutput',false);
    stats.covX = cellfun(@(y) nan(J*M,J*M),cell(T,1),'UniformOutput',false);
    stats.covBeta = cellfun(@(y) nan(J*M,J*M), cell(T,1),'UniformOutput',false);
    stats.varErrors = cellfun(@(y) nan(M,M), cell(1,T),'UniformOutput',false);
    
    for t = start_time:T % loop through time to run rolling regressions
        
        %create kernel weights
        kernelWeight = kernel_fun((timeIndex-t)*(1./(bandwidth*T))).*repmat(1./(bandwidth*T),T,1);
        sqrtKernel = sqrt(kernelWeight);
        cellKernel = mat2cell(sqrtKernel,T,ones(1,M));
        
        %create weighted X and Y, then put them in a cell
        weightedX = cellfun(@(weights) repmat(weights,1,J+1).*[ones(T,1) X],cellKernel,'UniformOutput',false);
        weightedY = mat2cell(sqrtKernel.*Y,T,ones(1,M));
        
        % Compute the least-squares coefficients, filling in zeros in elements corresponding
        % to rows of weightedX that were thrown out.
        [Q,R,perm] = cellfun(@qr_cell,weightedX,'UniformOutput',false); %qr decomposition
        
        if ~any(cellfun(@isempty,R))
        % check linear dependency
        p = cell2mat(cellfun(@(rr) sum(abs(diag(rr)) > max(T,J+1)*eps(rr(1)))<J+1,R,'UniformOutput',false));
        rank_weightedX = cellfun(@(x) rank(x),weightedX);
        if (max(p) ~= 0) && any(rank(X) > rank_weightedX) % rank deficient
            warning('rolling_regress:LinearDependent','Small bandwith made kernel regression rank deficient. Increase the bandwidth');
        elseif (max(p) ~= 0) && any(rank(X) <= rank_weightedX)  % some linear dependence
            warning('rolling_regress:LinearDependent','Some of the regressors are linearly dependent');
            % you could alternatively set coefficients to zero for linearly
            % dependent columns, as done in the case with length(bandwidth)==1
            b{t} = cell(1,M);
            b{t} = cellfun(@(x) nan(J+1,1),b{t},'UniformOutput',false);
            varX{t} = cell(1,M);
            varX{t} = cellfun(@(x) nan(J,J), varX{t} ,'UniformOutput',false);
            stats.varErrors{1,t} = nan(M,M);
        else % no linear dependence in X
            
            % get coefficients
            b{t} = cellfun(@(qq,rr,yy) rr \ (qq'*yy),Q,R,weightedY,'UniformOutput',false); %least-squares
            [~, invperm] =  cellfun(@(order) sort(order),perm,'UniformOutput',false); % reverse order of qr decomposition
            b{t} = cellfun(@(bb,order) bb(order'),b{t},invperm,'UniformOutput',false); %re-arrange in order of qr decomposition
            
            % compute some statistics
            if nargout == 3 
                
                % find variance of X
                normalization = num2cell(sum(kernelWeight));
                driftX{t} = cellfun(@(weights,nn) weights'*X/nn,mat2cell(kernelWeight,T,ones(1,M)),normalization,'UniformOutput',false);
                demeanedX = cellfun(@(mm) bsxfun(@minus,X,mm),driftX{t},'UniformOutput',false);
                weightedDemeanedX = cellfun(@(weights,xx) repmat(weights,1,J).*xx,cellKernel,demeanedX,'UniformOutput',false);
                varX{t} = cellfun(@(xx,nn) xx'*xx/nn,weightedDemeanedX,normalization,'UniformOutput',false);
                
                % find covariance of X
                weightedDemeanedX_cov = reshape(cell2mat(weightedDemeanedX'),T,J*M);
                cov_normalization = sqrt(cell2mat(normalization)'*cell2mat(normalization));
                stats.covX{t} = (weightedDemeanedX_cov'*weightedDemeanedX_cov)./repmat(cov_normalization,J,J);
                
                % find residuals' variance-covariance matrix
                weightedResid = cell2mat(cellfun(@(yy,xx,bb) yy-xx*bb,weightedY,weightedX,b{t},'UniformOutput',false)); %re-arrange in order of qr decomposition
                stats.varErrors{1,t} = (weightedResid'*weightedResid)./cov_normalization;
                
                % find variance-covariance matrix of beta
                kernelMatrix = sqrt(bandwidth'*bandwidth);
                stats.covBeta{t} = kappa2*inv_repeat_rows(stats.covX{t}).*repmat(stats.varErrors{1,t}./(T*kernelMatrix),J,J);               
            end
        end
        end
    end %end time iteration
    
    % Separate intercept and slope, put in desired form, 1-by-M cell for
    % coefficients and T-by-M matrix for slopes
    b = mat2cell(cell2mat(cellfun(@cell2mat,b,'UniformOutput',false)),T*(J+1),ones(1,M)); % group portfolios
    b = cellfun(@(bb) reshape(bb,J+1,T)',b,'UniformOutput',false); % put time dimension on rows
    alpha = cell2mat(cellfun(@(aa) aa(:,1),b,'UniformOutput',false)); % get intercept
    beta = cellfun(@(aa) aa(:,2:end),b,'UniformOutput',false); % get slopes
    
    % finish computing statistics
    if nargout == 3 
        % put varX in desired size format
        stats.varX = mat2cell(mat2cell(cell2mat(cellfun(@(xx) cell2mat(xx),varX,'UniformOutput',false)),J*ones(T,1),J*ones(1,M)),T,ones(1,M));
        
        % find overall residual
        stats.resid = cellfun(@(yy,bb) yy-sum([ones(T,1) X].*bb,2),mat2cell(Y,T,ones(1,M)),b,'UniformOutput',false);
        
        % put kernel and bandwidth into stats
        stats.kernel = kernel;
        stats.bandwidth = bandwidth;
    end
end

% sub-functions

    function [Q,R,perm] = qr_cell(regressor)
        % Use the rank-revealing QR to remove dependent columns of weightedX.
        [Q,R,perm] = qr(regressor,0);
        [TT,JJ] = size(regressor);
        p = sum(abs(diag(R)) > max(TT,JJ)*eps(R(1)));
        if p < JJ
            warning(message('stats:regress:RankDefDesignMat'));
            R = R(1:p,1:p);
            Q = Q(:,1:p);
            perm = perm(1:p);
        end
    end % end qr_cell function

    function Y = inv_repeat_rows(X)
        % computes the inverse of a matrix with repeated rows by first
        % deleting rows and then assuming repeated rows have the same
        % entries in the inverse -- coincides with \ operator when X has no
        % repeated rows
        
        %Xround = round(X*1e14)/1e14; % round to 14 significant figuers to eliminate floating-point errors
        Xround = cast(X,'single');% round to single precision to eliminate floating-point errors
        [~, index, invIndex]=unique(Xround,'rows');
        small_inv = inv_posdef(X(index,index)); % X(index,index)\eye(size(index,1),size(index,1));
        Y = small_inv(invIndex,invIndex);
    end % end inv_repeat_rows function

end % end main function