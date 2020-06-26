function [pcs vals share_var xs_mean xs_std coeff]= PCA_unbalanced(X,k,threshold,max_iter,max_runs,constant_flag,zscore_flag)
%PCA_UNBALANCED Computes principal components of a matrix with missing
%observations (such as an unbalanced panel)
%
%   PCS = PCA_UNBALANCED(X) returns the T-by-1 vector corresponding to the first principal
%   component of the T-by-M matrix X. Missing values of X are indicated by
%   NaN's. If M=1, PCS is X minus its mean. 
%
%   PCS = PCA_UNBALANCED(X,k) returns the T-by-k matrix corresponding to the first k principal
%   components of the T-by-M matrix X.
%
%   PCS = PCA_UNBALANCED(X,k,threshold) also specifies a threshold for
%   convergence in the algorithm (if missing, default is 1e-4).
%
%   PCS = PCA_UNBALANCED(X,k,threshold,max_iter) also specifies the maximum number
%   of itrations attempted for convergence (if missing, default is
%   100,000).
%
%   PCS = PCA_UNBALANCED(X,k,threshold,max_iter,max_runs) also specifies
%   the maximum number of runs, where the first run uses a random initial
%   guess for the principal components and the subsuquent runs use the
%   result of the previous run as an initial guess (if missing, defaul is
%   5).
%
%   PCS = PCA_UNBALANCED(X,k,threshold,max_iter,max_runs,constant_flag) also specifies
%   whether to include a constant in the regression used to "fill in"
%   missing observations. A value of 0 (the default) means no constant, a
%   value of 1 means include the constant.
%
%   PCS = PCA_UNBALANCED(X,k,threshold,max_iter,max_runs,constant_flag,zscore_flag)
%   normalizes X if zscore is 1 and does not normalize if zscore is 0 (the
%   default)
%
%   [PCS VALS] = PCA_UNBALANCED(...) also returns the 1-by-M matrix of
%   eigenvalues of the covariance of X.
%
%   [PCS VALS SHARE_VAR] = PCA_UNBALANCED(...) also returns the 1-by-M vector with the share of
%   the variance of X explained by each of the principal components.
%
%   [PCS VALS SHARE_VAR XS_MEAN] = PCA_UNBALANCED(...) also returns the T-by-1 vector
%   with the cross-sectional mean of the non-missing observations of X.
%
%   [PCS VALS SHARE_VAR XS_MEAN XS_STD] = PCA_UNBALANCED(...) also returns the T-by-1 vector
%   with the cross-sectional standard deviation of the non-missing observations of X.

%% set up default parameters if not specified by user
if nargin<7
    zscore_flag = 0;
end
if nargin<6
    constant_flag = 0;
end
if nargin<5
    max_runs=5;
end
if nargin<4
   max_iter=100000; 
end
if nargin<3
    threshold=1e-4;
end
if nargin<2
    k=1;
end     

%% get PCA

% replace rows with consensus when too few are forecasting
% fcnans=sum(isnan(X),2); % number of missing observations
% if sum(fcnans>2*size(Y,2)/3)>0
%   X(fcnans>2*size(X,2)/3,:)=globalConsensus(fcnans>2*size(X,2)/3)*ones(1,size(X,2));
% end

% cross-sectional mean and std dev
xs_mean=nanmean(X,2);
xs_std=nanstd(X,0,2);

% if input is a single time-series (M=1) then PCS is X minus its mean
if size(X,2)==1
    nanX = isnan(X);
    pcs=X-nanmean(X);
    pcs(nanX)=NaN;
    vals = 1;
    share_var = 1;
    return
end

% initialize variables
F = cell(max_runs,1);
F_old = cell(max_runs,1);
if constant_flag
    lambda=NaN(k+1,size(X,2)); %third dimension is 2 for regressing on factor + constant
else
    lambda=NaN(k,size(X,2)); 
end

for runn = 1:max_runs
    rng(1); % seed random number generator for same results across runs
    F{runn}=randn(size(X(:,1:k))); % random initial guess
    F{runn}=F{runn}-repmat(nanmean(F{runn},1),size(F{runn},1),1); % de-meaning random starting values
    % F{runn}(:,1)=nanmean(Y,2); % mean is initial guess
    F_old{runn}=zeros(size(F{runn}));
    Xe=X; iter=0;
    
    % find PC's
    while threshold<=sqrt(nanmean(nanmean((F{runn}-F_old{runn}).^2))) && iter<=max_iter
        iter=iter+1;
        for i=1:size(X,2)
            if constant_flag
                lambda(:,i)=regress(Xe(~isnan(X(:,i)),i),[F{runn}(~isnan(X(:,i)),:) ones(size(F{runn}(~isnan(X(:,i)),1)))]);
                Xe(isnan(X(:,i)),i)=lambda(:,i)' * [F{runn}(isnan(X(:,i)),:) ones(size(F{runn}(isnan(X(:,i)),1)))]'; % transposing here because I store lambda and F incorrectly (transposed) initially
            else
                lambda(:,i)=regress(Xe(~isnan(X(:,i)),i),F{runn}(~isnan(X(:,i)),:));
                Xe(isnan(X(:,i)),i)=lambda(:,i)' * F{runn}(isnan(X(:,i)),:)'; % transposing here because I store lambda and F incorrectly (transposed) initially
            end

        end
        F_old{runn}=F{runn};
        if zscore_flag == 1
            [coeff,PCs]=pca(zscore(Xe));
            coeff = coeff(:,1:k);
        else
            [coeff,PCs]=pca(Xe);
            coeff = coeff(:,1:k);    
        end
        F{runn}=PCs(:,1:k);
        
        % normalize weights to 1
        F{runn} = F{runn}./repmat(sum(coeff(:,1:k)),size(F{runn},1),1);
        
        % for debugging purposes
        % figure(2);
        % plot(zscore([consensus F{runn}(:,1) F_old{runn}(:,1)]));
        % pause;
        % if ~mod(iter,500), disp(iter); end;
        
        sqrt(nanmean(nanmean((F{runn}-F_old{runn}).^2)));
    end
    
    % store shares in terms of eigenvalues
    [~,vals] = eig(cov(Xe));
    [vals,~] = sort(diag(vals),'descend');
    share_var = vals./sum(vals);
    
%     outputting convergence status
%     if iter<max_iter
%         fprintf('convergence reached in %i steps \n',iter);
%     else
%         fprintf('WARNING: max iterations reached \n');
%     end
end % for runn


pcs = F{end};

% plots
%{
        figure(1);
        subplot(k,1,1);
        plot([F{1}(:,1) F{2}(:,1) F{3}(:,1) F{4}(:,1) F{5}(:,1)]);% F{6}(:,1) F{7}(:,1) F{8}(:,1) F{9}(:,1) F{10}(:,1)]);
        if k>1,
            subplot(k,1,2);
            plot([F{1}(:,2) F{2}(:,2) F{3}(:,2) F{4}(:,2) F{5}(:,2)]);% F{6}(:,2) F{7}(:,2) F{8}(:,2) F{9}(:,2) F{10}(:,2)]);
            if k>2,
            subplot(k,1,3);
            plot([F{1}(:,3) F{2}(:,3) F{3}(:,3) F{4}(:,3) F{5}(:,3)]);% F{6}(:,3) F{7}(:,3) F{8}(:,3) F{9}(:,3) F{10}(:,3)]);
            end,
        end,
%}
% figure(1);
% plot(zscore([consensus F{1}(:,1) F{2}(:,1) F{3}(:,1) F{4}(:,1) F{5}(:,1)]));% F{6}(:,1) F{7}(:,1) F{8}(:,1) F{9}(:,1) F{10}(:,1)]);
% eval(sprintf('title(''Q%i %s PC1 and Consensus'');',q,variable_array{v}));
% eval(sprintf('saveas(1,''Q%i_%s_fill1.fig'');',q,variable_array{v}));
% if k>1,
%     figure(2);
%     plot(zscore([volatility F{1}(:,2) F{2}(:,2) F{3}(:,2) F{4}(:,2) F{5}(:,2)]));% F{6}(:,1) F{7}(:,1) F{8}(:,1) F{9}(:,1) F{10}(:,1)]);
%     eval(sprintf('title(''Q%i %s PC2 and volatility'');',q,variable_array{v}));
%     eval(sprintf('saveas(2,''Q%i_%s_fill2.fig'');',q,variable_array{v}));
% end


