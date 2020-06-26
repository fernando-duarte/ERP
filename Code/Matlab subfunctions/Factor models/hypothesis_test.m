function [pvalue point_estimate conf_intvl st] = hypothesis_test(alpha,beta,stats,test,significance)
% HYPOTHESIS_TEST tests hypothesis about the coefficients obtained by rolling_regress
%
%   PVALUE = HYPOTHESIS_TEST(ALPHA,BETA,STATS,TEST) gives p-value for the
%   hypothesis test defined by the variable TEST using the estimates of
%   ALPHA, BETA and STATS returned by rolling_regress. ALPHA is a T-by-M matrix of
%   intercepts. BETA is a 1-by-M cell, each containing a T-by-J matrix of regression
%   coefficients. STATS is a structure returned by the function rolling_regress.
%   The different hypothesis tests are given by the string TEST:
%
%       'Conditional Beta'   : Test the null that each beta(t) is equal to zero
%                               against the alternative that beta(t) is not zero
%                               for each beta(t) individually
%       'Conditional Alpha'  : Test the null that alpha(t) is equal to zero
%                               against the alternative that alpha(t) is not zero
%                               for each alpha(t) individually
%       'LR Beta'            : Test the null that the long-run beta is
%                               equal to zero against the alternative that
%       'All LR Beta'        : Test the null that all long-run betas are
%                               jointly equal to zero against the alternative that
%                              any of them is not zero
%       'LR Alpha'           : Test the null that the long-run alpha is
%                               equal to zero against the alternative that
%                               the long-run alpha is not zero
%       'All LR Alpha'       : Test the null that all long-run alphas are
%                               jointly equal to zero against the
%                               alternative that any of them is not zero
%       'Constant Beta'      : Test the null that beta(t) is constant and
%                               equal to its time average
%       'Constant Alpha'     : Test the null that alpha(t) is constant and
%                               equal to its time average
%       'Constant Alpha Zero': Test the null that alpha(t) is constant and
%                               equal to zero
%
%   PVALUE = HYPOTHESIS_TEST(ALPHA,BETA,STATS,TEST,SIGNIFICANCE) defines
%   the statistical significance of the test. SIGNIFICANCE must be a number
%   between zero and one. If not specified, the default is 0.05 (ie. 5%)
%
%   [PVALUE POINT_ESTIMATE] = HYPOTHESIS_TEST(...) gives the point
%   estimates of the estimator used in the test.
%
%   [PVALUE POINT_ESTIMATE CONF_INTVL] = HYPOTHESIS_TEST(...) if SIGNIFICANCE
%   is not specified, gives the 95% confidence interval. If SIGNIFICANCE is
%   specified, gives the 100*(1-significance) confidence interval.
%
%   [PVALUE POINT_ESTIMATE CONF_INTVL ST] = HYPOTHESIS_TEST(...) returns
%   a structure with additional statistics. ST.TESTSTAT is the
%   test statistic used to run the test. ST.CRITICAL99,
%   ST.CRITICAL95, ST.CRITICAL90 are the 99%, 95% and 90%
%   critical values of the test statistic. ST.H is 1 if the null is
%   rejected, 0 if the test fails to reject the null.
%
%   References:
%      [1] Ang, A. and Kristensen, D. (2011) Testing Conditional Factor
%      Models
%
%   Copyright 2011 Fernando M. Duarte
%   Revision: 0.1   Date: 02/03/2012
%% Preliminary checks and setup

if nargin < 4
    error('hypothesis_test:TooFewInputs','Not enough input arguments');
elseif nargin == 4
    significance = 0.95;
    if ~ischar(test)
        error('hypothesis_test:NoTest','Variable TEST must be a string');
    end
elseif nargin > 4
    if significance < 0 || significance > 1
        error('hypothesis_test:BadSignificance','Significance must be between zero and one');
    end
end

% set up kernel
[~, kappa2,convK2]= find_kernel(stats.kernel);

% find sizes
[T,M] = size(alpha);
J = size(beta{1},2);
%% Run tests
switch lower(test)
    case 'conditional beta'
        point_estimate = beta;
        varBeta = cell(1,M);
        for m = 1:M % loop through left hand side variables
            % compute variance of beta(t)
            varBeta{m} = cellfun(@(varX,varErrors) kappa2*diag(inv_posdef(varX)*varErrors(m,m))'/(T*stats.bandwidth(m)),stats.varX{m},stats.varErrors','UniformOutput',false);
        end
        % compute the tstats
        tstatBeta = cellfun(@(bb,vv) bb./sqrt(cell2mat(vv)),point_estimate,varBeta,'UniformOutput',false);
        % do the test
        [st.h,pvalue,~,st.teststat] = cellfun(@(zz) ztest(zz(:),0,1,1-significance,'both',2),tstatBeta,'UniformOutput',false);
        % get confidence interval
        interval = norminv([significance/2 1-significance/2]);
        conf_intvl=cellfun(@(bb,vv) [bb+interval(1)*sqrt(cell2mat(vv)); bb+interval(2)*sqrt(cell2mat(vv))],beta,varBeta,'UniformOutput',false);    
        % put in desired shape
        pvalue = cellfun(@(xx) reshape(xx,T,J),pvalue,'UniformOutput',false);
        conf_intvl=cellfun(@(x) mat2cell(reshape(x,T,2*J),T,repmat(2,J,1)),conf_intvl,'UniformOutput',false);
        st.h = cellfun(@(xx) reshape(xx,T,J),st.h,'UniformOutput',false);
        st.teststat = cellfun(@(xx) reshape(xx,T,J),st.teststat,'UniformOutput',false);
        % find critical values
        st.critical99 = norminv(1-0.01/2);
        st.critical95 = norminv(1-0.05/2);
        st.critical90 = norminv(1-0.1/2);
    case 'all conditional beta'
        % no point value or confidence interval since it's a test of joint
        % significance
        point_estimate = nan;
        conf_intvl=nan;
        % find F-statistic
        betaCell = mat2cell(reshape(cell2mat(beta'),T,J*M),ones(T,1),J*M);
        Fstat = cell2mat(cellfun(@(x,v) x*inv_posdef(v)*x',betaCell,stats.covBeta,'UniformOutput',false));
        st.teststat = Fstat;
        % get p-value
        pvalue = 1-chi2cdf(Fstat,M*J);
        % perform test (h=1 means reject the null)
        st.h = pvalue<=significance;
        %find critical values
        st.critical99 = chi2inv(1-0.01,M*J);
        st.critical95 = chi2inv(1-0.05,M*J);
        st.critical90 = chi2inv(1-0.1,M*J);          
    case 'conditional alpha'
        point_estimate = alpha;
        % compute variance of alpha(t)
        varAlpha = cell2mat(cellfun(@(varErrors)  kappa2*diag(varErrors)'./(T*stats.bandwidth),stats.varErrors','UniformOutput',false));
        % compute the tstats
        tstatAlpha = alpha./sqrt(varAlpha);
        % do the test
        [st.h,pvalue,~,st.teststat] = ztest(tstatAlpha(:),0,1,1-significance,'both',2);
        % get confidence interval
        interval = norminv([significance/2 1-significance/2]);
        conf_intvl=[alpha+interval(1)*sqrt(varAlpha); alpha+interval(2)*sqrt(varAlpha)];    
        % put in desired shape
        pvalue = reshape(pvalue,T,M);
        conf_intvl = mat2cell(reshape(conf_intvl,T,2*M),T,repmat(2,M,1));
        st.h = reshape(st.h,T,M);
        st.teststat = reshape(st.teststat,T,M);
        % find critical values
        st.critical99 = norminv(1-0.01/2);
        st.critical95 = norminv(1-0.05/2);
        st.critical90 = norminv(1-0.1/2);
    case 'lr beta'
        % compute long-run beta
        point_estimate = cellfun(@(bb) mean(bb),beta,'UniformOutput',false);
        % find variance
        varBeta_t = cell(1,M);
        for m = 1:M % loop through left hand side variables
            % compute variance of beta(t) and the long-run beta
            varBeta_t{m} = cellfun(@(varX,varErrors) diag(inv_posdef(varX)*varErrors(m,m))'/T,stats.varX{m},stats.varErrors','UniformOutput',false);
        end
        varBetaLR = cellfun(@(vv) mean(cell2mat(vv)),varBeta_t,'UniformOutput',false);
        % compute the tstats
        tstatBetaLR = cellfun(@(bb,vv) bb./sqrt(vv),point_estimate,varBetaLR,'UniformOutput',false);
        % do the test
        [st.h,pvalue,~,st.teststat] = cellfun(@(zz) ztest(zz(:),0,1,1-significance,'both',2),tstatBetaLR,'UniformOutput',false);
        % get confidence interval
        interval = norminv([significance/2 1-significance/2]);
        conf_intvl=cellfun(@(bb,vv) [bb+interval(1)*sqrt(vv); bb+interval(2)*sqrt(vv)],point_estimate,varBetaLR,'UniformOutput',false);    
        % put in desired form
        conf_intvl=cellfun(@(x) mat2cell(reshape(x,1,2*J),1,repmat(2,J,1)),conf_intvl,'UniformOutput',false);       
        % find critical values
        st.critical99 = norminv(1-0.01/2);
        st.critical95 = norminv(1-0.05/2);
        st.critical90 = norminv(1-0.1/2);
    case 'all lr beta'
        point_estimate = nan;
        
        % compute long-run beta
        betaLR = cellfun(@(bb) mean(bb),beta,'UniformOutput',false);
        betaLR = cell2mat(betaLR);
        
        % find variance-covariance matrix of LR betas
        varBeta = cell(M,M);
        for k1=1:M
            for k2=1:M
                tempVarBeta_t = cell(1,T);
                for t=1:T
                    tempVarBeta_t{1,t} = inv_posdef(stats.covX{k1,k2}{t})*stats.varErrors{t}(k1,k2)/T;
                end
                varBeta{k1,k2} = mean(reshape(cell2mat(tempVarBeta_t),J,J,T),3);
            end
        end 
        varBetaLR = cell2mat(varBeta);

        % compute the test statistic
        st.teststat = T*betaLR*(varBetaLR\betaLR'); % Wald statistic, distributed Chi-squared with M*J degrees of freedom
        
        % do the test
        critical_value = chi2inv(significance,M*J);
        if st.teststat > critical_value
            st.h = 1; %reject the null
        else
            st.h = 0;
        end
        pvalue = 1 - chi2cdf(st.teststat,M*J);
        conf_intvl = NaN; % there is no confidence interval since there is no single point estimate
        % find critical values
        st.critical99 = chi2inv(0.99,M*J);
        st.critical95 = chi2inv(0.95,M*J);
        st.critical90 = chi2inv(0.90,M*J);

    case 'lr alpha'
        point_estimate = mean(alpha);
        % compute variance of alpha(t) and long-run alpha
        varAlpha_t = cell2mat(cellfun(@(varErrors)  diag(varErrors)'./T,stats.varErrors','UniformOutput',false));
        varAlphaLR = mean(varAlpha_t);
        % compute the tstat
        tstatAlphaLR = point_estimate./sqrt(varAlphaLR);
        % do the test
        [st.h,pvalue,~,st.teststat] = ztest(tstatAlphaLR,0,1,1-significance,'both',1);
        % get confidence interval
        interval = norminv([significance/2 1-significance/2]);
        conf_intvl = [point_estimate+interval(1)*sqrt(varAlphaLR);point_estimate+interval(2)*sqrt(varAlphaLR)];
        % put in desired form
        conf_intvl = mat2cell(reshape(conf_intvl,1,2*M),1,repmat(2,1,M));
        % find critical values
        st.critical99 = norminv(1-0.01/2);
        st.critical95 = norminv(1-0.05/2);
        st.critical90 = norminv(1-0.1/2);
    case 'all lr alpha'
        point_estimate = NaN; % there is no point estimate
        % compute variance-covariance of LR alpha
        varAllAlphaLR = mean(reshape(cell2mat(stats.varErrors),M,M,T),3);
        % compute the tstat
        st.teststat = T*mean(alpha)*(varAllAlphaLR\mean(alpha)'); % Wald statistic, distributed Chi-squared with M degrees of freedom
        % do the test
        critical_value = chi2inv(significance,M);
        if st.teststat > critical_value
            st.h = 1; %reject the null
        else
            st.h = 0;
        end
        pvalue = 1 - chi2cdf(st.teststat,M);
        conf_intvl = NaN; % there is no confidence interval since there is no single point estimate
        % find critical values
        st.critical99 = chi2inv(0.99,M);
        st.critical95 = chi2inv(0.95,M);
        st.critical90 = chi2inv(0.90,M);
    case 'constant beta'
        point_estimate = cellfun(@(xx) mean(xx),beta,'UniformOutput',false);% is beta(t) constant and equal to its long run mean?
        varAlpha = cell2mat(cellfun(@(varErrors)  diag(varErrors)',stats.varErrors','UniformOutput',false)); % variance of alpha(t), will be needed later
        varBeta = cell(1,M);
        demeanBeta = cell(1,M);
        Wk = nan(1,M);
        for m = 1:M % loop through left hand side variables
            % compute variance of beta(t)
            varBeta{m} = cellfun(@(varX,varErrors) kappa2*inv_posdef(varX)*(varErrors(m,m)')/(T*stats.bandwidth(m)),stats.varX{m},stats.varErrors','UniformOutput',false);
            % remove means from beta
            demeanBeta{m} = mat2cell(bsxfun(@minus,beta{m},point_estimate{m}),ones(T,1),J);
            % compute the mean squared error, weighted by the std dev
            Wk(m) = mean(cellfun(@(bb,vv) bb*vv*bb',demeanBeta{m},varBeta{m}).*varAlpha(:,m).^(-1));
        end
        % do the test
        mBeta = kappa2*J./(stats.bandwidth*T);
        v2Beta = 2*J*convK2./(stats.bandwidth*T^3);
        st.teststat = (Wk-mBeta)./sqrt(v2Beta); % is distributed N(0,1)       
        [st.h,pvalue,~,st.teststat] = ztest(st.teststat',0,1,1-significance,'both',2);
        % adjust confidence interval of tstat to confidence interval of
        % variable and put in desired shape
        conf_intvl = NaN; % there is no confidence interval
        % find critical values
        st.critical99 = norminv(1-0.01/2);
        st.critical95 = norminv(1-0.05/2);
        st.critical90 = norminv(1-0.1/2);
        
    case {'constant alpha','constant alpha zero'}
        % compute variance of alpha(t)
        varAlpha = cell2mat(cellfun(@(varErrors)  diag(varErrors)',stats.varErrors','UniformOutput',false));
        if strcmpi(test,'constant alpha')
            point_estimate = mean(alpha); % is alpha(t) constant and equal to its long run mean?
            % compute the mean squared error, weighted by the std dev
            Wk = mean(bsxfun(@minus,alpha,point_estimate).^2./varAlpha);
        elseif strcmpi(test,'constant alpha zero')
            point_estimate = zeros(1,M); % is alpha(t) constant and equal to zero?
            % compute the mean squared error, weighted by the std dev
            Wk = mean(bsxfun(@minus,alpha,point_estimate).^2./varAlpha);
        end
        % do the test
        mAlpha = kappa2./(stats.bandwidth*T);
        v2Alpha = 2*convK2./(stats.bandwidth*T^3);
        st.teststat = (Wk-mAlpha)./sqrt(v2Alpha); % is distributed N(0,1)       
        [st.h,pvalue,~,st.teststat] = ztest(st.teststat',0,1,1-significance,'both',2);
        % adjust confidence interval of tstat to confidence interval of
        % variable and put in desired shape
        conf_intvl = NaN; % there is no confidence interval
        % find critical values
        st.critical99 = norminv(1-0.01/2);
        st.critical95 = norminv(1-0.05/2);
        st.critical90 = norminv(1-0.1/2);
    otherwise
        error('hypothesis_test:NoTest',['Test' test 'not defined']);
end

% subfunctions
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

end