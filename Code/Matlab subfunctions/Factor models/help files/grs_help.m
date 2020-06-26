%% grs_test
% Gibbons-Ross-Shanken test of zero pricing errors in a factor model

%% Syntax
%       grs = grs_test(returns,factors)
%       [grs alpha] = grs_test(returns,factors)
%       [grs alpha beta] = grs_test(returns,factors)
%       [grs alpha beta stats] = grs_test(returns,factors)
%       [...] = grs_test(returns,factors,method)
%       [...] = grs_test(returns,factors,'GMM',lags)

%% Model
% grs runs the linear time series regression
%%
% $R_t = \alpha + \beta F_t + \varepsilon_t$
%%
% to test the hypothesis that all elements of the vector $\alpha$ are zero. $R_t$ is the $N\times 1$ vector of excess
% stock returns at time $t$ and $F_t$ is the $K\times 1$ vector of factor returns
% (assumed to be excess returns of zero-cost portfolios). grs reports the
% statistic
%
% $$grs = \frac{T-N-K}{N} (1 + \hat{\mu}_F' \hat{\Sigma}_F^{-1}
% \hat{\mu}_F)^{-1}  \hat{\alpha}' \hat{\Sigma}_\varepsilon^{-1}
% \hat{\alpha}$$

%%
% and uses its finite sample distribution $F(N,T-N-K)$ to produce p-values
% and other statistics. $T$ is the total number of observations and
%
% $$\hat{\mu} = \frac{1}{T} \sum\limits_{t=1}^T R_t$
%
% $$\hat{\mu}_F = \frac{1}{T} \sum\limits_{t=1}^T F_t$
%
% $$\hat{\Sigma}_F = \frac{1}{T-K} \sum\limits_{t=1}^T (F_t-\hat{\mu}_F)(F_t-\hat{\mu}_F)'$
%
% $$\hat{\Sigma}_\varepsilon = \frac{1}{T} \sum\limits_{t=1}^T
% (R_t-\hat{\alpha}-\hat{\beta}F_t)(R_t-\hat{\alpha}-\hat{\beta}F_t)'$.

%% Description
% |grs = grs_test(returns,factors)| gives the grs statistic (a number) for
% the null hypothesis that that all elements of the vector $\alpha$ are
% jointly zero. |Factors| is a $T\times K$ matrix of $T$ observations of
% $K$ pricing factors. |returns| is a $T\times N$ matrix of $T$ observations 
% of the returns of $N$ assets. Note: the matrix of factors should not
% include a column of ones.
%
% |[grs alpha] = grs_test(returns,factors)| gives the $1\times N$
% vector of intercepts (or pricing errors) $\alpha$.
%
% |[grs alpha beta] = grs_test(returns,factors)| gives the $K\times N$
% matrix of regression coefficients $\beta$.
%
% |[grs alpha beta stats] = grs_test(returns,factors)| gives a structure
% with the following statistics
%%
% |stats.R2| mean $R^2$ of all the the regressions
%
% |stats.adjR2| mean adjusted $R^2$ of all the regressions
%
% |stats.absDev| average absolute intercept $\frac{1}{N}
% \sum\limits_{i=1}^N |\alpha_i|$
%
% |stats.stdErr| $1\times N$ vector of standard errors of $\alpha$
%
% |stats.tstat| $1\times N$ vector of t-statisitics for $\alpha$
%
% |stats.SR| "Sharpe Ratio" of the intercepts $(\hat{\alpha}' \hat{\Sigma}_\varepsilon^{-1}
% \hat{\alpha})^{\frac{1}{2}}$. Lewellen, Nagel and Shanken (2010) recommend
% reporting it
%
% |stats.p| p-value of the grs test
%
% |stats.sigmaF| $K\times K$ matrix of covariances of $F_t$
%
% |stats.sigmaEpsilon| $N\times N$ matrix of covariances of $\varepsilon_t$
%
% |stats.sharpeF| Sharpe ratio ($\hat{\mu}_F' \hat{\Sigma}_F^{-1}
% \hat{\mu}_F)^{\frac{1}{2}}$ of the ex-post tangency portfolio constructed
% from the $K$ factors only
%
% |stats.sharpeAll| Sharpe ratio of the ex-post tangency portfolio
% constructed using the $K$ factors and the $N$ portfolios
%
% |stats.criticalValue| $1$, $5$ and $10\%$ critical values for the grs test
%
% |[...] = grs_test(returns,factors,method)| uses the string _method_ to select alternative
% ways to estimate the model and test $\alpha=0$
% 
% |'GRS'| This is the default. The computations are those described above.
%
% |'Asymptotic GRS'| Does not have a finite sample adjustment, ie. the factor
% $\frac{T-N-K}{N}$ in the grs statistic is replaced by $T$ and the limiting
% distribution is $\chi^2_N$ instead of  $F(N,T-N-K)$.
%
% |'GMM'| Uses a GMM procedure with Newey-West corrected standard
% errors:
%
% $$grs = T\hat{\alpha}' var(\hat{\alpha})^{-1} \hat{\alpha} \sim \chi^2_N$$
%
% with $var(\hat{\alpha})$ computed by adjusting for heteroskedasticity and
% serial correlation using the estimator of Newey and West (1987). You can 
% specify the number of lags to use by using the input
% argument |lags|. The default is $floor(T^{\frac{1}{4}})$.



