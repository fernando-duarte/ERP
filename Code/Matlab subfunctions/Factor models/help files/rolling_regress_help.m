%% rolling_regress
% Estimates coefficients of rolling least-squares regressions.

%% Syntax
%       alpha = rolling_regress(Y,X)
%       [alpha beta] = rolling_regress(Y,X)
%       [...] = rolling_regress(Y,X,kernel)
%       [...] = rolling_regress(Y,X,kernel,bandwidth)
%% Description
% |alpha = rolling_regress(Y,X)| gives the intercepts of rolling weighted least-squares 
% regressions of the columns of Y on X using weights given by |kernel|.
% |X| is a $T\times J$ matrix of $T$ observations of
% $J$ independent variables. |Y| is a $T\times M$ matrix of $T$ observations 
% of $M$ dependent variables. |alpha| is a $T\times M$ matrix of
% intercepts. Note: |X| should not include a column of ones.
%
% |[alpha beta] = rolling_regress(Y,X)| also returns the a $1 \times M$ cell, each
% containing a $T\times J$ matrix of regression coefficients.
%
% |[...] = rolling_regress(Y,X,kernel)| specifies the kernel. If not
% specified, the default is 'gaussian'. Other options are 'one-sided-uniform',
% 'one-sided exponential' or a used-supplied function handle (the kernel must integrate to one).
%
% |[...] = rolling_regress(Y,X,kernel,bandwidth)| specifies the
% kernel's bandwidth. If |bandwidth| is a number, the same bandwidth is
% used for all columns of |Y|. If |bandwidth| is a $1\times K$ vector, then the
% $k^{th}$ element of the vector is used for the $k^{th}$ column of |Y|.
%
% The estimates of alpha and beta corresponding to the $k^{th}$ column of |Y| at
% time $t$ are given by:
%
% $[\hat{\alpha}_k(t), \hat{\beta}_k(t)] = \arg\max_{\alpha,\beta} \sum_{i=1}^n K_{h_kT}(t_i-t)(Y_{k,i}-\alpha-\beta X_i)^2$
%
% where $K_{h_kT}(z)=\frac{K(z/(h_kT))}{h_kT}$ is the weighting function and $h_kT$ the bandwidth.
% The possible kernels are:
%
% |'gaussian'|  $K(z)=\frac{1}{\sqrt{2\pi}} \mathrm{exp} \left(-\frac{z^2}{2} \right)$
%
% |'one-sided-uniform'| $K(z)= \mathbf{1}_{\{-1\leq z\leq 0 \}}$
%
% |'one-sided exponential'| $K(z)= \mathbf{1}_{\{z\leq 0\}} \mathrm{exp}(z)$
%
% If the kernel is a user-defined function handle, it must satisfy $\int_{-\infty}^\infty K(z)\,dz = 1$
%   References:
%      [1] Ang, A. and Kristensen, D. (2011) Testing Conditional Factor
%      Models
%
%   Copyright 2011 Fernando M. Duarte
%   $Revision: 0.1 $  $Date: 01/30/2012$
%




