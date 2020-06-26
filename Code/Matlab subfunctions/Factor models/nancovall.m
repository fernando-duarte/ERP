function c = nancovall(x,varargin)
%NANCOVALL Covariance matrix, ignoring NaNs.
%   C = NANCOVALL(X), if X is a vector, returns the sample variance of the
%   values in X, treating NaNs as missing values.  For matrices, where
%   each row is an observation and each column a variable, NANCOVALL(X) is
%   the covariance matrix computed dropping individual NaN's.  NANCOVALL(X,Y),
%   where X and Y are matrices with the same number of elements, is equivalent
%   to NANCOVALL([X(:) Y(:)]). NANCOVALL is different from NANCOV in that
%   it doesn't drop entire rows that contain at least one NaN, it just
%   drops single observations that are NaN's.
%
%   NANCOVALL(X) or NANCOVALL(X,Y) normalizes by (N-1) if N>1, where N is the
%   number of observations after removing missing values.  This makes
%   NANCOVALL(X) the best unbiased estimate of the covariance matrix if the
%   observations are from a normal distribution. For N=1, COV normalizes
%   by N.
%
%   NANCOVALL(X,1) or NANCOVALL(X,Y,1) normalizes by N and produces the second
%   moment matrix of the observations about their mean.  NANCOVALL(X,Y,0) is
%   the same as NANCOVALL(X,Y), and NANCOVALL(X,0) is the same as NANCOVALL(X).
%
%   NANCOVALL(X,2) or NANCOVALL(X,Y,2) normalizes by N^2 and produces the second
%   moment matrix of the observations about their mean. This is useful for
%   some time-series applications like Fama-MacBeth regressions.
%
%   Example:  Generate random data having non-zero covariance between
%             column 4 and the other columns.
%       test = [1 nan 2 3;nan 22 33 nan;44 55 66 nan;77 nan 88 nan]; %some data
%       x(1,2) = NaN; x(2,1) = NaN;x(4,2) = NaN;x(3:4,4) = NaN;      %introduce missing values
%       c = nancovall(x)                                             %compute sample covariance
%       d = nancov(x)                                                %compare c with d (d has all NaN's, c has no NaN's)
%
%   Class support for inputs X,Y:
%      float: double, single
%
%   See also COV, NANCOV
%   $Revision: 0.1 $  $Date: 2012/02/22$
%   Note: draws heavily from the function NANCOV provided by Matlab

if nargin<1
    error(message('stats:nancov:NotEnoughInputs'));
end

% Should we use the mle (divide by N) or unbiased estimate (N-1)?
domle = false;
if numel(varargin)>0
    temp = varargin{end};
    if isequal(temp,0) || isequal(temp,1) || isequal(temp,2)
        domle = temp;
        varargin(end) = [];
    end
end

if numel(varargin)>1
    error(message('stats:nancov:TooManyArgs'));
end

scalarxy = false; % nancov(scalar,scalar) is an ambiguous case
if numel(varargin)>0
    y = varargin{1};
    
    % Two inputs, convert to equivalent single input
    x = x(:);
    y = y(:);
    if length(x)~=length(y)
        error(message('stats:nancov:XYmismatch'));
    end
    scalarxy = isscalar(x) && isscalar(y);
    x = [x y];
elseif ndims(x)>2
    error(message('stats:nancov:InputDim'));
end

if isvector(x) && ~scalarxy
    x = x(:);
end

xnan = isnan(x);
[m,n] = size(x);

if isempty(x);
    if (m==0 && n==0)
        c = NaN(class(x));
    else
        c = NaN(n,class(x));
    end
    return;
end

% Compute variance using data ignoring NaN's
x(xnan) = 0;
colsize = sum(~xnan,1);
xmean = sum(x,1) ./ max(1,colsize);
xmean(colsize==0) = NaN;
xc = x - repmat(xmean,m,1);
xc(xnan) = 0;
c = nan(n,class(x));
for indn = 1:n
    for indm = 1:indn
        nobs = sum(~(xnan(:,indn) | xnan(:,indm)));
        if domle == 1
            denom = nobs;
        elseif domle == 2
            denom = nobs^2;
        else
            denom = max(1,nobs-1);
        end
        if nobs ~= 0
            c(indn,indm) = xc(:,indn)'*xc(:,indm)/denom;
        end
    end
end
%c(isnan(c)) = 0;
c= tril(c) + tril(c,-1)';

end

