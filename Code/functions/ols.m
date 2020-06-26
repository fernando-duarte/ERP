function [varargout] = ols(Y, X, includeMean)
%ols(Y, X, includeMean)
%Simple, no-nonsense OLS.  
%Model: Y=[1 X]b+u (observations in rows), with u~N(0,sig).
%If third argument is false do not append a constant regressor.
%
%Outputs:
%   b_hat, u_hat, sig, R^2, Yhat, var(vec(b_hat))
%
%Michael Abrahams, 2011



    if (nargin < 3)
        includeMean = true;
    end
    
    if includeMean 
        X = [ones(size(X,1), 1) X];
    end
    
    b = (X' * X) \ X' * Y; 
    u = Y - X*b;    
    sig = u' * u / (length(Y) - 1);
    

    
    switch nargout
        case 1
            varargout = {b};
        case 2
            varargout = {b, u};
        case 3
            varargout = {b, u, sig};
        case 4            
            SST = sum(demean(Y).^2);
            SSE = sum(u.^2);
            R2 = 1 - SSE./SST;
            varargout = {b, u, sig, R2};
        case 5
            SST = sum(demean(Y).^2);
            SSE = sum(u.^2);
            R2 = 1 - SSE./SST;
            Yhat = X * b;
            varargout = {b, u, sig, R2, Yhat};
        case 6
            SST = sum(demean(Y).^2);
            SSE = sum(u.^2);
            R2 = 1 - SSE./SST;
            Yhat = X * b;
            ups = X' * X;
            varBHat = kron(sig, inv(ups));
            varargout = {b, u, sig, R2, Yhat, varBHat};
        otherwise
            varargout = {b};
            
    end
    
end