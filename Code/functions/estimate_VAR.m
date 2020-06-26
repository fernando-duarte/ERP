% Function estimate_VAR delivers OLS estimates of VAR coefficients and a
% matrix of errors
%
% Syntax
% [mu, Phi, Omega, u] = estimate_VAR(Z, nlags, const)
%
% copyright: Emanuel Moench, Humboldt University Berlin, January 2005
% Email : moench@wiwi.hu-berlin.de

function [mu, Phi, Omega, u] = estimate_VAR(Z, nlags, const)

[T, k] = size(Z);
p = nlags;

mu = zeros(k,1);
Phi = zeros(k,k,p);
u = zeros(T-p,k);

for j = 1:k,
    z_j_t = Z(p+1:end,j);
    z_j_lags = [];
    lag = 1;    while lag <= p,
        if lag == 1, 
            z_j_lags = Z(1:end-lag,:); 
        else 
            z_j_lags = [z_j_lags(2:end,:) Z(1:end-lag,:)];
        end
        lag = lag + 1;
                end;
    if const,
        z_j_lags = [ones(T-p,1) z_j_lags];
    end;   
    phi_all = (inv(z_j_lags'*z_j_lags)*z_j_lags'*z_j_t)';       
    if const,
        mu(j) = phi_all(1);
    end;
    lag = 1; while lag <= p,
        Phi(j,:,lag) = phi_all(const+(lag-1)*k+1:const+lag*k);
        lag = lag + 1;
    end;
    u(:,j) = z_j_t - z_j_lags*phi_all';
end;
Omega = 1/(T-1)*u'*u;


