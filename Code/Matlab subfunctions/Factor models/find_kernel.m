function [y kappa2 convK2 one_sided_flag]= find_kernel(kernel)
% FIND_KERNEL returns the kernel function described by KERNEL 
%
% Y = FIND_KERNEL(kernel) returns the function handle Y described the the
% string or function handle KERNEL. If KERNEL is a string, FIND_KERNEL 
% finds the corresponding kernel function. If KERNEL is a function handle and 
% the function integrates to one over the real line, then FIND_KERNEL
% returns the function handle KERNEL. If KERNEL is empty, the gaussian
% kernel is returned
%
% [Y KAPPA2] = FIND_KERNEL(kernel) returns the "kernel variance", which is
% the integral over the entire real line of the kernel squared
% [Y KAPPA2 CONVK2] = FIND_KERNEL(kernel) returns the integral over the real 
% line of the squared convolution of the kernel function with itself
% [Y KAPPA2 CONVK2 ONE_SIDED_FLAG] = FIND_KERNEL(kernel) returns true for
% ONE_SIDED_FLAG if the kernel is a one-sided kernel and false otherwise

grid = linspace(-1e2,1e2,1e4); % grid to integrate kernel over the real line

one_sided_flag = false;
if isempty(kernel)
    y = @(x) sqrt(2*pi)^(-1)*exp(-x.^2/2);
elseif strcmpi(kernel,'gaussian')
    y = @(x) sqrt(2*pi)^(-1)*exp(-x.^2/2);
elseif strcmpi(kernel,'one_sided_uniform')
    y = @(x) double(x>=-1 & x<=0);
    one_sided_flag = true;
elseif strcmpi(kernel,'two_sided_uniform')
    y = @(x) double(x>=-1 & x<=1);
elseif strcmpi(kernel,'one_sided_exponential')
    y = @(x) double(x<=0).*min(exp(x),1e50);
    one_sided_flag = true;
elseif strcmpi(kernel,'lag_one_sided_uniform')
    y = @(x) double(x>=-1 & x<0);
    one_sided_flag = true;
elseif strcmpi(kernel,'lag_one_sided_exponential')
    y = @(x) double(x<0).*min(exp(x),1e50);
    one_sided_flag = true;
elseif ~isa(kernel, 'function_handle')
    error('rolling_regress:BadKernel','User-defined kernel must be a function handle');
else 
    y = kernel;
    
    % check if it's one-sided
    one_sided_flag = y(eps)==0;
    
    %check that user supplied kernel adds up to one
    if abs(y(grid(2:end))*diff(grid')-1)>1e-3
        error('rolling_regress:BadKernel','User-defined kernel must integrate to one');
    end
end

% compute kernel "variance"
kappa2 = y(grid(2:end)).^2*diff(grid');
if nargout >= 3
    convK2 = conv(y(grid(2:end)),y(grid(2:end)),'same').^2*diff(grid').^3;
end