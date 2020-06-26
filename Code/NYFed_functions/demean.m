% function standardize.m standardizes a data matrix such that each column
% has zero mean and unit variance

function [X_demeaned] = demean(X_nonstd);

[rows, cols] = size(X_nonstd);
X_demeaned = X_nonstd - kron(nanmean(X_nonstd),ones(rows,1));

