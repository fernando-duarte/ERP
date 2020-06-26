function [X, valid] = transform_(X, tf, frequency)
% TRANSFORM_ - Transform a vector or matrix from levels into specified
%              other specified unit
%
% [X, valid] = transform_(X, TF, FREQUENCY) transforms vector X (or all
% columns of matrix X) into units specified by tf. Assumes matrix X is
% all the same frequency of data (no irregular spacing, like you might
% get from fred.vint or fred.latest).
%
% To transform a matrix with different frequencies of data or to apply
% different transformations to different columns, see fred.transform (no
% underscore), which is a more general wrapper for this function.


  %% Define transformation types and check chosen transformation

    % Define transformations that compute changes from the previous
    % period, versus changes from the previous year
    chgPrevPd = {'chg', 'pch', 'cch', 'pca', 'cca'};
    chgPrevYr = {'ch1', 'pc1'};

    % Define series where we want to scale to annual
    scaleAnnual = {'pca', 'cca'};

    % If the transformation type is not valid, return
    if strcmp('lin', tf)
      valid = 1;
      return
    elseif ~any(strcmp([chgPrevPd, chgPrevYr], tf))
      valid = 0;
      warning(sprintf('%s is not a valid transformation type', tf));
      return
    else
      valid = 1;
    end

  %% Set periods per year based on the frequency

    switch upper(frequency)
      case 'Q'
        pdsPerYear = 4;
      case 'M'
        pdsPerYear = 12;
      case 'W'
        pdsPerYear = 52;
      case 'D'
        pdsPerYear = 365;
      otherwise
        pdsPerYear = [];
    end
    if isempty(pdsPerYear)
      error(['Periods per year not defined for frquency of type ' frequency]);
    end

  %% Set Nlags based on tf's chgPrev group

  switch tf
  case chgPrevPd
    Nlags = 1;
  case chgPrevYr
    Nlags = pdsPerYear;
  end

  %% Set scaling

    scaling = 1;
    if any(strcmp(tf, scaleAnnual))
      scaling = pdsPerYear;
    end

  %% Make lead and lag mats, which we will operate on

    Xlead = X(Nlags+1:end,:);
    Xlag  = X(1:end-Nlags,:);
    pad   = nan(Nlags,size(X,2));

  %% Compute
  switch tf
    case {'chg', 'ch1'}
      X = [pad; Xlead - Xlag];

    case {'pch', 'pc1','pca'}
      X = [pad; 100*(((Xlead ./ Xlag).^scaling)-1)];

    case {'cch', 'cca'}
      X = [pad; scaling*100*(log(Xlead)-log(Xlag))];

    case 'log'
      X = log(X);

  end

end
