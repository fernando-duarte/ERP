function [returned] = aggregate_(bob, X, dt, native, intended, method)


  %% Ordering

    % Define ordering of series from highest to lowest frqcy
    ordering  = {'d', 'w', 'm', 'q', 'a'};

    % Make sure we're going from higher to lower frequency, not the
    % other way around
    if find(strcmp(ordering,intended)) < find(strcmp(ordering, native))
      error(sprintf('Cannot disaggregate series from %s to %s frequency', native, intended))
    end

  %% Define thresholds

    % Define how many of the native frequency observations we need to
    % get to the intended frequency when at edge cases (beginning or end
    % of date vector)
    switch native
      case 'd'
        thresholds = [NaN 4 17 52 208];
      case 'w'
        thresholds = [NaN NaN 3 10 39];
      case 'm'
        thresholds = [NaN NaN NaN 3 12];
      case 'q'
        thresholds = [NaN NaN NaN NaN 4];
    end
    threshold = thresholds(strcmp(intended, ordering));


  %% Define aggregation method, whether we need worry about threshold

    if ~exist('method', 'var')
      method = 'avg';
    end
    switch method
      case 'avg'
        agg = @(x) nanmean(x);
      case 'sum'
        agg = @(x) nansum(x);
      case 'eop'
        agg = @(x) x( find(~isnan(x),1,'last') );
    end

  %% Loop over columns of X and aggregate

    % Take the given native frequency dates, and compute the
    % corresponding unique intended-frequency dates
    allIntended = fred.dtGivenFrequency_(intended, dt(1), dt(end));
    Ntagg = length(allIntended);

    % Loop over columns and do one at a time
    Y = nan(Ntagg,size(X,2));
    for n = 1:size(X,2)

      % Index of last obs for this column
      lastRow = find(~isnan(X(:,n)), 1, 'last');

      % Loop over lower, intended frequency observation dates
      for t = 2:Ntagg

        % For lower, intended frequency time period t, define the
        % starting and stopping dates of the corresponding native higher
        % frequency time periods
        start = allIntended(t-1);
        stop  = allIntended(t);

        % Get the rows that are within the window
        rowMatch = find((dt > start) & (dt <= stop));

        % Check if this is an edge case for the column
        edgeCase = ( (t==2) | any(rowMatch >= lastRow) | isempty(rowMatch));

        % If not an edge case or it passes the threshold, compute and store
        if ~edgeCase | (edgeCase & (length(rowMatch) >= threshold))
          Y(t,n) = agg(X(rowMatch,n));
        end
      end
    end

    %% Method
    % You have the native frequency dates
    % Take the first and last date
    % Construct all the date ranges for that frequency
    % Loop over those possible dates and see if you can fill them in
    returned.value = Y;
    returned.date = allIntended;

    % Don't know what how the fuck fred is aggregating shit
    % Need to figure this out. Maybe start by verifying how
    %
    % Ignore daily for now
    % Weekly to monthly, quarterly, annual


end
