function [calendar] = releaseCalendar(series, varargin)

  %% Parse options and fill in defaults if need be

    % Get some global options
    optGlobal = fred.GlobalOptions();

    % Parse input
    opt = fred.parseVarargin_({'frequency', 'parworkers', 'realtime_start', 'realtime_end'}, varargin{:});

    % Fill in defaults if none given for start and end dates of calendar
    fill = {'realtime_start', 'realtime_end'};
    for n = 1:length(fill)
      fld = fill{n};

      % If opt.realtime_{start,end} is 0, fill default
      if isnumeric(opt.(fld)) && ~opt.(fld)
        opt.(fld) = optGlobal.(fld);

      % Make sure it's in the right fred format form
      else
        opt.(fld) = fred.dtstr(opt.(fld));
      end
    end
    opt.api         = optGlobal.api;
    opt.max_attempt = optGlobal.max_attempt;


  %% Get the frequency of the series given

    % If no native frequency given for series, fetch them
    if isnumeric(opt.frequency) && ~opt.frequency
      fprintf('Fetching data frequencies...')

      infoURL = @(series) sprintf([...
                'https://api.stlouisfed.org/fred/series?' ...
                'series_id=%s' ...
                '&api_key=%s' ...
                '&file_type=json'],...
                series, ...
                opt.api);
      toDispatch = @(series) fred.ReadFredURL_(infoURL(series), 1, opt.max_attempt);

      fromFred_frequencyInfo = fred.dispatch_(0, opt.parworkers, toDispatch, series);
      unpacked = [fromFred_frequencyInfo.seriess]';
      frequency = cellfun(@(s) s.frequency_short, unpacked, 'un', 0);

      fprintf('done\n')

    % If given frequencies, use them
    else
      frequency = opt.frequency;
    end

    frequency = upper(frequency);

  %% Prune daily series bc you don't need a calendar for that shit

    rem = find(strcmp(frequency, 'D'));
    if length(rem)
      daily          = series(rem);
      series(rem)    = [];
      frequency(rem) = [];
      warning(['Removing daily series: ' sprintf(' %s,', daily{:})])
    end

  %% Download the release ids of all series given

    fprintf('Matching series to Fred release code...')

    % Fetch the release IDs for the series given
    fromFred_releaseIDs = fred.releaseID(series);

    % Extract the release IDs from Fred query results
    seriesReleaseIDs_all = [fromFred_releaseIDs.release_id]';

    % Determine unique IDs that we should download release dates for
    seriesReleaseIDs_unique = unique(seriesReleaseIDs_all);

    fprintf('done\n')

  %% Download the release dates for those releases


    fprintf('Fetching release dates...')
    datesURL = @(id) sprintf([...
            'https://api.stlouisfed.org/fred/release/dates?' ...
            'release_id=%d' ...
            '&api_key=%s' ...
            '&realtime_start=%s' ...
            '&realtime_end=%s' ...
            '&include_release_dates_with_no_data=true' ...
            '&file_type=json'],...
            id, ...
            opt.api, ...
            opt.realtime_start, ...
            opt.realtime_end ...
            );
    toDispatch = @(id) fred.ReadFredURL_(datesURL(id), 1, opt.max_attempt);

    % Note: the "series" field will be the release number. That's just
    % an artifact of the way dispatch_ and multipleSeries_ was set up.
    % Will probably be changed in the future.
    fromFred_releaseDates = ...
      fred.dispatch_(0, opt.parworkers, toDispatch, num2cell(seriesReleaseIDs_unique));

    fprintf('done\n')

  %% Extract release dates and IDs for those dates; collapse to unique

    % Transpose nonempty entries of release_dates field so we can stack
    % things without a dimension mismatch warning from Matlab
    for n = 1:length(fromFred_releaseDates)
      if isempty(fromFred_releaseDates(n).release_dates)
        fromFred_releaseDates(n).release_dates = fromFred_releaseDates(n).release_dates';
      end
    end
    stacked = [fromFred_releaseDates.release_dates]';

    % Extract release dates and releaseIDs
    releaseDates    = cellfun(@(s) fred.dtnum(s.date), stacked);
    releaseDatesIDs = cellfun(@(s) s.release_id, stacked);

    % Sort dates and ids by release dates
    [releaseDates, order] = sort(releaseDates);
    releaseDatesIDs       = releaseDatesIDs(order);

    % Collapse release into unique dates
    [releaseDates_unq, ~, newIndex] = unique(releaseDates);
    releaseDatesIDs_unq       = cell(length(releaseDates_unq),1);
    for n = 1:length(releaseDatesIDs_unq)
      grab = find(newIndex == n);
      releaseDatesIDs_unq{n} = releaseDatesIDs(grab)';
    end


  %% Group into data structure

    Nreleases = length(releaseDates_unq);
    calendar  = struct('release_date',  num2cell(releaseDates_unq), ...
                       'series',        cell(Nreleases,1), ...
                       'release_id',    cell(Nreleases,1), ...
                       'frequency',     cell(Nreleases,1), ...
                       'date',          cell(Nreleases,1));

    % Loop over release dates
    for n = 1:Nreleases

      % The release IDs that were released on day "n"
      releasedIDs = releaseDatesIDs_unq{n};

      % The series that belong to those release IDs
      releasedSeries = find(arrayfun(@(sid) any(sid == releasedIDs), seriesReleaseIDs_all));

      calendar(n).series     = series(releasedSeries);
      calendar(n).release_id = seriesReleaseIDs_all(releasedSeries);
      calendar(n).frequency  = frequency(releasedSeries);

      % Assume that the observation date is the month or quarter before
      % the release date
      [y,m,d,~,~,~] = datevec(calendar(n).release_date);
      dt.M = datenum(y,m,1);
      dt.Q = datenum(y,ceil(m/3),1);
      dt.A = datenum(y,1,1);

      calendar(n).date = nan(length(releasedSeries),1);
      for s = 1:length(calendar(n).frequency)
        calendar(n).date(s) = dt.(calendar(n).frequency{s});
      end

    end

end
