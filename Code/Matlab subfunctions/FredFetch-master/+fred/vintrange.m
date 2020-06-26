function [vintdata] = vintrange(series, realtime_start, realtime_end, varargin)

  %% Handle units and frequency changes
  %
  % For whatever reason, fred freaks out if you try to agggreate to a
  % lower frequency or change the units. So we want strip those out of
  % the optional args we pass to the fred api and handle them ourselves
  % below
  [opt, toPass] = fred.parseVarargin_({'frequency', 'units'}, varargin{:});

  series = upper(series);
  realtime_start = fred.dtstr(realtime_start);
  realtime_end   = fred.dtstr(realtime_end);

  %% Try to grab the data
  fprintf('Downloading %s...\n', series);
  [query, success] = fred.ReadFredData_('series_id', series, 'realtime_start', realtime_start, 'realtime_end', realtime_end, toPass{:});

  %% Return on errors
  if ~success
    vintdata.info      = query;
    vintdata.series    = series;
    vintdata.frequency = '';
    vintdata.units     = '';
    vintdata.pseudo    = [];
    vintdata.realtime  = [];
    vintdata.date      = [];
    vintdata.value     = [];
    return
  end

  %% Parse the data
  vintdata.info      = vertcat(query.info.seriess{:});
  vintdata.series    = series;
  vintdata.frequency = query.info.seriess{end}.frequency_short;
  vintdata.units     = query.obs.units;
  obs                = vertcat(query.obs.observations{:});
  Nall               = length(obs);

  %% Put the dates and values into arrays
  flds = {'realtime_start', 'realtime_end', 'date', 'value'};
  isdt = [1 1 1 0];
  for n = 1:length(flds)
    if isdt(n)
      all.(flds{n}) = fred.dtnum({obs(:).(flds{n})}',1);
      unq.(flds{n}) = unique(all.(flds{n}));
    else
      all.(flds{n}) = arrayfun(@(t) str2double(obs(t).(flds{n})), 1:Nall)';
    end
  end

  %% Define a field for all realtime dates

  realtime = unique([unq.realtime_start; unq.realtime_end]);
    % Note: Used to fill the data matrix by looping over observations
    % (looping over rows) and filling the data from realtime start to
    % realtime start (filling left to right, overwriting vintage dates)
    % in the vintage data matrix.
    %
    % This was maybe problematic. In particular, if there are two
    % realtime starts (01-Jan-2000 and 01-Jan-2010) with two separate
    % values, 10 and 20, we used to put a value of 10 from 01-Jan-2000
    % until 31-Dec-2009. Then 20 from there. That's maybe wrong if
    % there's a realtime_end for value 10 of, say, 01-Jan-2005.
    %
    % This could happen if a data provider decides to drop old
    % observations because maybe they were wrong for some reason. To be
    % fair, I have no concrete evidence this has occurred (except some
    % vague recollection of seeing something like this at some point).
    % Regardless, the new method implemented below is super safe and
    % guaranteed to be correct, while the old method was not.

  %% Loop over realtime dates and fill in what you have at that time

    % Set up vintdata structure that will hold everything
    vintdata.pseudo   = nan(length(realtime),1);
    vintdata.realtime = realtime;
    vintdata.date     = unq.date;
    Nobs  = length(vintdata.date);
    Nvint = length(vintdata.realtime);
    vintdata.value = nan(Nobs, Nvint);

    % Loop over vintdates
    date2row   = arrayfun(@(dt) find(dt == vintdata.date), all.date); % Matches all dates to their row index in vintdata.date
    deleteCol  = zeros(1,Nvint);
    match_last = [];
    for n = 1:Nvint

      % Indices within "all" that are relevant for this vintage date
      match = find((realtime(n) >= all.realtime_start) .* (realtime(n) <= all.realtime_end));

      % If match for this vintage data is exactly the same as the last,
      % mark column for delection and move on
      if length(match) == length(match_last) && ~any(match ~= match_last)
        deleteCol(n) = 1;

      else % If genuinely different, store the data
        fillrows = date2row(match);
        vintdata.value(fillrows,n) = all.value(match);
        match_last = match;
      end
    end

    % Kill the columns that are identical to others
    vintdata.realtime(find(deleteCol)) = [];
    vintdata.value(:,find(deleteCol))  = [];
    vintdata.pseudo(find(deleteCol))   = [];


    % Transform the data if the user wanted to
    if opt.units
      vintdata = fred.transform(vintdata, opt.units);
    end

end
