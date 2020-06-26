function [ release_dates ] = ...
  query_releasedates(api, series, frequency, start_date, end_date, varargin)

  if ~iscell(series), series = {series}; end


  %% Construct a query string and to pull the release dates
  url_start = 'http://api.stlouisfed.org/fred/series/observations?';
  url = @(srs) ...
         sprintf('%sseries_id=%s&api_key=%s&realtime_start=%s&realtime_end=%s&file_type=json&output_type=4', ...
                 url_start, srs, api, start_date, end_date);


  %% Loop over series and pull
  for s = 1:length(series)
    from_fred = loadjson(urlread(url(series{s})));

    obs_dates     = cellfun(@(obs) datenum(obs.date), from_fred.observations)';
    obs_dates     = arrayfun(@(dt) addtodate(dt, 2, 'month'), obs_dates); % Add 2 months so the obs dates is the first day of the last month in the quarter
    release_dates = cellfun(@(obs) datenum(obs.realtime_start), from_fred.observations)';

    if nargin > 5
      saving = varargin{1};
      f = fopen([saving series{s} '.txt'], 'w');
      fprintf(f, 'obs_date\trelease_date\n');
      arrayfun(@(row) fprintf(f, '%9.0f\t%9.0f\n', obs_dates(row), release_dates(row)), 1:length(obs_dates));
      fclose(f);
    end
  end

end
