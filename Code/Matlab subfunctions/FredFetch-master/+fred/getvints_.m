function [returned] = getvints_(series, max_attempt)

  %% Options
  opt = fred.GlobalOptions();
  if ~exist('max_attempt', 'var')
    max_attempt = opt.max_attempt;
  end

  %% Construct a query string to pull the release dates
  url = sprintf([...
    'https://api.stlouisfed.org/fred/series/vintagedates?'...
    'series_id=%s' ...
    '&api_key=%s' ...
    '&file_type=json'],...
    series,...
    opt.api);

  try
    vintdates = jsonlab.loadjson(urlread(url));
    returned.series    = series;
    returned.vintdates = fred.dtnum(vintdates.vintage_dates,1)';
    returned.success   = 1;
  catch
    if max_attempt - 1
      pause(5);
      returned = fred.getvints_(series, max_attempt-1);
      return
    else
      returned.series    = series;
      returned.vintdates = [];
      returned.success   = 0;
    end
  end

end
