function [url] = MakeFredURL_(urltype, varargin)

  %% Get global options
  opt = fred.GlobalOptions();

  %% Start of URL
  %  This is specific to the type of pull you want to
  %  make (data observations, series url, etc.)
  switch urltype
    % Information about a series
    case 'info'
      url = 'https://api.stlouisfed.org/fred/series?';
    case 'obs'
      url = 'https://api.stlouisfed.org/fred/series/observations?';
    otherwise
      error(sprintf('No url specificed to fetch data of type ''%s''', urltype));
  end


  %% Add api and json file type
  url = sprintf('%sapi_key=%s&file_type=json', url, opt.api);


  %% Walk through and append additional options passed
  date_flds = {'observation_end'; 'observation_start'; 'realtime_end'; 'realtime_start'; 'vintage_dates'};
  Nextra    = length(varargin);
  for a = 1:2:Nextra

    % If a datefield is a datenum, convert to Fred format
    if any(strcmp(date_flds, varargin{a})) && isnumeric(varargin{a+1})
      varargin{a+1} = fred.dtstr(varargin{a+1});
    end

    url = sprintf('%s&%s=%s', url, varargin{a}, varargin{a+1});
  end

end
