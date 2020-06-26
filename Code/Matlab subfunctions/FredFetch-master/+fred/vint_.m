function [vintdata] = vint_(series, vint_date, pseudo, varargin)

  series = upper(series);

  %% Try to grab the data
  [query, success] = fred.ReadFredData_('series_id', series, 'realtime_start', vint_date, 'realtime_end', vint_date, varargin{:});


  %% Return on errors
  if ~success

    % Maybe construct a pseudo vintage
    if pseudo
      vintdata = fred.vintsFromAll_(series, vint_date, 1, varargin{:});
    else
      vintdata.info      = query;
      vintdata.series    = series;
      vintdata.frequency = '';
      vintdata.units     = '';
      vintdata.pseudo    = [];
      vintdata.realtime  = [];
      vintdata.date      = [];
      vintdata.value     = [];
    end

    return
  end

  %% Parse the data
  vintdata.info   = query.info.seriess{:};
  vintdata.series = series;

  frqcy_ind = find(strcmp('frequency', varargin));
  if ~isempty(frqcy_ind)
    vintdata.frequency = upper(varargin{frqcy_ind+1});
  else
    vintdata.frequency = query.info.seriess{end}.frequency_short;
  end

  vintdata.units  = query.obs.units;
  vintdata.pseudo = NaN;

  obs               = vertcat(query.obs.observations{:});
  Nobs              = length(obs);
  vintdata.realtime = fred.dtnum(vint_date,1);
  vintdata.date     = nan(Nobs,1);
  vintdata.value    = nan(Nobs,1);
  for t = 1:Nobs
    vintdata.date(t) = fred.dtnum(obs(t).date,1);
    val = str2num(obs(t).value);
    if ~isempty(val)
      vintdata.value(t) = val;
    end
  end

end
