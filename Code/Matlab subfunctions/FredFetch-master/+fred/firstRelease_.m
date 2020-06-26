function [rlsdata] = firstRelease_(series, varargin)

  %% Grab the series

    % If 'units' is in the varargin, take that out and handle later bc
    % there are issues pulling that from FRED
    [opt, toPass] = fred.parseVarargin_({'units'}, varargin{:});
    if isnumeric(opt.units) && ~opt.units
      units = 'lin';
    else
      units = opt.units;
      if iscell(units)
        units = units{:};
      end
    end

    % Download
    vintdata = fred.vintall(series, toPass{:});

    % Handle errors
    if isempty(vintdata.value)
      rlsdata.info      = vintdata.info;
      rlsdata.series    = vintdata.series;
      rlsdata.frequency = '';
      rlsdata.units     = '';
      rlsdata.date      = [];
      rlsdata.released  = [];
      rlsdata.value     = [];
      return
    end

    % Transform
    Nvint = length(vintdata.realtime);
    [vintdata.value, tfValid] = fred.transform_(vintdata.value, units, vintdata.info(end).frequency_short);

  %% Replace all the values with only the first releases

    rlsdata.info      = vintdata.info;
    rlsdata.series    = vintdata.series;
    rlsdata.frequency = vintdata.info(end).frequency_short;
    if tfValid
      rlsdata.units = units;
    else
      rlsdata.units = 'lin';
    end
    rlsdata.date      = vintdata.date;

    [Nobs,Nvint] = size(vintdata.value);
    rlsdata.released = nan(Nobs,1);
    rlsdata.value    = nan(Nobs,1);
    rlsdata.latest   = nan(Nobs,1);
    keepRow = zeros(Nobs,1);
    for t = 1:Nobs
      first = find(~isnan(vintdata.value(t,:)),1); % First non-nan entry in row t
      if ~isempty(first)
        rlsdata.value(t)    = vintdata.value(t,first);
        rlsdata.latest(t)   = vintdata.value(t,end);
        rlsdata.released(t) = vintdata.realtime(first);
        keepRow(t)  = 1;
      end
    end
    keepRow = find(keepRow);
    rlsdata.date     = rlsdata.date(keepRow);
    rlsdata.value    = rlsdata.value(keepRow);
    rlsdata.latest   = rlsdata.latest(keepRow);
    rlsdata.released = rlsdata.released(keepRow);

end

