function [data] = latest(series, varargin)

  opt = fred.parseVarargin_({'parworkers', 'pseudo', 'units'}, varargin{:});
  data = fred.dispatch_(opt.toDatasetByVint, opt.parworkers, @fred.latest_, series);

  % Maybe transform
  if iscell(opt.units) | ischar(opt.units)
    data = fred.transform(data, opt.units);
  end

end
