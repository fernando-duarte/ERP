function [vintdata] = dispatch_(toDatasetByVint, parworkers, fcn, series, varargin)


  if iscell(series) && (length(series) == 1)
    series = series{1};
  end
  if isstr(series)
    vintdata = feval(fcn, series, varargin{:});
  else
    vintdata = fred.multipleSeries_(toDatasetByVint, parworkers, fcn, series, varargin{:});
  end


end
