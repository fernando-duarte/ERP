function [vintdata] = vintall(series, varargin)

  [opt, toPass] = fred.parseVarargin_({'parworkers', 'pseudo'}, varargin{:});

  % Dispatch the call to different function depending upon whether one
  % or multiple series are specified
  very_first = '1776-07-04';
  very_last  = '9999-12-31';

  vintdata = fred.dispatch_(0, opt.parworkers, @fred.vintrange, series, very_first, very_last, toPass{:});

  % Add publication delays
  if ischar(series)
    series = {series};
  end
  for s = 1:length(series)
    if ~isempty(vintdata(s).value)
      vintdata(s).publag = fred.computePublag_(vintdata(s).date, vintdata(s).realtime, vintdata(s).value);
    else
      vintdata(s).publag = repmat(struct(), 0, 1);
    end
  end

end
