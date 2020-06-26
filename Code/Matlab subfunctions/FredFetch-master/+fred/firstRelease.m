function [returned] = firstRelease(series, varargin)

  [opt, toPass] = fred.parseVarargin_({'parworkers', 'pseudo'}, varargin{:});
  returned = fred.dispatch_(0, opt.parworkers, @fred.firstRelease_, series, toPass{:});

end
