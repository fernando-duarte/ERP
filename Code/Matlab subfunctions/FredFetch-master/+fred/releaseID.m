function [returned] = releaseID(series, varargin)

  % Check for parallel workers
  [opt, toPass] = fred.parseVarargin_({'parworkers'}, varargin{:});

  % Dispatch the call
  returned = fred.dispatch_(0, opt.parworkers, @fred.releaseID_, series);

end
