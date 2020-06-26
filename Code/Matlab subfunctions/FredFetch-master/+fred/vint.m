%% vint.m
%
% Get a dataset for some vintage date. Allows users to query single and
% multiple series.
%
% Input Arguments
% - series      A single series name like 'GDPC1' or a cell of series
%               names like {'GDPC1'; 'NAPM';}
% - vint_date   A vintage date; can be a datenum or datestr (that can be
%               converted by function datenum() to a datenum
% - varargin    Any FRED API options you want to add to the URL,
%               specified in order (option_name, option_value).
%
%               Example:
%                 vint('GDPC1', '2015-01-01', 'observation_start', '1991-01-01', 'units', 'chg')
%
%               This will append '&observation_start=1991-01-01&units=chg'
%               onto the URL.
%
%               NOTE: This library is designed so that you don't
%               have to worry about specifying these extra options, but
%               if you've read the FRED API documentation and you want
%               to make heaver use of it, you have that option.
%
function [vintdata] = vint(series, vint_date, varargin)

  [opt, toPass] = fred.parseVarargin_({'pseudo', 'parworkers'}, varargin{:});

  % Make vint dates a datenum and sort
  vint_date = sort(fred.dtnum(vint_date));
  if any(vint_date < datenum(1991,1,1)) && ~opt.pseudo
    warning('Early vintage date; data might not exist for some or all series.')
  end

  % Call to different functions depending upon whether one or multiple
  % vint dates are specified
  if length(vint_date) > 1
    vintdata = fred.dispatch_(opt.toDatasetByVint, opt.parworkers, @fred.vintsFromAll_, series, vint_date, opt.pseudo, toPass{:});
  else
    vintdata = fred.dispatch_(opt.toDatasetByVint, opt.parworkers, @fred.vint_, series, vint_date, opt.pseudo, toPass{:});
  end


end
