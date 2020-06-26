function [opt] = GlobalOptions()

  opt.max_attempt = 3;
  opt.realtime_start = '1776-07-04'; % Very first available
  opt.realtime_end   = '9999-12-31'; % Very last available

  % Whether leading and trailing NaNs on a data matrix should be trimmed
  % when possible
  opt.trimLeadTrailNaN = 1;


  try
    api = textread(['api.txt'], '%s');
    opt.api = api{:};
  catch
    error('Please supply an API Key in +fred/GlobalOptions.m')
  end

end
