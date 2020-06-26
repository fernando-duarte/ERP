%% ReadFredURL.m
%
% This function is the most general "Get/Read something from FRED"
% function. Most other functions that get data from FRED are actually
% just tailored wrappers of this guy.
%
% Behavior:
% ---------
% - Accepts a URL
% - Try at most max_attempt times to download that URL
% - Return either the error message or the URL contents (which are json)
%   in a struct
%
function [returned, success] = ReadFredURL_(url, json, max_attempt)

  %% Get max attempt if not set
  if ~exist('max_attempt', 'var')
    opt = fred.GlobalOptions();
    max_attempt = opt.max_attempt;
  end

  %% Try max_attempt times to download; if error, return error
  try
    if json
      returned = jsonlab.loadjson(urlread(url));
    else
      returned = urlread(url);
    end
    success  = 1;
  catch
    if max_attempt - 1
      [returned, success] = fred.ReadFredURL_(url, json, max_attempt-1);
      return
    else
      returned = lasterror();
      returned.url = url;
      success  = 0;
    end
  end

end
