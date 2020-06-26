function [outdate] = dtnum(indate, fred_fmt)
% FRED.DTNUM datenum for arguments of type char or cell
%  OUTDATE = FRED.DTNUM(INDATE) returns datenum(s) for each date in INDATE.

  if exist('fred_fmt', 'var') && fred_fmt
    fmt = {'yyyy-mm-dd'};
  else
    fmt = {};
  end

  if ischar(indate)
    indate = {indate};
  end
  if iscell(indate)
    outdate = cellfun(@(in) datenum(in, fmt{:}), indate);
  else
    outdate = indate;
  end

end
