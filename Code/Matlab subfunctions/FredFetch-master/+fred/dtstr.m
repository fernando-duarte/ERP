function [outdate] = dtstr(indate, forceCell)

  if ~exist('forceCell', 'var')
    forceCell = 0;
  end

  fred_fmt = 'yyyy-mm-dd';

  if isnumeric(indate)
    outdate = arrayfun(@(in) datestr(in, fred_fmt), indate, 'un', 0);
    if length(outdate) == 1
      outdate = outdate{1};
    end
  else
    outdate = datestr(indate, fred_fmt);
  end

  if ischar(outdate) && forceCell
    outdate = {outdate};
  end

end
