function [varargout] = dtfld(dt, fld)

  switch fld
  case {'year', 'month'}
    [ret.year, ret.month, ret.day] = datevec(dt);

    varargout{1} = ret.(fld);
  case {'quarter'}
    q = datestr(dt, 'QQ');
    varargout{1} = str2num(q(:,2));

  case {'week', 'day'}
    firstInYear = datenum(fred.dtfld(dt, 'year'), 1, 1);
    varargout{1} = str2num(q(:,2));
  end

end
