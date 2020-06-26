function [dataByVint] = reshapeByVint(individualBySeries)
% reshapeByVint - Take a cell array where each element is a structure
%                 with a different *series* (over possibly different
%                 vintage dates along columns of the "value" field) and
%                 reshape it into a structure array with different
%                 *vintages* in each element of the array and different
%                 series along the columns of the "value" field.
%
% Effectively, this puts merged, possibly-mixed frequency *vintage
% datasets* in each element of the returned structure array, which is
% maybe more useful than one series per element.

  opt = fred.GlobalOptions();

  % Make sure the individualbyseries argument is an array structure; if
  % it's a cell of structs, stack them
  if iscell(individualBySeries)
    dataBySeries = vertcat(individualBySeries{:});
  else
    dataBySeries = individualBySeries;
  end
  Nseries = length(dataBySeries);

  % Vintages
  vints = unique([dataBySeries.realtime]);
  Nvint = length(vints);

  % Loop over vintages and build up the dataByVint array struct
  for v = 1:Nvint

    % Current vintage date
    vint = vints(v);
    toStore.realtime = vint;
    fprintf('Creating vintage dataset for %s...\n', datestr(vint));

    % Observation dates
    toStore.date = sort(unique(vertcat(dataBySeries.date)));

    % Add info about the series
    carry_over = {'info', 'series', 'frequency', 'units'};
    for n = 1:length(carry_over)
      toStore.(carry_over{n}) = {dataBySeries.(carry_over{n})}';
    end

    % Fill in pseudo and data across different series for this vintage
    toStore.pseudo = nan(Nseries,1);
    toStore.value  = nan(length(toStore.date), Nseries);
    for n = 1:Nseries

      % For series n, find the index among its different vintage dates
      % that matches vint
      vintMatchInd = find(vint == dataBySeries(n).realtime);

      if ~isempty(vintMatchInd)
        % Store pseudo info
        toStore.pseudo(n) = dataBySeries(n).pseudo(vintMatchInd);

        % Fill in the data, and data, accounting for mixed frequency
        insert = arrayfun(@(t) find(toStore.date==t), dataBySeries(n).date);
        toStore.value(insert,n) = dataBySeries(n).value(:,vintMatchInd);
      end
    end

    % Trim the leading and trailing nans
    if opt.trimLeadTrailNaN
      [toStore.value, rem] = fred.RemLeadTrailNaN_(toStore.value, 'all');
      toStore.date = toStore.date(~rem);
    end

    toStore = orderfields(toStore, [carry_over, {'pseudo', 'realtime', 'date', 'value'}]);
    dataByVint(v) = toStore;
  end

end
