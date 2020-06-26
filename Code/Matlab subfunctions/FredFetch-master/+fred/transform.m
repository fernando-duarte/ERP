function [data] = transform(data, tf, frequency)


  if ischar(tf), tf = {tf}; end
  if exist('frequency', 'var')
    if ischar(frequency), frequency = {frequency}; end
  end

  %% Handle the different types of data structures that might come in

  % If data is a struct that would come out of calls to fred.latest or
  % fred.vint
  if isstruct(data)

    % If multiple series stacked in an array struct.  Each element of
    % the struct array should have data for the same series in "value"
    % field (whether value is simply a column vector or a matrix of the
    % series at different vintage dates).
    %
    % This if block operates on what you get back from calls to
    % fred.vint with multiple series *and* vintages or calls fred.vint
    % and fred.latest with toDatasetByVint = 0.
    if length(data) > 1

      Nseries = length(data);
      fromUnits = data(n).units;
      toUnits   = tf{n};
      trivial   = strcmp(fromUnits, toUnits);
      for n = 1:Nseries
        if ~trivial && strcmp(fromUnits, 'lin')
          [data(n).value, valid] = fred.transform_(data(n).value, toUnits, data(n).frequency);
          if valid
            data(n).units = toUnits;
          end

        elseif ~trivial
          warning(sprintf('Cannot go from %s to %s units', fromUnits, toUnits))
        end
      end

    % data.value is same series over multiple vintage dates in the cols
    elseif ischar(data.series)
      [data.value, valid] = fred.transform_(data.value, tf{:}, data.frequency);
      if valid
        data.units = tf{:};
      end

    % Dataset merges multiple series, possible frequency-mixing
    elseif iscell(data.series)

      % Loop over series and make the transformation
      Nseries = length(data.series);
      for n = 1:Nseries

        fromUnits = data.units{n};
        toUnits   = tf{n};
        trivial   = strcmp(fromUnits, toUnits);

        % Only transform if going from lin to a pct change or diff
        if ~trivial && strcmp(fromUnits, 'lin')
          notNaN = ~isnan(data.value(:,n));
          [data.value(notNaN,n), valid] = fred.transform_(data.value(notNaN,n), toUnits, data.frequency{n});

          if valid
            data.units{n} = toUnits;
          end
        elseif ~trivial
          warning(sprintf('Cannot go from %s to %s units', fromUnits, toUnits))
        end
      end
    end

  % If data is just an array, not a data struct from one of the programs
  elseif isnumeric(data)

    Nseries = size(data,2);
    for n = 1:Nseries
      data(:,n) = fred.transform_(data(:,n), tf{n}, frequency{n});
    end


  else
    error('I don''t know how to transform that')
  end

end
