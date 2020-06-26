function [data] = multipleSeries_(toDatasetByVint, parworkers, dl_fcn, series, varargin)


  %% Maybe open up parallel pool
  if parworkers
    poolcheck = gcp('nocreate');
    if isempty(poolcheck)
      poolobj  = parpool(min([parworkers length(series)]));
      poolkill = 1;
    else
      poolobj  = poolcheck;
      poolkill = 0;
    end
  end


  %% Download the series individually using dl_fcn

    Nseries    = length(series);
    individual = cell(Nseries, 1);

    % Set up indexing of extra arguments, and if some are cells, we want
    % to loop over them in the calls to dl_fcn
    cells = cellfun(@iscell, varargin);
    toPass = cell(Nseries, length(varargin));
    toPass(:,~cells) = repmat(varargin(~cells), Nseries, 1);
    toPass(:,cells) = [varargin{cells}];

    if parworkers
      parfor s = 1:Nseries
        individual{s} = feval(dl_fcn, series{s}, toPass{s,:});
      end
    else
      for s = 1:Nseries
        individual{s} = feval(dl_fcn, series{s}, toPass{s,:});
      end
    end


  %% Merge or Join the downloaded objects

    % Just stack the returned arrays
    if ~toDatasetByVint
      data = vertcat(individual{:});
      for n = 1:length(series)
        data(n).series = upper(series{n});
      end

    % Reshape into array structure where each element is a different
    % vintage data and the data for different series are merged into a
    % single data matrix of possibly mixed frequency
    else

      fprintf('Reshaping by vintage...\n')
      data = fred.reshapeByVint(individual);
    end


  %% Shut down parpool
  if parworkers && poolkill
    delete(poolobj);
  end

end
