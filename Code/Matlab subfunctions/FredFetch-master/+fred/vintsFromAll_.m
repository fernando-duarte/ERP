function [returned] = vintsFromAll_(series, vintdates, pseudo, varargin)

  % Convert all vintage dates to datenums
  vintdates = fred.dtnum(vintdates,1);

  % If pseudo specified and you don't actually need a pseudo vintage (bc
  % a vintage is available), turn it off
  if pseudo
    available = fred.getvints(series);
    if available.success && vintdates(1) >= available.vintdates(1)
      pseudo = 0;
    end
  end

  %% Pull all vintages from beginning to end of dates

    if pseudo
      start = '1776-07-04';
      stop  = '9999-12-31';
    else
      start = vintdates(1);
      stop  = vintdates(end);
    end
    vintall = fred.vintrange(series, start, stop, varargin{:});


  %% If didn't download correctly quit

    if isempty(vintall.value)
      returned = vintall;
      return
    end

  %% If it did download, get values for given vintage dates and return

    % Get publication lags for pseudo vintage
    if pseudo
      publag = fred.computePublag_(vintall.date, vintall.realtime, vintall.value);
    end

    % Select data column for each given vintage date from the matrix of
    % all vintages
    transfer = {'info', 'series', 'frequency', 'units'};
    for n = 1:length(transfer)
      returned.(transfer{n}) = vintall.(transfer{n});
    end
    returned.pseudo   = nan(length(vintdates),1);
    returned.realtime = vintdates;
    returned.date = vintall.date;

    Nobs  = length(returned.date);
    Nvint = length(returned.realtime);
    returned.value = nan(Nobs, Nvint);
    for n = 1:Nvint
      vint = returned.realtime(n);
      col = find(vint >= vintall.realtime, 1, 'last');
      if ~isempty(col)
        returned.value(:,n) = vintall.value(:,col);
      elseif pseudo
        % What would have been available: Anything >= publag days after
        % the observation date
        available = ((vint-vintall.date) >= publag.median);
        returned.value(find(available),n) = vintall.value(find(available),1);
        returned.pseudo(n) = vintall.realtime(1); % First available date
      end
    end

  %% Chop off trailing nan rows, which arise for pseudo vintages

  if pseudo
    [returned.value, trailNaNRows] = fred.RemLeadTrailNaN_(returned.value, 'trail');
    returned.date(find(trailNaNRows)) = [];
  end



end
