function [dt] = dtGivenFrequency_(frequency, start, stop)
% dtgivenfrequency_ - Produce all period dates at a given frequency
%                     between a starting and stopping date

  switch frequency
    case 'w'
      % Every Friday in between the dates
      dow   = 6; % 6 is the matlab code for friday when you run weekday()
      dtall = [start-7:stop+7]';
      dt    = dtall(find(weekday(dtall) == dow));

    case 'm'
      dtY = fred.dtfld([start:stop]', 'year');
      dtM = fred.dtfld([start:stop]', 'month');
      dt = unique(datenum(dtY, dtM, 1));

    case 'q'
      dtY = fred.dtfld([start:stop]', 'year');
      dtQ = fred.dtfld([start:stop]', 'quarter');
      dt = unique(datenum(dtY, 3*dtQ-2, 1));

    case 'a'
      dt = datenum(unique(fred.dtfld([start:stop]', 'year')), 1, 1);

  end

end
