% OVERIVEW - gen_xday_vints.m
%
% Function to generate a cell of Matlab dates in format YYYY-MM-DD that
% correspond to every Xday (Friday, Monday, etc.) between startDate and today
% or some (optional) specified endDate.
%
% Can specify a frequency_code which will only return the first Xday of
% each year, each quarter, each month, etc.
%
% dow codes:
% 1 Sunday
% 2 Monday
% 3 Tuesday
% 4 Wednesday
% 5 Thursday
% 6 Friday
% 7 Saturday
%
function [ vintdates ] = genVintdates(startDate, endDate, as_datestr, frequency_code, dow)

  helper_functions;
  if ~exist('frequency_code', 'var')
    frequency_code = 0;
  end

  % All possible days between start and end date
  allPossible = [datenum(startDate):datenum(endDate)]';
  if ~exist('dow', 'var')
    dow = 6;
  end
  dowDays = allPossible(find(weekday(allPossible) == dow));

  switch frequency_code
    case -1 % Just start and end date
      vintDates = allPossible([1 end]);

    case 0 % All
      vintdates = allPossible;

    case 1 % Each dow per week
      vintdates = dowDays;

    case 2 % First dow per month
      [Y,M,~,~,~,~] = datevec(dowDays);
      YM = 100*Y + M; % Store year and month
      use = arrayfun(@(ym) find(YM == ym, 1), unique(YM));
      vintdates = dowDays(use); % Picks out the indices of the first friday of each month

    case 3 % First dow in each quarter
      quarters    = datestr(dowDays, 'yyyyqq');
      quartersAll = str2num(quarters(:,1:4)) + (str2num(quarters(:,end))-1)/4;
      quartersUnq = unique(quartersAll);
      use = arrayfun(@(q) find( quartersAll == q ,1), quartersUnq);
      vintdates = dowDays(use);

    case 4 % First dow in each year
      yearsAll = str2num(datestr(dowDays, 'yyyy'));
      yearsUnq = unique(yearsAll);
      use = arrayfun(@(q) find( yearsAll == q ,1), yearsUnq);
      vintdates = dowDays(use);
  end

  if exist('as_datestr', 'var') && as_datestr
    vintdates = datestrFred(vintdates, 1);
  end

end
