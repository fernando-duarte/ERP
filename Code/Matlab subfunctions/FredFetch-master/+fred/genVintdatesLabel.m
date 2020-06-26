function [ lab ] = genVintdatesLabel(frequency_code)
%% genVintDatesLabel.m - Labels for the codes in genVintdates.m
%
% Given a letter code that you would pass to gen_xday_vints, this will
% generate a letter label like "m", "d", "y" based on whether the the
% frequency code would generate monthly, daily, yearly, etc. vintage
% dates.

  switch frequency_code
  case 0
    lab = 'd'; % daily
  case 1
    lab = 'w'; % Weekly
  case 2
    lab = 'm'; % monthly
  case 3
    lab = 'q'; % quarterly
  case 4
    lab = 'y'; % yearly
  end

end

