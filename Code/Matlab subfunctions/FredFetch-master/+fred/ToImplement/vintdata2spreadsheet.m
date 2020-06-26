% OVERIVEW - vintdata2spreadsheet.m
%
% Takes a vintage structure as produced by align_vintdata.m and writes it to a
% csv file with the series names and transformations (if they have been
% applied).
%
% Not done automatically when pulling the data bc a vintage file for a single
% data in 2014 is roughly 1.1MB,  too big to do automatically, especially since
% we have the matfile
%
% Now takes argument num_nan, which is the number of nan rows to append
% to the end of the vintdata spreadsheet
function [ ] = vintdata2spreadsheet(outfile, vintdata, num_nan)

  %vintdata
  %load(vintdata)

  %% Possibly extend data:
  % - Appendiing num_nan rows of nans  to the end of spreadsheet and
  %   extending date vector by num_nan months
  if exist('var', 'num_nan')
    vintdata.values = [vintdata.values; nan(num_nan,size(vintdata.values,2))];
    vintdata.dates  = [vintdata.dates;  arrayfun(@(n) addtodate(vintdata.dates(end), n, 'month'), 1:num_nan)];
  end


  f = fopen(outfile, 'w');

  % Write series names
  fprintf(f, 'series,');
  fprintf(f, '%s,', vintdata.series{:});
  fprintf(f, '\n');

  % Write transformation
  if vintdata.transf_ind
    fprintf(f, 'transformation,');
    fprintf(f, '%s,', vintdata.transf{:});
    fprintf(f, '\n');
  end

  % Write data
  for dt = 1:length(vintdata.dates)
    fmt = [datestr(vintdata.dates(dt)), repmat(',%14.5f', 1, size(vintdata.values,2))];
    fprintf(f, fmt, vintdata.values(dt,:) );
    fprintf(f, '\n');
  end

  fclose(f);

end
