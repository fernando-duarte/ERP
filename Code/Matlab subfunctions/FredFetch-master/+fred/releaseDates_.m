function [returned] = releaseDates_(input)
% [returned] = releaseDates_(series) will return release dates for the
% given series.
%
% [returned] = releaseDates_(releaseNumber) will return dates for the
% given Fred release number.


  %% Get the release ID for the given series

    % It's a series name; get the release id number
    if ischar(input)
      release_id = fred.releaseID_(input);
    else
      release_id = input;
    end

  %% Download the release dates for that release ID

    opt = fred.GlobalOptions();
    datesURL = sprintf([...
            'https://api.stlouisfed.org/fred/release/dates?' ...
            'release_id=%d' ...
            '&api_key=%s' ...
            '&include_release_dates_with_no_data=true' ...
            '&file_type=json'],...
            release_id, ...
            opt.api);
    fromFred = fred.ReadFredURL_(datesURL, 1, opt.max_attempt);
    fromFred = [fromFred.release_dates{:}];


    returned.release_id = release_id;
    returned.date       = fred.dtnum({fromFred.date}, 1)';

end
