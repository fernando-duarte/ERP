function [opt, toPass] = parseVarargin_(nameFlds, varargin)
% PARSE_VARARGIN - Check for various flags in varargin and return stuff
%                  to pass as additional options to FRED API
%
% [PAR, PSEUDO, VARAGOUT] = parse_varargin(VARARGIN) will parse the
% varargin and check for arguments of type boolean or logical, which
% would follow a string denoting either 'parallel' (to download series
% in parallel) or 'pseudo' (to construct pseudo vintages).
%
% The rest of varargin should just be strings that will be passed along
% as additional options to the FRED API. They are returned as VARARGOUT
% so you can use them (without the par and pseudo flags)

  %% toDatasetByVint option

    % If the first entry is a straight number, no options, it represents
    % "toDatasetByVint" flag
    opt.toDatasetByVint = 1;
    if ~isempty(varargin) && isnumeric(varargin{1})
      opt.toDatasetByVint = varargin{1};
      varargin = varargin(2:end);
    end


  %% Now check for parallel and pseudo-vintage flags

    rem = zeros(1,2*length(varargin));
    for n = 1:length(nameFlds)
      ind = find(strcmp(nameFlds{n}, varargin));
      if isempty(ind)
        opt.(nameFlds{n}) = 0;
      else
        opt.(nameFlds{n}) = varargin{ind+1};
        rem([ind ind+1]) = 1;
      end
    end

    % Return varargin, with the parallel and pseudo shit stripped out
    toPass = varargin(setdiff(1:length(varargin), find(rem)));

end
