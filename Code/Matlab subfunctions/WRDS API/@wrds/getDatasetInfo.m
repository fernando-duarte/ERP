function info = getDatasetInfo(wrds, libdataname)
% GETDATASETNAMES Retrieve dataset names within a library
%
%    INFO = GETDATASETINFO(WRDS, LIBDATANAME) LIBDATANAME should be a 
%                              string with the SAS library and dataset in 
%                              the format <libref>.<data set>, 
%                              e.g. 'CRSPA.MSI'.
%               
% See also: WRDS, GETLIBREFS, GETDATASETNAME

% Library and datasetname
tmp    = regexp(libdataname, '\.','split');
libref = tmp{1};
dtname = upper(tmp{2});

% Sanitize input
try
    allLib = wrds.getLibrefs;
    idx    = strcmpi(libref, allLib);
    libref = allLib{idx};
catch ME
end

if wrds.isVerbose, fprintf('Retrieving dataset info for ''%s''.\n', libdataname), end

% SAS command
sascmd = sprintf(['FILENAME out "~/tmp/cmd.lst";',...
    'PROC PRINTTO print=out;',...
    'RUN;',...
    'PROC DATASETS LIBRARY=%s NOLIST;',...
    'CONTENTS DATA=%s;',...
    'RUN;'],...
    libref, dtname);
% UNIX command
cmd    = sprintf(['touch ~/tmp/cmd.sas;',...                % Create file
    'printf ''%s'' > ~/tmp/cmd.sas;',...                    % Write sas command
    'sas ~/tmp/cmd.sas -log ~/tmp/cmd.log;',...             % Execute sas
    'cat ~/tmp/cmd.lst;'],...                               % Print file
    sascmd);

% Execute through ssh
if wrds.isVerbose, fprintf('Request submitted to WRDS servers.\n'), end
[~,result] = wrds.cmd(cmd,false);
info       = char(result);

% Cleanup
wrds.cmd('rm ~/tmp/cmd.*');

end