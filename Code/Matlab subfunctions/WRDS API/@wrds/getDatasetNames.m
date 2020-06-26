function dtnames = getDatasetNames(wrds, libref)
% GETDATASETNAMES Retrieve dataset names within a library
%
%    DTNAMES = GETDATASETNAMES(WRDS, LIBREF) It retrieves the data set parts 
%                                    in <libref>.<data set>, e.g. 
%                                    CRSPA.STOCKNAMES, etc...
%               
% See also: WRDS, GETLIBREFS, SAS2CSV

% Sanitize input
try
    allLib = wrds.getLibrefs;
    idx    = strcmpi(libref, allLib);
    libref = allLib{idx};
catch ME
end

try
    dtnames = wrds.Libdatasets.(libref);
catch
    if wrds.isVerbose, fprintf('Retrieving dataset names for ''%s''.\n', libref), end
    
    % SAS command
    sascmd = sprintf(['FILENAME out "~/tmp/cmd.lst";',...
        'PROC PRINTTO print=out;',...
        'RUN;',...
        'PROC SQL;',...
        'SELECT DISTINCT memname',...
        ' FROM dictionary.columns WHERE libname="%s";'],...
        libref);
    % UNIX command
    cmd = sprintf(['touch "~/tmp/cmd.sas";',...                 % Create file
        'printf ''%s'' > ~/tmp/cmd.sas;',...                    % Write sas command
        'sas ~/tmp/cmd.sas -log ~/tmp/cmd.log;',...             % Execute sas
        'grep ''^ *[A-Z_0-9]* *$'' ~/tmp/cmd.lst | ',...
        'sed ''s/ *//g'';'],...                                 % Parse .lst
        sascmd);
    
    % Execute through ssh
    if wrds.isVerbose, fprintf('Request submitted to WRDS servers.\n'), end
    [~,result] = wrds.cmd(cmd, false);
    
    % Cleanup
    wrds.cmd('rm ~/tmp/cmd.*');
    
    % Store in wrds
    dtnames = sort(result);
    wrds.Libdatasets.(libref) = dtnames;
end
end