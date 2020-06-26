function librefs = getLibrefs(wrds)
% GETLIBDATANAMES Retrieve all SAS library names
%
%    LIBREFS = GETLIBDATANAMES(WRDS) It retrieves the libref part 
%                                    in <libref>.<data set>, e.g. 
%                                    CRSPA.
%               
% See also: WRDS, SAS2CSV

librefs = wrds.Librefs;

if isempty(librefs)
    if wrds.isVerbose, fprintf('Retrieving SAS library names (libref).\n'), end
    
    % Build command
    sascmd = 'libname _all_ list;';
    cmd = sprintf(['touch "~/tmp/cmd.sas";',...                 % Create file
        'printf ''%s'' > ~/tmp/cmd.sas;',...                    % Write sas command
        'sas ~/tmp/cmd.sas -log ~/tmp/cmd.log;',...             % Execute sas
        'grep ''NOTE: Libref='' ~/tmp/cmd.log | ',...
        'sed ''s/^NOTE: Libref= *//'' | ',...                     
        'sed ''s/ *//g'';'],...                                  % Parse log for librefs
        sascmd);
    
    % Execute through ssh
    if wrds.isVerbose, fprintf('Request submitted to WRDS servers.\n'), end
    [~,result] = wrds.cmd(cmd,false);
    
    % Cleanup
    wrds.cmd('rm ~/tmp/cmd.*');
    
    % Store in wrds
    librefs = sort(result);
    wrds.Librefs = librefs;
end

end