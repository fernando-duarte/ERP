function str = pwencode(wrds,pass)
% PWENCODE Econdes the WRDS password with SAS
%
%   PWENCODE(PASS) PASS should be a string
%
%   STR = PWENCODE(...) Returns the encoded PASS
%
% See also: WRDS

% Import sas command
fid = fopen(fullfile(wrds.Fullpath, 'sas','pwencode.sas'));
str = fread(fid,'*char')';
fclose(fid);

% Create .zip temp name and 8char uuid for the fileref
fulluuid = ['f', strrep(char(java.util.UUID.randomUUID),'-','_')];
tmpfile  = sprintf('~/tmp/%s',fulluuid);

% Fill in dataset and temp file names
sascmd = sprintf(str, tmpfile, pass);
sascmd = regexprep(sascmd,'\*[^\n\r]*[\n\r]*','');      % strip comments
sascmd = regexprep(sascmd,'[ \t]*',' ');                % multiple spaces to one
sascmd = regexprep(sascmd,'[\n\r]*','\\n');             % newlines to literal \n
sascmd = regexprep(sascmd,'''','\\047');                % single quote ' to octal representation \047

% Build command
% mkdir -p tmp
cmd = sprintf(['rm ~/tmp/pwencode.sas;'...                  % Delete
    'touch ~/tmp/pwencode.sas;',...                         % Create
    'printf ''%s'' > ~/tmp/pwencode.sas;',...               % Write sas command
    'sas ~/tmp/pwencode.sas -log ~/tmp/report.log;',...     % Execute sas
    'rm ~/tmp/report.log;',...                              % Remove log
    'rm ~/tmp/pwencode.sas;',...                            % Cleanup sas 
    ],sascmd);

% Execute through ssh
wrds.cmd(cmd);

% Read the encoded pass
[~, str] = wrds.cmd(sprintf('cat %s', tmpfile));
str = str{1};

% Cleanup
wrds.cmd(sprintf('rm %s',tmpfile));
end
