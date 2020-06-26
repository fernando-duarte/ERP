function sas2csv(wrds, libdataname, outfile)
% SAS2CSV Download SAS data set as zipped CSV
%
%   SAS2CSV(CONN, LIBDATANAME) Where CONN should be a valid WRDS connection
%                              and LIBDATANAME should be a string with the 
%                              SAS library and dataset in the format
%                              <libref>.<data set>, e.g. 'CRSPA.MSI'.
%
%   SAS2CSV(..., OUTFILE)
%
% See also: WRDS, UNZIP, CSVREAD, READTABLE

if wrds.isVerbose, fprintf('Retrieving ''%s''.\n', libdataname), end

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

% Import sas command
fid = fopen(fullfile(wrds.Fullpath, 'sas','sas2csv.sas'));
str = fread(fid,'*char')';
fclose(fid);

% Create .zip temp name and 8char uuid for the fileref
fulluuid = ['f', strrep(char(java.util.UUID.randomUUID),'-','_')];
tmpzip   = sprintf('~/tmp/%s.zip',fulluuid);
uuid     = fulluuid(1:8);

% Fill in dataset and temp file names
sascmd = sprintf(str, uuid, tmpzip, libref, dtname, libdataname, uuid);
sascmd = regexprep(sascmd,'\*[^\n\r]*[\n\r]*','');      % strip comments
sascmd = regexprep(sascmd,'[ \t]*',' ');                % multiple spaces to one
sascmd = regexprep(sascmd,'[\n\r]*','\\n');             % newlines to literal \n
sascmd = regexprep(sascmd,'''','\\047');                % single quote ' to octal representation \047

% Build command
% mkdir -p tmp
cmd = sprintf(['rm tmp/sas2csv.sas;'...                     % Delete
    'touch "~/tmp/sas2csv.sas";',...                        % Create
    'printf ''%s'' > ~/tmp/sas2csv.sas;',...                % Write sas command
    'sas tmp/sas2csv.sas -log ~/tmp/report.log && ',...     % Execute sas
    'printf "@ -\\n@=%s.csv\\n" | zipnote -w %s',...        % Rename file in zip 
    ],sascmd, libdataname, tmpzip);

% Execute through ssh
if wrds.isVerbose, fprintf('Request submitted to WRDS servers.\n'), end
wrds.cmd(cmd,false);

% Transfer the data
try 
    wrds.getFile(tmpzip, outfile);
    ME = [];
catch ME
end

% Cleanup
wrds.cmd(sprintf('rm %s',tmpzip),false);

if ~isempty(ME)
    rethrow(ME)
end
end