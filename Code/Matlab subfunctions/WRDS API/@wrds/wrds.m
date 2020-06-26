classdef wrds < handle
    % WRDS Connect to Wharton Reasearch Data Services 
    %
    %   High level Matlab API that interacts with the WRDS Unix server 
    %   and the SAS data sets through SSH2.
    %   
    %   Requirements:
    %       - An account with WRDS of the type that admits SSH connections. 
    %         See <a href="http://wrds-web.wharton.upenn.edu/wrds/support/Additional%20Support/Account%20Types.cfm">account types</a> for details.
    %       - <a href="matlab: disp(['Java enabled: ' isempty(javachk('jvm'))+'0'])">Java enabled</a>
    %
    %   Syntax:
    %       WRDS(USERNAME, PASS) Supply USERNAME and PASS as strings. 
    %
    %       WRDS(..., HOST, PORT) Optionally, provide HOST and/or PORT 
    %                             which are respectively defaulted to 
    %                             'wrds.wharton.upenn.edu' and 22. 
    %   
    %
    %       W = WRDS(...) Connection to the server
    %
    %
    %   Examples:
    %       w = username('myuser','forgiveMeIfIDontTellYou');
    %       w.cmd('echo "Hello World!"')
    %   
    % See also: SSH2
    
    
    properties
        isVerbose@logical   = true;           % Toggle verbosity
        isConnected@logical = false;        % Track connection status
    end
    properties (Access=private)
        SSH2conn                            % SSH2 connection
        Fullpath                            % Full path to wrds folder
        Librefs
        Libdatasets
    end
    
    methods
        function obj = wrds(username, pass, host, port)
            % WRDS Constructor
             
            % Only check jvm, other errors delegated
            error(javachk('jvm'))
            
            % Defaults initialization
            if nargin < 1 || isempty(username), 
                tmp      = passdlg('ups'); 
                username = tmp.User{1};
                pass     = tmp.Pass{1}; 
            elseif nargin < 2 || isempty(pass), 
                tmp  = passdlg('ps'); 
                pass = tmp.Pass{1}; 
            end
            if nargin < 3 || isempty(host), host = 'wrds.wharton.upenn.edu';     end
            if nargin < 4 || isempty(port), port = 22;                           end
            
            % Establish ssh2 connection
            obj.SSH2conn = ssh2_config(host, username, pass, port);
            if ~isempty(username)
                obj.SSH2conn = ssh2_main(obj.SSH2conn);
            end
            
            % Record where the wrds path is
            obj.Fullpath = regexprep(fileparts(mfilename('fullpath')),'\@wrds','');
            
            % Check if connected
            obj.isConnected = ~isempty(obj.SSH2conn.connection);
            
            % Initializations
            if obj.isConnected
                fprintf('Connected.\n')
                obj.cmd('cd ~; mkdir -p tmp',false);
                obj.Librefs = getLibrefs(obj);
            else
                fprintf('Could not connect.\n')
            end
        end
        
        function [obj, result] = cmd(obj, cmdstr, isVerbose)
            % CMD Execute command on UNIX shell
            %
            % Example:
            %       w = username('myuser','forgiveMeIfIDontTellYou');
            %       w.cmd('echo "Hello World!"')
            if nargin < 3, isVerbose = obj.isVerbose; end

            [obj.SSH2conn, result] = ssh2_command(obj.SSH2conn,cmdstr,isVerbose);
        end
        
        function [obj, outfile] = getFile(obj, remotefile, outfile)
            % getFile Transfer file from remote host by Secure Copy
            
            if nargin < 3 || isempty(outfile)
                outfile = fullfile(obj.Fullpath,'data\');
            end
            
            % Process paths
            [rpath, rfname, rext] = fileparts(remotefile);
            [lpath, lfname, lext] = fileparts(outfile);
            if isempty(lpath)
                lpath = fullfile(obj.Fullpath,'data\');
            end
            if isempty(lfname)
                lfname = rfname;
                lext   = rext;
            end
                        
            if obj.isVerbose, fprintf('Downloading file. Please, wait.\n'), end
            
            % Download file
            obj.SSH2conn = scp_get(obj.SSH2conn, [rfname, rext], lpath, rpath);
            
            % Rename
            if ~strcmp(lfname, rfname)
                outfile = fullfile(lpath, [lfname, lext]);
                movefile(fullfile(lpath, [rfname, rext]), outfile)
            end
        end
        
        function obj = close(obj)
            if obj.isConnected
                if obj.isVerbose, fprintf('Closing connection.\n'), end
                obj.SSH2conn = ssh2_close(obj.SSH2conn);
                obj.isConnected = false;
            end
        end
        function delete(obj)
            close(obj);
        end
    end
end