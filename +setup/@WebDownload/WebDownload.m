classdef WebDownload
    %WEBDOWNLOAD
    %   Utility to perform file downloads from web.
    
    properties
        url = ''            % URL of the zip file to be downloaded.
        zipFilepath = ''    % Local filepath to download file.
        folderPath  = ''    % Local folder path to unzip downloaded file.
    end
    
    properties (Dependent)
        zipExists           % Checks if zip file exists.
        folderExists        % Checks if folder with unzipped files exist.
    end
    
    methods
        function self = WebDownload(varargin)
            %
            % Web download instance.
            %
            % INPUTS:
            %   1. url         - url of the file to download.
            %   2. zipFilepath - local filepath to download file.
            %   3. folderPath  - local folder path to unzip downloaded
            %   file.
            %
            narginchk(0,3);
            
            if ~isempty(varargin)
                if nargin >= 3
                    self.folderPath  = varargin{3};
                end
                
                if nargin >= 2
                    self.zipFilepath = varargin{2};
                end
                
                if nargin >= 1
                    self.url = varargin{1};
                end
            end
        end
        
        function self = set.url(self,newURL)
            self.url = newURL;
        end
        
        function self = set.zipFilepath(self,newFilepath)
            self.zipFilepath = newFilepath;
        end
        
        function self = set.folderPath(self,newFolderPath)
            self.folderPath = newFolderPath;
        end
        
        function exists = get.zipExists(self)
            %
            % Validate if zip file already exists.
            %
            
            % If filepath is empty, then assume that zip already exists
            if isempty(self.zipFilepath)
                exists = true;
                
            % Otherwise, verify if a file or folder with the zip name
            % already exists.
            elseif exist(self.zipFilepath,'file') && ~exist(self.zipFilepath,'dir')
                exists = true;
            else
                exists = false;
            end
        end
        
        function exists = get.folderExists(self)
            %
            % Validate if folder already exists.
            %

            % If folder path is empty, then assume that folder already 
            % exists
            if isempty(self.folderPath)
                exists = true;
                
            % Otherwise, verify if a folder with the folder name already
            % exists.
            elseif exist(self.folderPath,'dir')
                exists = true;
            else
                exists = false;
            end
        end
    end
    
    methods
        function download(self)
            %
            % Download file from a given web address. This method downloads
            % the zip source file and unzips it. If you do not want to
            % manage the file lifecycle, you should use this method.
            %
            
            % If url is not provided, throw error to user
            if isempty(self.url)
                error = MException('WebDownload:download:noURL', 'URL is not defined...');
                error.throw();
            end
                
            % If file has already been downloaded and unzipped, do nothing.
            if ~self.folderExists
                % Download file source if not yet downloaded.
                if ~self.zipExists
                    self.get();
                end

                % Unzip source and delete temporary files
                self.unzip();
                self.deleteTemporaryFiles();
            end
        end
        
        function get(self)
            %
            % Download source from url. The source should be in zipped
            % format. The recomended method to be used to manage the file
            % lifecycle is the download method.
            %

            % If url or local filepath are not provided, throw error to user
            if isempty(self.url) || isempty(self.zipFilepath)
                error = MException( 'WebDownload:get:insufficientInformation' ...
                                  , 'URL and/or local filepath are not defined...' );
                error.throw();
            end
            
            % For earlier Matlab releases than R2014b, use urlwrite.
            if verLessThan('matlab','8.4')
                dump = urlwrite(self.url,self.zipFilepath);
            else
                dump = websave(self.zipFilepath,self.url);
            end
        end
        
        function unzip(self)
            %
            % Unzip source file to a given folder. The recomended method to
            % be used to manage the file lifecycle is the download method.
            %
            
            % If local filepath and unzip folder are not provided, throw 
            % error to user
            if isempty(self.folderPath) || isempty(self.zipFilepath)
                error = MException( 'WebDownload:unzip:insufficientInformation' ...
                                  , 'Local folder path and/or filepath are not defined...' );
                error.throw();
            end
            
            unzip(self.zipFilepath,self.folderPath);
        end
        
        function deleteTemporaryFiles(self)
            %
            % Delete temporary files created to download source dependency.
            % The recomended method to be used to manage the file lifecycle
            % is the download method.
            %
            
            if exist(self.zipFilepath,'file') && ~exist(self.zipFilepath,'dir')
                delete(self.zipFilepath);
            end
        end
    end
end
