classdef Dependency < setup.WebDownload
    %DEPENDENCY
    %   Utility to download and add dependencies to matlab.
    
    properties
        sourcePath = {}         % Source paths to be included on matlab path
    end
    
    properties (Dependent)
        inPath                  % Checks if source paths are on matlab path
        numberSourcePaths       % Number of source paths
    end
    
    methods
        function self = Dependency(varargin)
            %
            % Dependency instance.
            %
            % INPUTS:
            %   1. url         - url of the file to download.
            %   2. zipFilepath - local filepath to download file.
            %   3. folderPath  - local folder path to unzip downloaded
            %   file.
            %   4. sourcePath  - local folder paths to be included in
            %   Matlab path. If no path is provided all folders in
            %   folderPath will be added to the Matlab path.
            %
            narginchk(0,4);
            
            % Create super class instance
            self@setup.WebDownload();
            
            if ~isempty(varargin)
                if nargin >= 4
                    self.sourcePath  = varargin{4};
                end
                
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
        
        function self = set.sourcePath(self,newSourcePath)
            self.sourcePath = newSourcePath;
        end
        
        function number = get.numberSourcePaths(self)
            number = length(self.sourcePath);
        end
        
        function inPath = get.inPath(self)
            %
            % Validate if package is already on Matlab path.
            %
            
            inPath     = true;
            matlabPath = path;
            if self.numberSourcePaths == 0
                % If source path is not available, search for parent folder
                % if it is provided. If it is not provided assume that the
                % path is on matlab path.
                if ~isempty(self.folderPath) && ...
                    isempty(strfind(matlabPath,self.folderPath))
                    inPath = false;
                end
            else
                % If source is available, search source folders
                for iSource = 1:self.numberSourcePaths
                    if isempty(strfind(matlabPath,self.sourcePath{iSource}))
                        inPath = false;
                        break
                    end
                end
            end
        end
    end
    
    methods
        function install(self)
            %
            % Install package. This installation process will download the
            % source file to a specific folder and add this folder to the
            % matlab path.
            %
            
            % If package is not in Matlab path, package is installed.
            if ~self.inPath
                % Download package if url is defined
                if ~isempty(self.url)
                    self.download();
                end
                
                % Add package to path
                self.addPackageToPath();
            end
        end
        
        function addPackageToPath(self)
            %
            % Add package to Matlab path. If a source path is defined add
            % only the source path, otherwise include all folders and
            % subfolders from package.
            %
            
            if self.numberSourcePaths == 0
                % Add all folders and subfolders if folder path is defined
                if ~isempty(self.folderPath)
                    packagePaths = genpath(self.folderPath);
                    addpath(packagePaths);
                end
            else
                % Add only the source folders
                for iSource = 1:self.numberSourcePaths
                    addpath(self.sourcePath{iSource});
                end
            end
            savepath;
        end
    end
end
