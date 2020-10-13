classdef Setup
    %SETUP
    %   Utility to setup package dependencies.
    
    properties (Constant, Hidden)
        ROOT_FILE = 'setup.dummy';                  % Root file to determine toolbox path
        CLASSNAME = 'setup.enums.dependencies.%s';  % Location of the package dependencies
    end
    
    properties
        rootPath                            % Root path of the computer vision toolbox
        dependenciesPath                    % Package dependencies path of the computer vision toolbox
        toolboxesPath                       % Local path to store the dependencies of the computer vision toolbox
        datasetsPath                        % Local path to store the datasets of the computer vision toolbox
    end
    
    properties (Dependent)
        packages                            % Packages of toolbox
        numberPackages                      % Number of packages on toolbox
    end
    
    methods
        function self = Setup()
            %
            % Setup instance.
            %
            narginchk(0,0);
            
            % This action returns the setup package path
            self.rootPath = fileparts(which(self.ROOT_FILE));
            
            % Define the dependencies path
            self.dependenciesPath = [self.rootPath filesep '+enums' filesep '+dependencies'];
            
            % Since we want the toolbox root path, remove setup from path
            self.rootPath = fileparts(self.rootPath);
            
            % Define toolboxes and datasets path
            self.toolboxesPath = [self.rootPath filesep 'toolboxes'];
            self.datasetsPath  = [self.rootPath filesep 'data'];
        end
        
        function packages = get.packages(self)
            %
            % Obtain dependency packages information.
            %
            packages = what(self.dependenciesPath);
            packages = packages.m;
        end
        
        function number = get.numberPackages(self)
            number = length(self.packages);
        end
    end
    
    methods
        function install(self)
            %
            % Create toolbox folder structure and install toolbox 
            % dependencies.
            %
            
            % Validate if dataset folder already exists. If not, create it.
            if ~exist(self.datasetsPath,'dir')
                mkdir(self.datasetsPath);
            end
            
            % Validate if toolbox folder already exists. If not, create it.
            if ~exist(self.toolboxesPath,'dir')
                mkdir(self.toolboxesPath);
            end
            
            % Install dependencies for toolbox.
            for iPackage = 1:self.numberPackages
                % Obtain package name
                package = self.packages{iPackage};
                [~,packageName,~] = fileparts(package); % Remove extension

                % Obtain package dependencies and install them
                dependencies = enumeration(sprintf(self.CLASSNAME,packageName));
                for iDependency = 1:length(dependencies)
                    dependency  = dependencies(iDependency);
                    dependency.install();
                end
            end
        end
    end
end

