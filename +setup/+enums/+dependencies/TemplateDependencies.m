classdef (Abstract) TemplateDependencies
    %TEMPLATEDEPENDENCIES
    %   Template package dependencies.
    
    properties (SetAccess = immutable)
        dependency      % Dependencies information for package.
    end
    
    methods
        function self = TemplateDependencies(url,filepath,packagePath,sourcePath)
            %
            % Template package dependencies instance.
            %
            
            % Obtain root and package path
            toolbox  = setup.Setup();
            filepath = [toolbox.toolboxesPath filesep filepath];
            packagePath = [toolbox.toolboxesPath filesep packagePath];

            % Obtain source paths to add to matlab path
            newSourcePath = cell(length(sourcePath),1);
            for iSource = 1:length(sourcePath)
                newSourcePath{iSource} = [toolbox.toolboxesPath filesep sourcePath{iSource}];
            end
                        
            % Create utility to install dependency
            self.dependency = setup.Dependency(url,filepath,packagePath,newSourcePath);
        end
    end
    
    methods 
        function install(self)
            %
            % Install dependency.
            %
            self.dependency.install();
        end
    end
end