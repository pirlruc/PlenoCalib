classdef (Abstract) TemplateDatasets
    %TEMPLATEDATASETS
    %   Template to download datasets from web.
    
    properties (SetAccess = immutable)
        dataset      % Datasets information for package.
    end
    
    methods
        function self = TemplateDatasets(url,zipFilepath,folderPath)
            %
            % Template datasets instance.
            %
            
            % Obtain root and dataset path
            toolbox  = setup.Setup();
            zipFilepath = [toolbox.datasetsPath filesep zipFilepath];
            folderPath  = [toolbox.datasetsPath filesep folderPath];
            
            % Create utility to download dataset files
            self.dataset = setup.WebDownload(url,zipFilepath,folderPath);
        end
    end
    
    methods
        function download(self)
            %
            % Download dataset.
            %
            self.dataset.download();
        end
    end
end