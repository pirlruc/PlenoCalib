classdef ImageSize < abstract.TemplateSize
    %IMAGESIZE
    %   Utility to define the image size structure.
    
    properties (Dependent)
        numberPixels_u      % Number of horizontal pixels in image
        numberPixels_v      % Number of vertical pixels in image
        numberChannels      % Number of channels in image
    end
    
    methods
        function self = ImageSize(varargin)
            %
            % Image size instance.
            %
            % INPUTS:
            %   1. imageSize - array with information of image size.
            %
            narginchk(0,1)
            
            self@abstract.TemplateSize(varargin{:});
        end
        
        function number = get.numberPixels_u(self)
            number = self.numberItemsInDimension(2);
        end
        
        function number = get.numberPixels_v(self)
            number = self.numberItemsInDimension(3);
        end
        
        function number = get.numberChannels(self)
            number = self.numberItemsInDimension(1);
        end
    end
end