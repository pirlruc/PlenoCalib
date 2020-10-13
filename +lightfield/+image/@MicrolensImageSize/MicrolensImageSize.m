classdef MicrolensImageSize < abstract.TemplateSize
    %MICROLENSIMAGESIZE
    %   Utility to define the microlens image size structure.
    
    properties (Dependent)
        numberPixels_i      % Number of horizontal pixels in image
        numberPixels_j      % Number of vertical pixels in image
        numberChannels      % Number of channels in image
    end
    
    methods
        function self = MicrolensImageSize(varargin)
            %
            % Microlens image size instance.
            %
            % INPUTS:
            %   1. microlensImageSize - array with information of microlens 
            %   image size.
            %
            narginchk(0,1)
            
            self@abstract.TemplateSize(varargin{:});
        end
        
        function number = get.numberPixels_i(self)
            number = self.numberItemsInDimension(2);
        end
        
        function number = get.numberPixels_j(self)
            number = self.numberItemsInDimension(3);
        end
        
        function number = get.numberChannels(self)
            number = self.numberItemsInDimension(1);
        end
    end
end