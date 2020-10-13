classdef ViewpointImageSize < abstract.TemplateSize
    %VIEWPOINTIMAGESIZE
    %   Utility to define the viewpoint image size structure.
    
    properties (Dependent)
        numberMicrolenses_k         % Number of horizontal microlenses
        numberMicrolenses_l         % Number of vertical microlenses
        numberChannels              % Number of channels in image
    end
    
    methods
        function self = ViewpointImageSize(varargin)
            %
            % Viewpoint image size instance.
            %
            % INPUTS:
            %   1. viewpointImageSize - array with information of viewpoint
            %   image size.
            %
            narginchk(0,1)
            
            self@abstract.TemplateSize(varargin{:});
        end
        
        function number = get.numberMicrolenses_k(self)
            number = self.numberItemsInDimension(2);
        end
        
        function number = get.numberMicrolenses_l(self)
            number = self.numberItemsInDimension(3);
        end
        
        function number = get.numberChannels(self)
            number = self.numberItemsInDimension(1);
        end
    end
end