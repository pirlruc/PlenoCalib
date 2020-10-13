classdef EpipolarPlaneImageSize < abstract.TemplateSize
    %EPIPOLARPLANEIMAGESIZE
    %   Epipolar plane image size properties.
    
    properties (Dependent)
        numberChannels   	% Number of channels in epipolar plane image
        numberPixels        % Number of pixels in epipolar plane image
        numberFrames        % Number of frames used to obtain the epipolar plane image
    end
    
    methods
        function self = EpipolarPlaneImageSize(varargin)
            %
            % Epipolar plane image size instance.
            %
            % INPUTS:
            %   1. epipolarPlaneImageSize - array with information of 
            %   epipolar plane image size.
            %
            narginchk(0,1)
            
            self@abstract.TemplateSize(varargin{:});
        end
        
        function number = get.numberChannels(self)
            number = self.numberItemsInDimension(1);
        end
        
        function number = get.numberPixels(self)
            number = self.numberItemsInDimension(2);
        end
        
        function number = get.numberFrames(self)
            number = self.numberItemsInDimension(3);
        end        
    end
end