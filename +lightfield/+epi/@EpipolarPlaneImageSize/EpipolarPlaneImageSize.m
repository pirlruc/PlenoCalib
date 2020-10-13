classdef EpipolarPlaneImageSize < abstract.TemplateSize
    %EPIPOLARPLANEIMAGESIZE
    %   Epipolar plane image size properties for lightfield.
    
    properties (Dependent)
        numberChannels   	% Number of channels
        numberMicrolenses   % Number of microlenses (k-th coordinate). Corresponds to
                            % the pixels in the epipolar plane image.
        numberPixels        % Number of pixels (i-th coordinate). Corresponds to the
                            % frames in the epipolar plane image.
    end
    
    methods
        function self = EpipolarPlaneImageSize(varargin)
            %
            % Epipolar plane image size instance for lightfield.
            %
            % INPUTS:
            %   1. epipolarPlaneImageSize - array with information of 
            %   epipolar plane image size for lightfield.
            %
            narginchk(0,1)
            
            self@abstract.TemplateSize(varargin{:});
        end
        
        function number = get.numberChannels(self)
            number = self.numberItemsInDimension(1);
        end
        
        function number = get.numberMicrolenses(self)
            number = self.numberItemsInDimension(2);
        end        
        
        function number = get.numberPixels(self)
            number = self.numberItemsInDimension(3);
        end
    end
end