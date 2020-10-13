classdef EpipolarPlaneImage < abstract.TemplateEpipolarPlaneImage
    %EPIPOLARPLANEIMAGE
    %   Epipolar plane image methods and properties.
    
    properties (Dependent)
        size            % Epipolar plane image size
    end

    methods
        function self = EpipolarPlaneImage(varargin)
            %
            % Create epipolar plane image instance.
            %
            % INPUTS:
            %   1. epiData - epipolar plane image data.
            %
            narginchk(0,1);
            
            % Create super class instance.
            self = self@abstract.TemplateEpipolarPlaneImage();
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
        
        function epiSize = get.size(self)
            % If epipolar plane image data is empty, consider the default 
            % epipolar plane image size values.
            if isempty(self.data)
                epiSize = image.epi.EpipolarPlaneImageSize();
            else
                epiSize = image.epi.EpipolarPlaneImageSize(size(self.data));
            end
        end
    end
end