classdef TemplatePlenopticCamera
    %TEMPLATEPLENOPTICCAMERA
    %   Template camera utility to decode lightfields and calibrate camera.

    properties (Constant)
        POSE_CODE_SCALE      = 100  % Scale for creating pose code
        LENS_CODE_SCALE      = 10   % Scale for creating lens code
        DIMENSION_CODE_SCALE = 1    % Scale for creating dimension code
    end
    
    properties
        whiteImagesFilepath = ''    % Local filepath for the white image file
    end
    
    properties (Abstract, Dependent)
        numberMicrolensesTypes      % Number of different microlenses types
    end
    
    methods
        function self = TemplatePlenopticCamera(varargin)
            %
            % Template plenoptic camera instance.
            %
            % INPUTS:
            %   1. whiteImageFilepath - filepath to white image file.
            %
            narginchk(0,1);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.whiteImagesFilepath = varargin{1};
                end
            end
        end
        
        function self = set.whiteImagesFilepath(self,newWhiteImagesFilepath)
            self.whiteImagesFilepath = newWhiteImagesFilepath;
        end
    end
end