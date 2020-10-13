classdef VideoSize < image.ImageSize
    %VIDEOSIZE
    %   Video size properties.
    
    properties (Dependent)
        numberFrames        % Number of frames in video
    end
    
    methods
        function self = VideoSize(varargin)
            %
            % Video size instance.
            %
            % INPUTS:
            %   1. videoSize - array with information of video size.
            %
            narginchk(0,1);
            
            self@image.ImageSize(varargin{:});
        end
        
        function numberFrames = get.numberFrames(self)
            numberFrames = self.numberItemsInDimension(4);
        end
    end
end