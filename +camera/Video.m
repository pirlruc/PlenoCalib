classdef Video < abstract.TemplateVideo
    %VIDEO
    %   Video methods and properties.
    
    methods
        function self = Video(varargin)
            %
            % Video instance.
            %
            % INPUTS:
            %   1. videoData - video data.
            %
            narginchk(0,1);
            
            % Create superclass instance
            self = self@abstract.TemplateVideo();
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
    end
    
    methods (Static)
        function convertVideoToFrames(videoFilepath)
            %
            % Create video instance from file.
            %
            % INPUTS:
            %   1. videoFilepath - video local filepath.
            %
            narginchk(1,1);
            
            % Obtain video object
            video = VideoReader(videoFilepath);
            
            % Obtain number frames
            numberFrames = floor(video.FrameRate * video.Duration) - 1;
            
            % Obtain root filename
            [root,videoFilename] = fileparts(videoFilepath);
            if ~isempty(root) 
                root = [root filesep];
            end
            for iFrame = 1:numberFrames
                % Obtain file to write image frame 
                frameFilepath = [root videoFilename sprintf('_%d',iFrame) '.png'];
                
                % Read frame data and write file
                imwrite(video.readFrame,frameFilepath);
            end
        end
        
        function self = VideoFromFile(videoFilepath)
            %
            % Create video instance from file.
            %
            % INPUTS:
            %   1. videoFilepath - video local filepath.
            %
            narginchk(1,1);
            
            % Read video data. Convert data to double.
            videoData = im2double(VideoReader(videoFilepath).read);
            
            % Obtain matlab format and decode
            matlabFormat = camera.enums.Formats.MATLAB_FORMAT();
            videoData    = matlabFormat.decode(videoData);
            
            % Create video instance
            self = camera.Video(videoData);
        end
    end
end
