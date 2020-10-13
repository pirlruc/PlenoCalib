classdef (Abstract) TemplateVideo
    %TEMPLATEVIDEO
    %   Video template.
    
    properties (Constant)
        MATLAB_FORMAT = camera.enums.Formats.MATLAB_FORMAT()
            % Toolbox format that allows to decode and encode the video data
    end
    
    properties
        data = []           % Video data
    end
    
    properties (Dependent)
        size                % Video size
    end
    
    methods
        function self = TemplateVideo(varargin)
            %
            % Video instance.
            %
            % INPUTS:
            %   1. videoData - video data.
            %
            narginchk(0,1);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
        
        function self = set.data(self,newVideoData)
            self.data = newVideoData;
        end
        
        function videoSize = get.size(self)
            if isempty(self.data)
                videoSize = camera.VideoSize();
            else
                videoSize = camera.VideoSize(size(self.data));
            end
        end
    end
    
    methods
        function self = addFrames(self, newFrames)
            %
            % Add frames to video.
            %
            % INPUTS:
            %   1. newFrames - frames to be included in the video.
            %
            narginchk(2,2);
            
            % Add frames to video
            self.data = cat(4,self.data,newFrames);
        end
        
        function self = rgb2gray(self)
            %
            % Convert RGB values to grayscale:
            %       0.2990 * R + 0.5870 * G + 0.1140 * B 
            %
            % The rgb2gray matlab method can only be used for images.
            %

            if isempty(self.data)
                self.data = [];
            else
                self.data = 0.299 * self.data(1,:,:,:) ...
                          + 0.587 * self.data(2,:,:,:) ...
                          + 0.114 * self.data(3,:,:,:);
            end
        end
        
        function write(self, outputFilepath, frameRate)
            %
            % Write vide to output filepath considering the specified
            % frame rate.
            %
            % INPUTS:
            %   1. outputFilepath - local filepath to write video image
            %   sequences.
            %   2. frameRate - frames per second considered to write video.
            %
            narginchk(1,3);
            
            % If filepath is not provided, assume default.
            if nargin <= 1
                outputFilepath = 'video_output.avi';
            end
            
            % If frame rate is not provided, assume default value
            if nargin <= 2
                frameRate = 1/30;
            end

            % Obtain fileparts
            [filepath,filename,extension] = fileparts(outputFilepath);
            % If extension is not provided consider avi extension.
            if isempty(extension)
                extension = '.avi';
            end
            % Obtain output filepath again
            outputFilepath = [filepath, filesep, filename, extension];
            
            % Create video writer object
            outputVideo           = VideoWriter(outputFilepath);
            outputVideo.FrameRate = frameRate;

            % Write frames to file. This function needs to have the data in
            % matlab format.
            open(outputVideo);
            outputVideo.writeVideo(self.MATLAB_FORMAT.encode(self.data));
            close(outputVideo);
        end
    end
end
