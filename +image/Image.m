classdef Image < abstract.TemplateImage
    %IMAGE
    %   Image methods and properties.

    properties (Dependent)
        size            % Image size
    end

    methods
        function self = Image(varargin)
            %
            % Image instance
            %
            % INPUTS:
            %   1. imageData - image data.
            %
            narginchk(0,1);
            
            % Create super class instance
            self = self@abstract.TemplateImage();
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
        
        function imageSize = get.size(self)
            % If image data is empty, consider the default image size
            % values.
            if isempty(self.data)
                imageSize = image.ImageSize();
            else
                imageSize = image.ImageSize(size(self.data));
            end
        end
    end
    
    methods 
        function self = flipPixels_u(self)
            %
            % Flip u-pixel coordinate of image.
            %
            matlabImage = self.MATLAB_FORMAT.encode(self.data);
            matlabImage = flip(matlabImage,self.MATLAB_FORMAT().pixelDimension_u);
            self.data   = self.MATLAB_FORMAT.decode(matlabImage);
        end
        
        function self = flipPixels_v(self)
            %
            % Flip v-pixel coordinate of image.
            %
            matlabImage = self.MATLAB_FORMAT.encode(self.data);
            matlabImage = flip(matlabImage,self.MATLAB_FORMAT().pixelDimension_v);
            self.data   = self.MATLAB_FORMAT.decode(matlabImage);
        end
    end
    
    methods (Static)
        function self = ImageFromVideo(videoInstance,framePosition)
            %
            % Create image instance from video data.
            %
            % INPUTS:
            %   1. videoInstance - video instance.
            %   2. framePosition - frame position in video.
            %
            narginchk(1,2);
            
            % Convert input data to instances
            videoInstance = utils.enums.Classes.VIDEO().convert(videoInstance);

            if nargin <= 1
                framePosition = round(videoInstance.size.numberFrames / 2);
            end
            
            % Create image instance with frame
            self = image.Image(videoInstance.data(:,:,:,framePosition));
        end

        function self = ImageFromPixels(imagePixels,imageSize)
            %
            % Create image instance from projected image pixels and
            % considering the image size.
            %
            % INPUTS:
            %   1. imagePixels - projected points on the image plane. These
            %      points should not be provided in homogeneous
            %      coordinates. Each pixel should be given in a different
            %      column.
            %   2. imageSize   - size of the image plane.
            %
            narginchk(2,2);
            
            % If pixels are provided in rows instead of columns, transpose
            % the pixels
            imagePixels = utils.enums.Classes.PIXEL().convert(imagePixels);
            imageSize   = utils.enums.Classes.IMAGE_SIZE().convert(imageSize);

            % Create black image
            imageData = zeros(imageSize.data);
            
            % Round pixels to the nearest integer
            imagePixels.data = round(imagePixels.data);
            
            % Obtain pixels that are valid for image size
            validPixels_u = imagePixels.u >= 1 & imagePixels.u <= imageSize.numberPixels_u;
            validPixels_v = imagePixels.v >= 1 & imagePixels.v <= imageSize.numberPixels_v;
            imagePixels   = imagePixels.obtainVectors(validPixels_u & validPixels_v);
            
            % If number vectors obtained are greater or equal to one, fill
            % information on pixels.
            if imagePixels.numberVectors > 0
                % Obtain linear indices for valid image pixels considering the
                % number of colors in image
                linearIndices = [];
                for iChannel = 1:imageSize.numberChannels
                    linearIndicesForColor = sub2ind( imageSize.data ...
                                                   , iChannel * ones(1, imagePixels.numberVectors ) ...
                                                   , imagePixels.u ...
                                                   , imagePixels.v );
                    linearIndices = cat(2,linearIndices,linearIndicesForColor);
                end

                % Fill image with pixels obtained from projection and create
                % image instance
                imageData(linearIndices) = 1;
            end
            
            % Create image instance
            self = image.Image(imageData);
        end
    end
end