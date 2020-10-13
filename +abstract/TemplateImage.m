classdef (Abstract) TemplateImage
    %TEMPLATEIMAGE
    %   Image template.

    properties (Constant)
        MATLAB_FORMAT = image.enums.Formats.MATLAB_FORMAT()
            % Toolbox format that allows to decode and encode the image data
    end
    
    properties
        data = []       % Image data
    end
    
    properties (Abstract, Dependent)
        size            % Image size
    end
    
    methods
        function self = TemplateImage(varargin)
            %
            % Template image instance
            %
            % INPUTS:
            %   1. imageData - image data.
            %
            narginchk(0,1);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
        
        function self = set.data(self,newImageData)
            self.data = newImageData;
        end
    end
    
    methods
        function self = correctColor(self, colorCorrectionMethod)
            %
            % Correct color from image using a given color correction
            % method.
            %
            % INPUTS:
            %   1. colorCorrectionMethod - color correction method. The
            %   default is the shades of gray method.
            %
            narginchk(1,2);
            
            if nargin <= 1
                colorCorrectionMethod = image.enums.ColorCorrectionMethods.SHADES_OF_GRAY();
            end
            
            % Normalize image color
            self = colorCorrectionMethod.method(self);
        end
        
        function self = rgb2gray(self)
            %
            % Convert RGB values to grayscale:
            %       0.2990 * R + 0.5870 * G + 0.1140 * B 
            %
            
            % If no data is provided, return empty list
            if isempty(self.data)
                self.data = [];
            elseif self.size.numberChannels > 1
                matlabImage = self.MATLAB_FORMAT.encode(self.data);
                self.data   = self.MATLAB_FORMAT.decode(rgb2gray(matlabImage));
            end
        end
        
        function self = applyHomography(self,homography,imageSize,showSamples)
            %
            % Obtain target image from applying a given homography mapping
            % to the source image (self object).
            %
            % INPUTS:
            %   1. homography  - homography matrix.
            %   2. imageSize   - resolution of target image.
            %   3. showSamples - show sampled pixels in source image.
            %
            narginchk(2,4);

            if nargin <= 2
                imageSize = self.size;
            end
            
            if nargin <= 3
                showSamples = false;
            end
            
            % Obtain number pixels according to type of image
            if isa(self,'lightfield.image.ViewpointImage')
                numberPixels_u = self.size.numberMicrolenses_k;
                numberPixels_v = self.size.numberMicrolenses_l;
            elseif isa(self,'lightfield.image.MicrolensImage')
                numberPixels_u = self.size.numberPixels_i;
                numberPixels_v = self.size.numberPixels_j;
            elseif isa(self,'image.epi.EpipolarPlaneImage')
                numberPixels_u = self.size.numberPixels;
                numberPixels_v = self.size.numberFrames;
            else
                numberPixels_u = self.size.numberPixels_u;
                numberPixels_v = self.size.numberPixels_v;
            end
            
            % Convert to homography object
            homography = utils.enums.Classes.HOMOGRAPHY().convert(homography);
            imageSize  = utils.enums.Classes.IMAGE_SIZE().convert(imageSize);
            
            % Extrapolation values are nan
            EXTRAPOLATION_VALUE = nan;
            
            % Obtain source pixel coordinates
            [indices_u,indices_v] = ndgrid( (1:numberPixels_u) - 0.5 ...
                                          , (1:numberPixels_v) - 0.5 );

            % Obtain target pixel coordinates
            [sampleIndices_u,sampleIndices_v] = ndgrid( (1:imageSize.numberPixels_u ) - 0.5 ...
                                                      , (1:imageSize.numberPixels_v) - 0.5 );
            samplePixels = image.Pixel([sampleIndices_u(:),sampleIndices_v(:)]',false);
            samplePixels = samplePixels.setHomogeneousCoordinates();

            % Obtain sample pixel coordinates in source image and remove
            % homogeneous coordinates
            samplePixels.data = mldivide(homography.homographyMatrix,samplePixels.data);
            samplePixels      = samplePixels.removeHomogeneousCoordinates();

            % Transform pixel coordinates to image like structure assuming matlab
            % format
            newIndices_u = reshape(samplePixels.u,[1,imageSize.numberPixels_u ...
                                                    ,imageSize.numberPixels_v]);
            newIndices_v = reshape(samplePixels.v,[1,imageSize.numberPixels_u ...
                                                    ,imageSize.numberPixels_v]);
            
            if showSamples == true
                figure; self.show;
                hold on; plot(newIndices_u(:),newIndices_v(:),'r.');
            end

            % Obtain intensity image after applying homography
            % Put channel information in the last dimension
            imageData = zeros(imageSize.data(:)');
            for iChannel = 1:self.size.numberChannels
                imageData(iChannel,:,:) = interpn( indices_u, indices_v ...
                                                 , permute(self.data(iChannel,:,:),[2,3,1]) ... 
                                                 , newIndices_u, newIndices_v ...
                                                 , image.enums.InterpolationMethods.CUBIC().method ...
                                                 , EXTRAPOLATION_VALUE );
            end
            self.data = imageData;
        end        
        
        function self = interpolate(self,newSize,interpolationMethod)
            %
            % Interpolate image information in order to avoid aliasing.
            % This is useful when filtering results are being multiplied.
            %
            % INPUTS:
            %   1. newSize - new size to be considered for image. The
            %   default size is 2 times the image size.
            %   2. interpolationMethod - interpolation method. For more
            %   information see InterpolationMethods enum or see the
            %   imresize help. The default interpolation method is bicubic.
            %
            narginchk(1,3);
            
            % Define default interpolation method
            if nargin <= 2
                interpolationMethod = image.enums.InterpolationMethods.BICUBIC;
            end
            
            % Define default size
            if nargin <= 1
                newSize = self.size.data * 2;
                newSize(image.enums.Formats.INTERNAL_FORMAT().channelDimension) = ...
                          self.size.numberChannels;
            end
            newSize = utils.enums.Classes.IMAGE_SIZE().convert(newSize);

            % If image data is empty, throw error
            if isempty(self.data)
                error = MException( 'TemplateImage:interpolate:noData' ...
                                  , 'Image data is not defined...' );
                error.throw();
            end
            
            % Encode to matlab format
            matlabImage = self.MATLAB_FORMAT.encode(self.data);
            
            % Interpolate image. imresize needs matlab format to resize
            % color images. Decode matlab format to internal format
            matlabImage = imresize( matlabImage ...
                                  , [newSize.numberPixels_v, newSize.numberPixels_u] ...
                                  , interpolationMethod.method );
            self.data   = self.MATLAB_FORMAT.decode(matlabImage);
        end
        
        function show(self,titleText)
            %
            % Show image data.
            %
            % INPUTS:
            %   1. titleText - title to be displayed with image.
            %
            narginchk(1,2);
            
            if nargin <= 1
                titleText = '';
            end
            
            % If there is not data, throw error
            if isempty(self.data)
                error = MException( 'TemplateImage:show:noData' ...
                                  , 'Image data is not defined...' );
                error.throw();
            end
            
            imshow(self.MATLAB_FORMAT.encode(self.data));
            
            % Display title if it is given as input
            if ~isempty(titleText)
                title(titleText);
            else
                set(gca,'position',[0 0 1 1],'units','normalized')
            end
            axis tight;
        end
        
        function write(self,outputFilepath)
            %
            % Write image data to file.
            %
            % INPUTS:
            %   1. outputFilepath - local filepath to write image data.
            %
            narginchk(1,2);
            
            if nargin <= 1
                outputFilepath = 'image_output.png';
            end
            
            % Obtain fileparts
            [filepath,filename,extension] = fileparts(outputFilepath);
            % If extension is not provided consider png extension.
            if isempty(extension)
                extension = '.png';
            end
            
            % If image data is empty, throw error
            if isempty(self.data)
                error = MException( 'TemplateImage:write:noData' ...
                                  , 'Image data is not defined...' );
                error.throw();
            end
            
            % If output filename corresponds to eps format
            if strcmp(extension,'.eps') == true
                % Supress figure
                set(gcf,'visible','off');
                
                % Show image
                self.show()
                
                % Set axis resolution based on image and avoid printing
                % axis
                axis image;
                axis off;
                
                % Save figure to eps
                options = struct('norepos',[]);
                utils.eps.save2eps([filepath, filesep, filename, extension],options);
                close(gcf);
                
            % For other file formats. The data for imwrite, must be 
            % formatted in matlab format.
            else
                imwrite(self.MATLAB_FORMAT.encode(self.data), [filepath, filesep, filename, extension]);
            end
        end
        
        function self = applyMedianFilter(self,medianFilterType,maskSize)
            %
            % Apply median filter to image. 
            % 
            % INPUTS:
            %   1. medianFilterType - median filter type. Default is 2D
            %   median filter.
            %   2. maskSize - mask size. If median filter is 3D, the mask 
            %   size must have 3 dimensions. The default is a 3 x 3 for 2D
            %   median filters and 3 x 3 x 3 for 3D median filters.
            %
            narginchk(1,3);
            
            % Default median filter
            if nargin <= 1
                medianFilterType = image.enums.MedianFilters.MEDIAN_2D;
            end
            
            if nargin <= 2
                if medianFilterType == image.enums.MedianFilters.MEDIAN_2D
                    maskSize = [1,3,3];
                else
                    maskSize = [3,3,3];
                end
            end

            % Convert mask size to image size
            maskSize = utils.enums.Classes.IMAGE_SIZE().convert(maskSize);
            
            % Transform image and mask size to matlab image format
            matlabImage    = self.MATLAB_FORMAT.encode(self.data);
            matlabMaskSize = maskSize.data(self.MATLAB_FORMAT.encodingOrder);
            
            % Apply median filter
            if medianFilterType == image.enums.MedianFilters.MEDIAN_2D
                % The mask must have two dimensions
                for iChannel = 1:self.size.numberChannels
                    matlabImage(:,:,iChannel) = medfilt2(matlabImage(:,:,iChannel),matlabMaskSize(1:2));
                end
            else
                matlabImage = medfilt3(matlabImage,matlabMaskSize);
            end
            
            % Update image data and convert back to internal format
            self.data = self.MATLAB_FORMAT.decode(matlabImage);
        end
        
        function self = smooth(self,smoothingFilter)
            %
            % Apply smoothing to image. 
            % 
            % INPUTS:
            %   1. smoothingFilter - filter to be used for smoothing image.
            %   The default filter is a 5 x 5 gaussian mask.
            %
            narginchk(1,2);
            
            % Default filter mask
            if nargin <= 1
                smoothingFilter = image.filters.Gaussian([5,5],1).gaussian;
            end
            
            % Convert to filter instance
            smoothingFilter = utils.enums.Classes.FILTER().convert(smoothingFilter);
            
            % Smooth image
            smoothedImage   = smoothingFilter.apply(self);
            self.data       = smoothedImage.data;
        end
    end
    
    methods (Static)
        function self = ImageFromFile(imageFilepath,convertToDouble)
            %
            % Create image instance from image file.
            %
            % INPUTS:
            %   1. imageFilepath   - image filepath.
            %   2. convertToDouble - flag to indicate if image should be
            %   converted to double. Default is true.
            %
            narginchk(1,2);
            
            if nargin <= 1
                convertToDouble = true;
            end
            
            % Obtain image data in matlab format
            if convertToDouble == true
                matlabImage = im2double(imread(imageFilepath));
            else
                matlabImage = imread(imageFilepath);
            end
            
            % Convert to internal format and create image instance
            self = image.Image(image.Image.MATLAB_FORMAT.decode(matlabImage));
        end
    end
end