classdef Lightfield
    %LIGHTFIELD
    %   Lightfield methods and properties.
    
    properties (Constant)
        TOOLBOX_FORMAT = lightfield.enums.Formats.DANSEREAU_FORMAT()
            % Toolbox format that allows to decode and encode the lightfield data
    end
    
    properties
        data = []           % Lightfield data
    end
    
    properties (Dependent)
        size                % Lightfield size
        vectorize           % Vectorize lightfield
    end
    
    methods
        function self = Lightfield(varargin)
            %
            % Create lightfield instance.
            %
            % INPUTS:
            %   1. lightfieldData - lightfield data values.
            %
            if ~isempty(varargin)
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
        
        function self = set.data(self,newLightfieldData)
            self.data = newLightfieldData;
        end
        
        function lightfieldSize = get.size(self)
            if isempty(self.data)
                lightfieldSize = lightfield.LightfieldSize();
            else
                lightfieldSize = lightfield.LightfieldSize(size(self.data));
            end
        end
        
        function vector = get.vectorize(self)
            % Vectorize lightfield. Lightfield has the size of an image
            % with:
            %               channels x viewpoints x microlenses.
            if isempty(self.data)
                vector = image.Image();
            else
                vector = image.Image( reshape( self.data, [ self.size.numberChannels ...
                                                          , self.size.numberPixels_i * self.size.numberPixels_j ...
                                                          , self.size.numberMicrolenses_k * self.size.numberMicrolenses_l ] ));
            end
        end
    end
    
    methods
        function self = rgb2gray(self)
            %
            % Convert RGB values to grayscale:
            %       0.2990 * R + 0.5870 * G + 0.1140 * B 
            %
            % The rgb2gray matlab method can only be used for images.
            %
            
            % If no data is provided, return empty list
            if isempty(self.data)
                self.data = [];
            elseif self.size.numberChannels > 1
                self.data = 0.299 * self.data(1,:,:,:,:) ...
                          + 0.587 * self.data(2,:,:,:,:) ...
                          + 0.114 * self.data(3,:,:,:,:);
            end
        end
        
        function self = equalizeHistogram(self)
            %
            % Adjust the contrast of the lightfield based on the histogram.
            %
            
            % If no data is provided, return empty list
            if isempty(self.data)
                self.data = [];
            else
                % Equalize histogram and obtain lightfield in Dansereau's
                % structure. This is a toolbox function that requires
                % Dansereau's format for lightfield.
                lightfieldData = LFHistEqualize(self.TOOLBOX_FORMAT.encode(self.data));

                % Convert structure back to the structure assumed 
                % internally
                self.data = self.TOOLBOX_FORMAT.decode(lightfieldData);
            end
        end
        
        function self = interpolate(self,newSize,interpolationMethod,imageType)
            %
            % Interpolate lightfield information.
            %
            % INPUTS:
            %   1. newSize - new size to be considered for lightfield. The
            %   default size is 2 times the lightfield size.
            %   2. interpolationMethod - interpolation method. For more
            %   information see InterpolationMethods enum or see the
            %   interpn help. The default interpolation method is linear.
            %   3. imageType - image type to consider the dimensions to
            %   perform interpolation. This allows to choose between
            %   viewpoints and microlenses and epipolar plane images.
            %   Default is epipolar plane images.
            %
            narginchk(1,4);
            
            % Define default interpolation method
            if nargin <= 2
                interpolationMethod = image.enums.InterpolationMethods.LINEAR;
            end
            
            if nargin <= 3
                imageType = lightfield.epi.EpipolarPlaneImage();
            end
            
            % Define default size
            if nargin <= 1
                newSize = self.size.data * 2;
                newSize(lightfield.enums.Formats.INTERNAL_FORMAT().channelDimension) = ...
                          self.size.numberChannels;
            end
            newSize = utils.enums.Classes.LIGHTFIELD_SIZE().convert(newSize);

            % If lightfield data is empty, throw error
            if isempty(self.data)
                error = MException( 'Lightfield:interpolate:noData' ...
                                  , 'Lightfield data is not defined...' );
                error.throw();
            end
            
            % For interpolation methods specific to images
            if interpolationMethod == image.enums.InterpolationMethods.BILINEAR ...
            || interpolationMethod == image.enums.InterpolationMethods.BICUBIC ...
            || interpolationMethod == image.enums.InterpolationMethods.BOX ...
            || interpolationMethod == image.enums.InterpolationMethods.LANCZOS2 ...
            || interpolationMethod == image.enums.InterpolationMethods.LANCZOS3
                if isa(imageType,'lightfield.epi.EpipolarPlaneImage')
                    self = self.interpolateUsingEpipolarPlaneImages(newSize,interpolationMethod);
                else
                    self = self.interpolateUsingViewpointsAndMicrolenses(newSize,interpolationMethod);
                end
                
            else
                self = self.interpolateUsingAllLightfield(newSize,interpolationMethod);
            end
        end
        
        function self = interpolateUsingAllLightfield(self,newSize,interpolationMethod)
            %
            % Interpolate lightfield information using all lightfield
            % coordinates.
            %
            % INPUTS:
            %   1. newSize - new size to be considered for lightfield. The
            %   default size is 2 times the lightfield size.
            %   2. interpolationMethod - interpolation method. For more
            %   information see InterpolationMethods enum or see the
            %   interpn help. The default interpolation method is linear.
            %
            narginchk(1,3);
            
            % Define default interpolation method
            if nargin <= 2
                interpolationMethod = image.enums.InterpolationMethods.LINEAR;
            end
            
            % Define default size
            if nargin <= 1
                newSize = self.size.data * 2;
                newSize(lightfield.enums.Formats.INTERNAL_FORMAT().channelDimension) = ...
                          self.size.numberChannels;
            end
            newSize = utils.enums.Classes.LIGHTFIELD_SIZE().convert(newSize);

            % If lightfield data is empty, throw error
            if isempty(self.data)
                error = MException( 'Lightfield:interpolate:noData' ...
                                  , 'Lightfield data is not defined...' );
                error.throw();
            end
            
            % Obtain original lightfield coordinates
            [pixels_i,pixels_j,microlenses_k,microlenses_l] = ...
                            ndgrid( 1:self.size.numberPixels_i ...
                                  , 1:self.size.numberPixels_j ...
                                  , 1:self.size.numberMicrolenses_k ...
                                  , 1:self.size.numberMicrolenses_l );

            % Obtain step size for new lightfield sampling
            pixelsStep_i      = (self.size.numberPixels_i - 1) / (newSize.numberPixels_i - 1);
            pixelsStep_j      = (self.size.numberPixels_j - 1) / (newSize.numberPixels_j - 1);
            microlensesStep_k = (self.size.numberMicrolenses_k - 1) / (newSize.numberMicrolenses_k- 1);
            microlensesStep_l = (self.size.numberMicrolenses_l - 1) / (newSize.numberMicrolenses_l - 1);

            % Obtain sampling for lightfield coordinates
            [newPixels_i,newPixels_j,newMicrolenses_k,newMicrolenses_l] = ...
                            ndgrid( 1:pixelsStep_i:self.size.numberPixels_i ...
                                  , 1:pixelsStep_j:self.size.numberPixels_j ...
                                  , 1:microlensesStep_k:self.size.numberMicrolenses_k ...
                                  , 1:microlensesStep_l:self.size.numberMicrolenses_l );

            % Interpolate for each channel
            newLightfieldData = zeros(newSize.data);
            for iChannel = 1:self.size.numberChannels
                % Obtain channel information
                colorLightfield = permute(self.data(iChannel,:,:,:,:),[2,3,4,5,1]);

                % Obtain interpolated color channel
                colorLightfield = interpn( pixels_i, pixels_j, microlenses_k, microlenses_l, colorLightfield ...
                                         , newPixels_i, newPixels_j, newMicrolenses_k, newMicrolenses_l ...
                                         , interpolationMethod.method );

                % Set interpolated color lightfield
                newLightfieldData(iChannel,:,:,:,:) = colorLightfield;
            end
            self.data = newLightfieldData;
        end
        
        function self = interpolateUsingEpipolarPlaneImages(self,newSize,interpolationMethod)
            %
            % Interpolate lightfield information using epipolar plane
            % images.
            %
            % INPUTS:
            %   1. newSize - new size to be considered for lightfield. The
            %   default size is 2 times the lightfield size.
            %   2. interpolationMethod - interpolation method. For more
            %   information see InterpolationMethods enum or see the
            %   interpn help. The default interpolation method is linear.
            %
            narginchk(1,3);
            
            % Define default interpolation method
            if nargin <= 2
                interpolationMethod = image.enums.InterpolationMethods.LINEAR;
            end
            
            % Define default size
            if nargin <= 1
                newSize = self.size.data * 2;
                newSize(lightfield.enums.Formats.INTERNAL_FORMAT().channelDimension) = ...
                          self.size.numberChannels;
            end
            newSize = utils.enums.Classes.LIGHTFIELD_SIZE().convert(newSize);

            % If lightfield data is empty, throw error
            if isempty(self.data)
                error = MException( 'Lightfield:interpolate:noData' ...
                                  , 'Lightfield data is not defined...' );
                error.throw();
            end
            
            % Obtain new size for (i,k) epipolar plane images if sizes have
            % changed.
            newLightfieldData = self.data;
            if newSize.numberPixels_i      ~= self.size.numberPixels_i ...
            || newSize.numberMicrolenses_k ~= self.size.numberMicrolenses_k 
                % Modify dimensions to have (i,k) epipolar plane image 
                % format in matlab image format
                newLightfieldData     = lightfield.epi.enums.EpipolarPlaneImageTypes.EPI_IK().encode(newLightfieldData);
                lightfieldData_matlab = image.Image.MATLAB_FORMAT.encode(newLightfieldData);

                % Obtain new size for (i,k) epipolar plane images
                lightfieldData_matlab = imresize( lightfieldData_matlab ...
                                                , [newSize.numberPixels_i, newSize.numberMicrolenses_k] ...
                                                , interpolationMethod.method );

                % Modify to internal format and to lightfield format
                newLightfieldData = image.Image.MATLAB_FORMAT.decode(lightfieldData_matlab);
                newLightfieldData = lightfield.epi.enums.EpipolarPlaneImageTypes.EPI_IK().decode(newLightfieldData);
            end

            % Obtain new size for (j,l) epipolar plane images if sizes have
            % changed.
            if newSize.numberPixels_j      ~= self.size.numberPixels_j ...
            || newSize.numberMicrolenses_l ~= self.size.numberMicrolenses_l
                % Modify dimensions to have (j,l) epipolar plane image 
                % format in matlab image format
                newLightfieldData     = lightfield.epi.enums.EpipolarPlaneImageTypes.EPI_JL().encode(newLightfieldData);
                lightfieldData_matlab = image.Image.MATLAB_FORMAT.encode(newLightfieldData);

                % Obtain new size for (j,l) epipolar plane images
                lightfieldData_matlab = imresize( lightfieldData_matlab ...
                                                , [newSize.numberPixels_j, newSize.numberMicrolenses_l] ...
                                                , interpolationMethod.method );

                % Modify to internal format and to lightfield format
                newLightfieldData = image.Image.MATLAB_FORMAT.decode(lightfieldData_matlab);
                newLightfieldData = lightfield.epi.enums.EpipolarPlaneImageTypes.EPI_JL().decode(newLightfieldData);
            end
            self.data = newLightfieldData;
        end
        
        function self = interpolateUsingViewpointsAndMicrolenses(self,newSize,interpolationMethod)
            %
            % Interpolate lightfield information using viewpoints and
            % microlenses images.
            %
            % INPUTS:
            %   1. newSize - new size to be considered for lightfield. The
            %   default size is 2 times the lightfield size.
            %   2. interpolationMethod - interpolation method. For more
            %   information see InterpolationMethods enum or see the
            %   interpn help. The default interpolation method is linear.
            %
            narginchk(1,3);
            
            % Define default interpolation method
            if nargin <= 2
                interpolationMethod = image.enums.InterpolationMethods.LINEAR;
            end
            
            % Define default size
            if nargin <= 1
                newSize = self.size.data * 2;
                newSize(lightfield.enums.Formats.INTERNAL_FORMAT().channelDimension) = ...
                          self.size.numberChannels;
            end
            newSize = utils.enums.Classes.LIGHTFIELD_SIZE().convert(newSize);

            % If lightfield data is empty, throw error
            if isempty(self.data)
                error = MException( 'Lightfield:interpolate:noData' ...
                                  , 'Lightfield data is not defined...' );
                error.throw();
            end
            
            % Obtain new size for microlenses images if sizes have
            % changed.
            newLightfieldData = self.data;
            if newSize.numberPixels_i ~= self.size.numberPixels_i ...
            || newSize.numberPixels_j ~= self.size.numberPixels_j
                % Modify dimensions to have microlens format in matlab
                % image format
                newLightfieldData     = permute(newLightfieldData,lightfield.image.MicrolensImage.ENCODING_ORDER);
                lightfieldData_matlab = lightfield.image.MicrolensImage.MATLAB_FORMAT.encode(newLightfieldData);

                % Obtain new size for microlenses images
                lightfieldData_matlab = imresize( lightfieldData_matlab ...
                                                , [newSize.numberPixels_j, newSize.numberPixels_i] ...
                                                , interpolationMethod.method );

                % Modify to internal format and to lightfield format
                newLightfieldData = lightfield.image.MicrolensImage.MATLAB_FORMAT.decode(lightfieldData_matlab);
                newLightfieldData = permute(newLightfieldData,lightfield.image.MicrolensImage.DECODING_ORDER);
            end

            % Obtain new size for viewpoint images if sizes have
            % changed.
            if newSize.numberMicrolenses_k ~= self.size.numberMicrolenses_k ...
            || newSize.numberMicrolenses_l ~= self.size.numberMicrolenses_l 
                % Modify dimensions to have viewpoint format in matlab
                % image format
                newLightfieldData     = permute(newLightfieldData,lightfield.image.ViewpointImage.ENCODING_ORDER);
                lightfieldData_matlab = lightfield.image.ViewpointImage.MATLAB_FORMAT.encode(newLightfieldData);

                % Obtain new size for viewpoint images
                lightfieldData_matlab = imresize( lightfieldData_matlab ...
                                                , [newSize.numberMicrolenses_l, newSize.numberMicrolenses_k] ...
                                                , interpolationMethod.method );

                % Modify to internal format and to lightfield format
                newLightfieldData = lightfield.image.ViewpointImage.MATLAB_FORMAT.decode(lightfieldData_matlab);
                newLightfieldData = permute(newLightfieldData,lightfield.image.ViewpointImage.DECODING_ORDER);
            end
            self.data = newLightfieldData;
        end
        
        function self = shear( self, dk_di, shearMethod, referenceViewpoint ...
                             , memoryEfficient )
            %
            % Shear lightfield. This allows to obtain a new lightfield by
            % sampling the coordinates of the original lightfield on the
            % microlenses. For a more general formulation of the shearing
            % please refer to SurfaceCameraImage object.
            %
            % INPUTS:
            %   1. dk_di - disparity change between consecutive viewpoint
            %   images. This slope relates with the depth of the virtual 
            %   focal plane.
            %   2. shearMethod - shear method. For more information see
            %   ShearMethods enum or see the interpn help. The default 
            %   interpolation method is frequency.
            %   3. referenceViewpoint - reference viewpoint considered to
            %   perform shearing. The default corresponds to the central
            %   viewpoint.
            %   4. memoryEfficient  - shearing using memory efficient
            %   implementation.
            %
            narginchk(1,5);
            
            % If disparity slope is not provided, consider the virtual
            % focal plane to be the real world focal plane.
            if nargin <= 1
                dk_di = 0;
            end
            
            % Define default shear method
            if nargin <= 2
                shearMethod = lightfield.enums.ShearMethods.FREQUENCY;
            end
            
            % Define default reference viewpoint. No need to transform
            % reference viewpoint to image ray, since we are not using it
            % in the method.
            if nargin <= 3
                referenceViewpoint = [ round(self.size.numberPixels_i/2) ...
                                     , round(self.size.numberPixels_j/2) ...
                                     , nan, nan];
            end
            
            if nargin <= 4
                memoryEfficient = true;
            end
            
            % Shear method using memory efficient Fourier transform
            if memoryEfficient == true ...
            && shearMethod == lightfield.enums.ShearMethods.FREQUENCY 
                self = self.shearUsingMemoryEfficientFrequencyOnViewpoints(dk_di,referenceViewpoint);

            % Shear method using Fourier transform
            elseif shearMethod == lightfield.enums.ShearMethods.FREQUENCY
                self = self.shearUsingFrequencyOnViewpoints(dk_di,referenceViewpoint);
                
            % Other shear methods using interpolation
            else
                self = self.shearUsingSpatialDomain(dk_di,shearMethod,referenceViewpoint);
            end
        end
        
        function self = shearUsingFrequencyOnViewpoints(self,dk_di,referenceViewpoint)
            %
            % Shear lightfield using frequency domain. This allows to 
            % obtain a new lightfield by sampling the coordinates of the 
            % original lightfield on the microlenses. For a more general 
            % formulation of the shearing please refer to 
            % SurfaceCameraImage object.
            %
            % INPUTS:
            %   1. dk_di - disparity change between consecutive viewpoint
            %   images. This slope relates with the depth of the virtual 
            %   focal plane.
            %   2. referenceViewpoint - reference viewpoint considered to
            %   perform shearing. The default corresponds to the central
            %   pixel.
            %
            narginchk(1,3);
            
            % If disparity slope is not provided, consider the virtual
            % focal plane to be the real world focal plane.
            if nargin <= 1
                dk_di = 0;
            end
            
            % Define default reference viewpoint
            if nargin <= 2
                referenceViewpoint = [ round(self.size.numberPixels_i/2) ...
                                     , round(self.size.numberPixels_j/2) ...
                                     , nan, nan];
            end
            
            % Transform reference viewpoint to image ray instance
            referenceViewpoint = utils.enums.Classes.IMAGE_RAY().convert(referenceViewpoint);
            
            % Obtain lightfield coordinates
            [pixels_i,pixels_j,microlenses_k,microlenses_l] = ...
                            ndgrid( 1:self.size.numberPixels_i ...
                                  , 1:self.size.numberPixels_j ...
                                  , 1:self.size.numberMicrolenses_k ...
                                  , 1:self.size.numberMicrolenses_l );
            coordinates = lightfield.ImageRay( [ pixels_i(:), pixels_j(:) ...
                                               , microlenses_k(:), microlenses_l(:) ]' ...
                                             , false );
                                    
            % Shearing is only performed on the microlenses.
            % The following steps are done in an unique command to prevent
            % memory overflow.
            % Compute discrete fourier transform of viewpoint images in
            % matlab format. For that, transform lightfield internal format
            % to viewpoint images format and then to matlab format.
            frequencyViewpointImagesMatlab = fft2( ...                                    % Perform DFT
                lightfield.image.ViewpointImage.MATLAB_FORMAT.encode( ...                 % Transform images to matlab format
                    permute(self.data,lightfield.image.ViewpointImage.ENCODING_ORDER) )); % Transform to viewpoint image

            % Obtain frequency associated with each position of the
            % discrete fourier transform
            frequencies_k = ifftshift( -fix(self.size.numberMicrolenses_k/2) ...
                                     : ceil(self.size.numberMicrolenses_k/2) - 1 );
            frequencies_l = ifftshift( -fix(self.size.numberMicrolenses_l/2) ...
                                     : ceil(self.size.numberMicrolenses_l/2) - 1 );
            [frequencies_k,frequencies_l] = meshgrid(frequencies_k,frequencies_l);

            % Obtain sub-pixel displacement for each viewpoint image. This
            % considers the internal format for the lightfield.
            delta_k = dk_di .* (coordinates.i - referenceViewpoint.i);
            delta_l = dk_di .* (coordinates.j - referenceViewpoint.j);

            % Obtain lightfield size not considering the channel 
            % information consider the internal format for the lightfield
            lightfieldSize = [ self.size.numberPixels_i ...
                             , self.size.numberPixels_j ...
                             , self.size.numberMicrolenses_k ...
                             , self.size.numberMicrolenses_l ];

            % Transform coordinates to lightfield like structure recovering
            % the channel information and convert to matlab format
            delta_k = permute(reshape(delta_k,lightfieldSize),[5,1,2,3,4]);
            delta_l = permute(reshape(delta_l,lightfieldSize),[5,1,2,3,4]);
            delta_k = permute(delta_k,lightfield.image.ViewpointImage.ENCODING_ORDER);
            delta_l = permute(delta_l,lightfield.image.ViewpointImage.ENCODING_ORDER);
            delta_k = lightfield.image.ViewpointImage.MATLAB_FORMAT.encode(delta_k);
            delta_l = lightfield.image.ViewpointImage.MATLAB_FORMAT.encode(delta_l);

            % Obtain size to replicate structure of frequencies. The 
            % replication is performed when performing the frequency 
            % product to avoid memory overflow.
            lightfieldSizeMatlab = self.size.data(lightfield.image.ViewpointImage.ENCODING_ORDER);
            lightfieldSizeMatlab = [ lightfieldSizeMatlab(image.enums.Formats.MATLAB_FORMAT.encodingOrder) ...
                                   , lightfieldSizeMatlab(length(image.enums.Formats.MATLAB_FORMAT.encodingOrder) + 1:end) ];
            lightfieldSizeMatlab(image.enums.Formats.MATLAB_FORMAT.pixelDimension_u) = 1;
            lightfieldSizeMatlab(image.enums.Formats.MATLAB_FORMAT.pixelDimension_v) = 1;

            % Shift viewpoint images on the frequency domain.
            % Considering that X(k) = F(x(n)), then x(n+m) has the Fourier
            % transform given by :
            %       X(k) * \exp^{(j 2 \pi m k/N)}
            frequencyViewpointImagesMatlab = frequencyViewpointImagesMatlab ...
                                          .* exp(1i * 2 * pi .* ...
                                                    ( delta_k .* repmat(frequencies_k./self.size.numberMicrolenses_k,lightfieldSizeMatlab) ...
                                                    + delta_l .* repmat(frequencies_l./self.size.numberMicrolenses_l,lightfieldSizeMatlab) ));

            % The following steps are done in an unique command to prevent
            % memory overflow.
            % Obtain viewpoint images again and transform to lightfield 
            % format.
            self.data = permute( lightfield.image.ViewpointImage.MATLAB_FORMAT.decode( ... % Transform to viewpoint format
                                          abs(ifft2(frequencyViewpointImagesMatlab)) ) ... % Perform IDFT
                               , lightfield.image.ViewpointImage.DECODING_ORDER );         % Transform to internal format
        end
        
        function self = shearUsingMemoryEfficientFrequencyOnViewpoints(self,dk_di,referenceViewpoint)
            %
            % Shear lightfield using frequency domain. This allows to 
            % obtain a new lightfield by sampling the coordinates of the 
            % original lightfield on the microlenses. The memory efficient 
            % procedure is obtained by performing 2D shearings. For a more 
            % general formulation of the shearing please refer to 
            % SurfaceCameraImage object.
            %
            % INPUTS:
            %   1. dk_di - disparity change between consecutive viewpoint
            %   images. This slope relates with the depth of the virtual 
            %   focal plane.
            %   2. referenceViewpoint - reference viewpoint considered to
            %   perform shearing. The default corresponds to the central
            %   pixel.
            %
            narginchk(1,3);
            
            % If disparity slope is not provided, consider the virtual
            % focal plane to be the real world focal plane.
            if nargin <= 1
                dk_di = 0;
            end
            
            % Define default reference viewpoint
            if nargin <= 2
                referenceViewpoint = [ round(self.size.numberPixels_i/2) ...
                                     , round(self.size.numberPixels_j/2) ...
                                     , nan, nan];
            end
            
            % Transform reference viewpoint to image ray instance
            referenceViewpoint = utils.enums.Classes.IMAGE_RAY().convert(referenceViewpoint);
            
            % Shearing is only performed on the microlenses.
            % Transform to viewpoint images and convert to matlab internal
            % format.
            viewpointImagesMatlab = lightfield.image.ViewpointImage.MATLAB_FORMAT.encode( ...               % Transform images to matlab format
                                        permute(self.data,lightfield.image.ViewpointImage.ENCODING_ORDER )); % Transform to viewpoint image

            % Obtain frequency associated with each position of the
            % discrete fourier transform
            frequencies_k = ifftshift( -fix(self.size.numberMicrolenses_k/2) ...
                                     : ceil(self.size.numberMicrolenses_k/2) - 1 );
            frequencies_l = ifftshift( -fix(self.size.numberMicrolenses_l/2) ...
                                     : ceil(self.size.numberMicrolenses_l/2) - 1 );
            [frequencies_k,frequencies_l] = meshgrid(frequencies_k,frequencies_l);

            % Obtain size to replicate structure of frequencies. The
            % replication is performed when performing the frequency 
            % product to avoid memory overflow.
            lightfieldSizeMatlab = [1,1,self.size.numberChannels];

            % Consider a 2D shearing to prevent memory overflow.
            for pixel_i = 1:self.size.numberPixels_i
                for pixel_j = 1:self.size.numberPixels_j
                    % Obtain sub-pixel displacement
                    delta_k = dk_di * (pixel_i - referenceViewpoint.i);
                    delta_l = dk_di * (pixel_j - referenceViewpoint.j);

                    % Compute discrete fourier transform
                    frequencyViewpointImageMatlab = fft2(viewpointImagesMatlab(:,:,:,pixel_i,pixel_j));

                    % Shift viewpoint images on the frequency domain.
                    % Considering that X(k) = F(x(n)), then x(n+m) has 
                    % the Fourier transform given by :
                    %       X(k) * \exp^{(j 2 \pi m k/N)}
                    frequencyViewpointImageMatlab = frequencyViewpointImageMatlab ...
                                                  .* exp(1i * 2 * pi .* ...
                                                            ( delta_k .* repmat(frequencies_k./self.size.numberMicrolenses_k,lightfieldSizeMatlab) ...
                                                            + delta_l .* repmat(frequencies_l./self.size.numberMicrolenses_l,lightfieldSizeMatlab) ));

                    % Obtain viewpoint images again and transform to 
                    % lightfield format.
                    self.data(:,pixel_i,pixel_j,:,:) = ...
                        permute( lightfield.image.ViewpointImage.MATLAB_FORMAT.decode( ... % Transform to image internal format
                                           abs(ifft2(frequencyViewpointImageMatlab)) ) ... % Perform IDFT
                               , lightfield.image.ViewpointImage.DECODING_ORDER );         % Transform to lightfield format
                end
            end
        end
        
        function self = shearUsingFrequencyOnMicrolenses(self,dk_di,referenceMicrolens)
            %
            % Shear lightfield using frequency domain. This allows to 
            % obtain a new lightfield by sampling the coordinates of the 
            % original lightfield on the viewpoints.
            %
            % INPUTS:
            %   1. dk_di - disparity change between consecutive viewpoint
            %   images. This slope relates with the depth of the virtual 
            %   focal plane.
            %   2. referenceMicrolens - reference microlens considered to
            %   perform shearing. The default corresponds to the central
            %   microlens.
            %
            narginchk(1,3);
            
            % If disparity slope is not provided, consider the virtual
            % focal plane to be the real world focal plane.
            if nargin <= 1
                dk_di = 0;
            end
            
            % Define default reference microlens
            if nargin <= 2
                referenceMicrolens = [ nan, nan ...
                                     , round(self.size.numberMicrolenses_k/2) ...
                                     , round(self.size.numberMicrolenses_l/2) ];
            end
            
            % Transform reference microlens to image ray instance
            referenceMicrolens = utils.enums.Classes.IMAGE_RAY().convert(referenceMicrolens);
            
            % Obtain lightfield coordinates
            [pixels_i,pixels_j,microlenses_k,microlenses_l] = ...
                            ndgrid( 1:self.size.numberPixels_i ...
                                  , 1:self.size.numberPixels_j ...
                                  , 1:self.size.numberMicrolenses_k ...
                                  , 1:self.size.numberMicrolenses_l );
            coordinates = lightfield.ImageRay( [ pixels_i(:), pixels_j(:) ...
                                               , microlenses_k(:), microlenses_l(:) ]' ...
                                             , false );
                                    
            % The following steps are done in an unique command to prevent
            % memory overflow.
            % Compute discrete fourier transform of microlens images in
            % matlab format. For that, transform lightfield internal format
            % to microlens images format and then to matlab format.
            frequencyMicrolensImagesMatlab = fft2( ...                                    % Perform DFT
                lightfield.image.MicrolensImage.MATLAB_FORMAT.encode( ...                 % Transform images to matlab format
                    permute(self.data,lightfield.image.MicrolensImage.ENCODING_ORDER) )); % Transform to microlens image

            % Obtain frequency associated with each position of the
            % discrete fourier transform
            frequencies_i = ifftshift( -fix(self.size.numberPixels_i/2) ...
                                     : ceil(self.size.numberPixels_i/2) - 1 );
            frequencies_j = ifftshift( -fix(self.size.numberPixels_j/2) ...
                                     : ceil(self.size.numberPixels_j/2) - 1 );
            [frequencies_i,frequencies_j] = meshgrid(frequencies_i,frequencies_j);

            % Obtain sub-pixel displacement for each microlens image. 
            % This considers the internal format for the lightfield.
            delta_i = (coordinates.k - referenceMicrolens.k) ./ dk_di;
            delta_j = (coordinates.l - referenceMicrolens.l) ./ dk_di;
                
            % Obtain lightfield size not considering the channel 
            % information consider the internal format for the lightfield
            lightfieldSize = [ self.size.numberPixels_i ...
                             , self.size.numberPixels_j ...
                             , self.size.numberMicrolenses_k ...
                             , self.size.numberMicrolenses_l ];

            % Transform coordinates to lightfield like structure recovering
            % the channel information and convert to matlab format
            delta_i = permute(reshape(delta_i,lightfieldSize),[5,1,2,3,4]);
            delta_j = permute(reshape(delta_j,lightfieldSize),[5,1,2,3,4]);
            delta_i = permute(delta_i,lightfield.image.MicrolensImage.ENCODING_ORDER);
            delta_j = permute(delta_j,lightfield.image.MicrolensImage.ENCODING_ORDER);
            delta_i = lightfield.image.MicrolensImage.MATLAB_FORMAT.encode(delta_i);
            delta_j = lightfield.image.MicrolensImage.MATLAB_FORMAT.encode(delta_j);
                
            % Obtain size to replicate structure of frequencies. The 
            % replication is performed when performing the frequency 
            % product to avoid memory overflow.
            lightfieldSizeMatlab = self.size.data(lightfield.image.MicrolensImage.ENCODING_ORDER);
            lightfieldSizeMatlab = [ lightfieldSizeMatlab(image.enums.Formats.MATLAB_FORMAT.encodingOrder) ...
                                   , lightfieldSizeMatlab(length(image.enums.Formats.MATLAB_FORMAT.encodingOrder) + 1:end) ];
            lightfieldSizeMatlab(image.enums.Formats.MATLAB_FORMAT.pixelDimension_u) = 1;
            lightfieldSizeMatlab(image.enums.Formats.MATLAB_FORMAT.pixelDimension_v) = 1;

            % Shift microlens images on the frequency domain.
            % Considering that X(k) = F(x(n)), then x(n+m) has the Fourier
            % transform given by :
            %       X(k) * \exp^{(j 2 \pi m k/N)}
            frequencyMicrolensImagesMatlab = frequencyMicrolensImagesMatlab ...
                                          .* exp(1i * 2 * pi .* ...
                                                    ( delta_i .* repmat(frequencies_i./lightfieldInstance.size.numberPixels_i,lightfieldSizeMatlab) ...
                                                    + delta_j .* repmat(frequencies_j./lightfieldInstance.size.numberPixels_j,lightfieldSizeMatlab) ));

            % The following steps are done in an unique command to prevent
            % memory overflow.
            % Obtain microlens images again and transform to lightfield 
            % format.
            self.data = permute( lightfield.image.MicrolensImage.MATLAB_FORMAT.decode( ... % Transform to microlens format
                                          abs(ifft2(frequencyMicrolensImagesMatlab)) ) ... % Perform IDFT
                               , lightfield.image.MicrolensImage.DECODING_ORDER );         % Transform to internal format
        end
        
        function self = shearUsingMemoryEfficientFrequencyOnMicrolenses(self,dk_di,referenceMicrolens)
            %
            % Shear lightfield using frequency domain. This allows to 
            % obtain a new lightfield by sampling the coordinates of the 
            % original lightfield on the viewpoints. The memory efficient 
            % procedure is obtained by performing 2D shearings.
            %
            % INPUTS:
            %   1. dk_di - disparity change between consecutive viewpoint
            %   images. This slope relates with the depth of the virtual 
            %   focal plane.
            %   2. referenceMicrolens - reference microlens considered to
            %   perform shearing. The default corresponds to the central
            %   microlens.
            %
            narginchk(1,3);
            
            % If disparity slope is not provided, consider the virtual
            % focal plane to be the real world focal plane.
            if nargin <= 1
                dk_di = 0;
            end
            
            % Define default reference microlens
            if nargin <= 2
                referenceMicrolens = [ nan, nan ...
                                     , round(self.size.numberMicrolenses_k/2) ...
                                     , round(self.size.numberMicrolenses_l/2) ];
            end
            
            % Transform reference microlens to image ray instance
            referenceMicrolens = utils.enums.Classes.IMAGE_RAY().convert(referenceMicrolens);
            
            % Transform to microlens images and convert to matlab internal
            % format.
            microlensImagesMatlab = lightfield.image.MicrolensImage.MATLAB_FORMAT.encode( ...                % Transform images to matlab format
                                        permute(self.data,lightfield.image.MicrolensImage.ENCODING_ORDER )); % Transform to microlens image

            % Obtain frequency associated with each position of the
            % discrete fourier transform
            frequencies_i = ifftshift( -fix(self.size.numberPixels_i/2) ...
                                     : ceil(self.size.numberPixels_i/2) - 1 );
            frequencies_j = ifftshift( -fix(self.size.numberPixels_j/2) ...
                                     : ceil(self.size.numberPixels_j/2) - 1 );
            [frequencies_i,frequencies_j] = meshgrid(frequencies_i,frequencies_j);

            % Obtain size to replicate structure of frequencies. The
            % replication is performed when performing the frequency 
            % product to avoid memory overflow.
            lightfieldSizeMatlab = [1,1,self.size.numberChannels];

            % Consider a 2D shearing to prevent memory overflow.
            for microlens_k = 1:self.size.numberMicrolenses_k
                for microlens_l = 1:self.size.numberMicrolenses_l
                    % Obtain sub-pixel displacement
                    delta_i = (microlens_k - referenceMicrolens.k) ./ dk_di;
                    delta_j = (microlens_l - referenceMicrolens.l) ./ dk_di;

                    % Compute discrete fourier transform
                    frequencyMicrolensImageMatlab = fft2(microlensImagesMatlab(:,:,:,microlens_k,microlens_l));

                    % Shift microlens images on the frequency domain.
                    % Considering that X(k) = F(x(n)), then x(n+m) has 
                    % the Fourier transform given by :
                    %       X(k) * \exp^{(j 2 \pi m k/N)}
                    frequencyMicrolensImageMatlab = frequencyMicrolensImageMatlab ...
                                                 .* exp(1i * 2 * pi .* ...
                                                           ( delta_i .* repmat(frequencies_i./self.size.numberPixels_i,lightfieldSizeMatlab) ...
                                                           + delta_j .* repmat(frequencies_j./self.size.numberPixels_j,lightfieldSizeMatlab) ));

                    % Obtain microlens images again and transform to 
                    % lightfield format.
                    self.data(:,:,:,microlens_k,microlens_k) = ...
                        permute( lightfield.image.MicrolensImage.MATLAB_FORMAT.decode( ... % Transform to image internal format
                                           abs(ifft2(frequencyMicrolensImageMatlab)) ) ... % Perform IDFT
                               , lightfield.image.MicrolensImage.DECODING_ORDER );         % Transform to lightfield format
                end
            end
        end
        
        function self = shearUsingSpatialDomain(self,dk_di,shearMethod,referenceViewpoint)
            %
            % Shear lightfield using spatial domain. This allows to obtain 
            % a new lightfield by sampling the coordinates of the original 
            % lightfield on the microlenses. For a more general formulation
            % of the shearing please refer to SurfaceCameraImage object.
            %
            % INPUTS:
            %   1. dk_di - disparity change between consecutive viewpoint
            %   images. This slope relates with the depth of the virtual 
            %   focal plane.
            %   2. shearMethod - shear method. For more information see
            %   ShearMethods enum or see the interpn help. The default 
            %   interpolation method is linear.
            %   3. referenceViewpoint - reference viewpoint considered to
            %   perform shearing. The default corresponds to the central
            %   pixel.
            %
            narginchk(1,4);
            
            % If disparity slope is not provided, consider the virtual
            % focal plane to be the real world focal plane.
            if nargin <= 1
                dk_di = 0;
            end
            
            % Define default shear method
            if nargin <= 2
                shearMethod = lightfield.enums.ShearMethods.LINEAR;
            end
            
            % Define default reference viewpoint
            if nargin <= 3
                referenceViewpoint = [ round(self.size.numberPixels_i/2) ...
                                     , round(self.size.numberPixels_j/2) ...
                                     , nan, nan];
            end
            
            EXTRAPOLATION_VALUE = nan;
            
            % Transform reference viewpoint to image ray instance
            referenceViewpoint = utils.enums.Classes.IMAGE_RAY().convert(referenceViewpoint);
            
            % Obtain lightfield coordinates
            [pixels_i,pixels_j,microlenses_k,microlenses_l] = ...
                            ndgrid( 1:self.size.numberPixels_i ...
                                  , 1:self.size.numberPixels_j ...
                                  , 1:self.size.numberMicrolenses_k ...
                                  , 1:self.size.numberMicrolenses_l );
            
            % Transform lightfield coordinates to vector
            coordinates = lightfield.ImageRay( [ pixels_i(:), pixels_j(:) ...
                                               , microlenses_k(:), microlenses_l(:) ]' ...
                                             , false );
                                             
            % Shearing on microlenses
            % Obtain new lightfield size
            newLightfieldSize = lightfield.LightfieldSize([1,size(pixels_i)]);

            % Obtain new lightfield considering the disparity on 
            % the viewpoints
            newMicrolenses_k = coordinates.k + dk_di .* (coordinates.i - referenceViewpoint.i);
            newMicrolenses_l = coordinates.l + dk_di .* (coordinates.j - referenceViewpoint.j);
            newPixels_i      = coordinates.i;
            newPixels_j      = coordinates.j;

            % Transform coordinates to lightfield like structure
            newPixels_i      = reshape(newPixels_i,newLightfieldSize.data);
            newPixels_j      = reshape(newPixels_j,newLightfieldSize.data);
            newMicrolenses_k = reshape(newMicrolenses_k,newLightfieldSize.data);
            newMicrolenses_l = reshape(newMicrolenses_l,newLightfieldSize.data);

            % If there is information in lightfield, perform shearing
            if ~isempty(self.data)
                % Add channel information to lightfield size
                newLightfieldSize.data(lightfield.enums.Formats.INTERNAL_FORMAT.channelDimension) = self.size.numberChannels;

                % Obtain sheared lightfield
                newLightfield = zeros(newLightfieldSize.data);
                for iChannel = 1:self.size.numberChannels
                    colorLightfield = permute(self.data(iChannel,:,:,:,:),[2,3,4,5,1]);
                    colorLightfield = interpn( pixels_i, pixels_j, microlenses_k, microlenses_l, colorLightfield ...
                                             , newPixels_i, newPixels_j, newMicrolenses_k, newMicrolenses_l ...
                                             , shearMethod.symbol, EXTRAPOLATION_VALUE );
                    newLightfield(iChannel,:,:,:,:) = colorLightfield;
                end
                self.data = newLightfield;
            end
        end
        
        function renderedImage = render(self,aperture,normalizeBorderPixels)
            %
            % Obtain image from rendering the entire lightfield.
            %
            % INPUTS: 
            %   1. aperture - define viewpoints to be added to the rendered
            %   image.
            %   2. normalizeBorderPixels - border pixels are not defined
            %   for every viewpoint image. Normalize the border pixels
            %   intensity by dividing by the number of valid rays captured
            %   by the lightfield.
            %
            narginchk(1,3);
            
            if nargin <= 1
                aperture = [1,self.size.numberPixels_i;1,self.size.numberPixels_j;nan,nan;nan,nan];
            end
            
            if nargin <= 2
                normalizeBorderPixels = true;
            end

            % Create image ray instance from aperture
            aperture = utils.enums.Classes.IMAGE_RAY().convert(aperture);
            
            if isempty(self.data)
                renderedImage = image.Image();
            else
                % Obtain aperture coordinates
                minimum_i = min(aperture.i);
                maximum_i = max(aperture.i);
                minimum_j = min(aperture.j);
                maximum_j = max(aperture.j);
                
                % Count the number of valid rays for each pixel. Only one
                % channel can be analyzed.
                if normalizeBorderPixels == true
                    apertureArea = ~isnan(self.data(1,minimum_i:maximum_i,minimum_j:maximum_j,:,:));
                    apertureArea = sum(sum(apertureArea, 2),3);
                else
                    apertureArea = (maximum_i - minimum_i + 1) * (maximum_j - minimum_j + 1);
                end
                
                % Integrate aperture images to render image. Do not
                % consider NaN values
                imageData = nansum( nansum( self.data(:,minimum_i:maximum_i,minimum_j:maximum_j,:,:), 2), 3) ...
                         ./ apertureArea;
                imageData = permute(imageData,lightfield.image.ViewpointImage.ENCODING_ORDER);
                
                % Obtain image instance
                renderedImage = image.Image(imageData);
            end
        end
        
        function self = show(self)
            %
            % Visualize the lightfield data.
            %
            
            % If no data is provided, return error
            if isempty(self.data)
                error = MException( 'Lightfield:show:noData' ...
                                  , 'Lightfield data not provided...' );
                error.throw();
            end
            
            try
                % This is a toolbox function that needs Dansereau's
                % lightfield format
                LFDispMousePan(self.TOOLBOX_FORMAT.encode(self.data),2); 
            catch 
                % Do nothing. The display function gives an error when the
                % the display is closed.
            end
        end
        
        function spiral_ij(self)
            %
            % Create gif presenting the images of the lightfield in spiral.
            %
            % Obtain number of cameras in each direction
            numberCameras_i = self.size.numberPixels_i;
            numberCameras_j = self.size.numberPixels_j;
            center_ij       = round([numberCameras_i;numberCameras_j]/2);
            
            % Reorganize LF to display a video from larger to shorter
            % distance
            [pixels_i,pixels_j]   = ndgrid(1:numberCameras_i,1:numberCameras_j);
            pixels_ij             = [pixels_i(:),pixels_j(:)]';
            distance_ij           = max(abs(pixels_ij - center_ij));
            [distance_ij,indices] = sort(distance_ij,'descend');
            uniqueDistances       = unique(distance_ij,'stable');
            
            % Intialize and fill video
            video  = zeros( numberCameras_i * numberCameras_j,self.size.numberChannels ...
                          , self.size.numberMicrolenses_k, self.size.numberMicrolenses_l );
            iFrame = 0;
            for iDistance = 1:length(uniqueDistances)
                % Obtain viewpoints with same distance
                distanceIndices = distance_ij == uniqueDistances(iDistance);
                cameras_ij      = pixels_ij(:,indices(distanceIndices));
                written         = false(1,size(cameras_ij,2));
                % Obtain minimum and maximum pixel i and j
                min_i = min(cameras_ij(1,:));
                max_i = max(cameras_ij(1,:));
                min_j = min(cameras_ij(2,:));
                max_j = max(cameras_ij(2,:));
                
                % Left to right frames (top)
                pixel_j = min_j;
                for pixel_i = min_i:max_i
                    % Validate if viewpoint exists
                    index = cameras_ij(1,:) == pixel_i & cameras_ij(2,:) == pixel_j;
                    if sum(index) > 0 && written(index) == false
                        % Obtain frame for viewpoint
                        iFrame              = iFrame + 1;
                        video(iFrame,:,:,:) = self.data(:,cameras_ij(1,index),cameras_ij(2,index),:,:);
                        written(index)      = true;
                    end
                end
                
                % Top to bottom frames (right)
                pixel_i = max_i;
                for pixel_j = min_j:max_j
                    % Validate if viewpoint exists
                    index = cameras_ij(1,:) == pixel_i & cameras_ij(2,:) == pixel_j;
                    if sum(index) > 0 && written(index) == false
                        % Obtain frame for viewpoint
                        iFrame              = iFrame + 1;
                        video(iFrame,:,:,:) = self.data(:,cameras_ij(1,index),cameras_ij(2,index),:,:);
                        written(index)      = true;
                    end
                end
                
                % Right to left frames (bottom)
                pixel_j = max_j;
                for pixel_i = max_i:-1:min_i
                    % Validate if viewpoint exists
                    index = cameras_ij(1,:) == pixel_i & cameras_ij(2,:) == pixel_j;
                    if sum(index) > 0 && written(index) == false
                        % Obtain frame for viewpoint
                        iFrame              = iFrame + 1;
                        video(iFrame,:,:,:) = self.data(:,cameras_ij(1,index),cameras_ij(2,index),:,:);
                        written(index)      = true;
                    end
                end
                
                % Bottom to top frames (left)
                pixel_i = min_i;
                for pixel_j = max_j:-1:min_j
                    % Validate if viewpoint exists
                    index = cameras_ij(1,:) == pixel_i & cameras_ij(2,:) == pixel_j;
                    if sum(index) > 0 && written(index) == false
                        % Obtain frame for viewpoint
                        iFrame              = iFrame + 1;
                        video(iFrame,:,:,:) = self.data(:,cameras_ij(1,index),cameras_ij(2,index),:,:);
                        written(index)      = true;
                    end
                end
                
                if all(written) == false
                    fprintf('Something is wrong...\n');
                end
            end
            
            % Save video as a gif
            filename = 'lf_spiral.gif';
            fig   = figure();
            axis tight manual;
            video = permute(video,[4,3,2,1]);
            for iFrame = 1:size(video,4)
                imshow(video(:,:,:,iFrame));
                
                % Capture the figure as an image
                frame      = getframe(fig);
                img        = frame2im(frame);
                [imind,cm] = rgb2ind(img,256);
                
                % Write to gif file
                if iFrame == 1
                    imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',0.05);
                else
                    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.05);
                end
            end
            close(fig);
        end
    end
    
    methods (Static)
        function self = LightfieldFromSurfaceCameraCollection(surfaceCameraCollection)
            % 
            % Create lightfield instance from surface camera collection.
            %
            % INPUTS:
            %   1. surfaceCameraCollection - surface camera collection
            %   instance.
            %
            narginchk(1,1);
            
            % Obtain lightfield size
            lightfieldSize = lightfield.LightfieldSize( [ surfaceCameraCollection.data{1,1}.size.numberChannels ...
                                                        , surfaceCameraCollection.data{1,1}.size.numberPixels_u ...
                                                        , surfaceCameraCollection.data{1,1}.size.numberPixels_v ...
                                                        , surfaceCameraCollection.size.numberMicrolenses_k ...
                                                        , surfaceCameraCollection.size.numberMicrolenses_l ] );

            % Initialize lightfield data
            lightfieldData = nan(lightfieldSize.data);
            for microlens_k = 1:lightfieldSize.numberMicrolenses_k
                for microlens_l = 1:lightfieldSize.numberMicrolenses_l
                    lightfieldData(:,:,:,microlens_k,microlens_l) = surfaceCameraCollection.data{microlens_k,microlens_l}.data;
                end
            end

            % Create lightfield instance
            self = lightfield.Lightfield(lightfieldData);
        end
        
        function self = LightfieldFromImageRays(imageRays_ijkl,lightfieldSize)
            %
            % Create lightfield instance from image rays considering the
            % lightfield size.
            %
            % INPUTS:
            %   1. imageRays_ijkl - image rays obtained from projection of
            %      rays in the image plane. These points should not be
            %      provided in homogeneous coordinates. Each image ray
            %      should be provided in different columns.
            %   2. lightfieldSize - size of the lightfield.
            %
            narginchk(2,2);

            % Create image rays structure and lightfield size
            imageRays_ijkl = utils.enums.Classes.IMAGE_RAY().convert(imageRays_ijkl);
            lightfieldSize = utils.enums.Classes.LIGHTFIELD_SIZE().convert(lightfieldSize);
            
            % Create black lightfield
            lightfieldData = zeros(lightfieldSize.data);
            
            % Round pixels and microlenses to the nearest integer
            imageRays_ijkl.data = round(imageRays_ijkl.data);
            
            % Obtain pixels and microlenses that are valid for lightfield
            % size
            validPixels_i       = imageRays_ijkl.i >= 1 & imageRays_ijkl.i <= lightfieldSize.numberPixels_i;
            validPixels_j       = imageRays_ijkl.j >= 1 & imageRays_ijkl.j <= lightfieldSize.numberPixels_j;
            validMicrolenses_k  = imageRays_ijkl.k >= 1 & imageRays_ijkl.k <= lightfieldSize.numberMicrolenses_k;
            validMicrolenses_l  = imageRays_ijkl.l >= 1 & imageRays_ijkl.l <= lightfieldSize.numberMicrolenses_l;
            imageRays_ijkl      = imageRays_ijkl.obtainVectors(  validPixels_i & validPixels_j ...
                                                               & validMicrolenses_k & validMicrolenses_l );

            % If number vectors obtained are greater or equal to one, fill
            % information on lightfield.
            if imageRays_ijkl.numberVectors > 0
                % Obtain linear indices for valid image rays considering the
                % number of colors in image
                linearIndices = [];
                for iChannel = 1:lightfieldSize.numberChannels
                    linearIndicesForColor = sub2ind( lightfieldSize.data ...
                                                   , iChannel * ones(1, imageRays_ijkl.numberVectors) ...
                                                   , imageRays_ijkl.i ...
                                                   , imageRays_ijkl.j ...
                                                   , imageRays_ijkl.k ...
                                                   , imageRays_ijkl.l );
                    linearIndices = cat(2,linearIndices,linearIndicesForColor);
                end

                % Fill lightfield with image rays obtained from projection and 
                % create lightfield instance
                lightfieldData(linearIndices) = 1;
            end
            
            % Create lightfield instance
            self = lightfield.Lightfield(lightfieldData);
        end
        
        function disparities = shear2disparity(shearValues)
            %
            % Obtain disparity on viewpoint images from shear values. The 
            % conversion is based on the shearing formula presented by Tao 
            % et al:
            %
            %       Tao, Michael W., et al. "Depth from combining defocus 
            %       and correspondence using light-field cameras." 
            %       Proceedings of the IEEE ICCV. 2013.
            %
            % In this formula, the new microlens coordinate k' is corrected
            % according to the current microlens k and pixel i coordinates: 
            %           
            %       k' = k + i (1 - 1 / \alpha)
            %
            % where \alpha is the shear value. Since the disparity on the
            % viewpoint images corresponds to dk/di, we have:
            %
            %       dk_di = 1 - 1 / \alpha
            %
            narginchk(1,1);
            
            disparities = 1 - 1./shearValues;
        end
        
        function shearValues = disparity2shear(disparities)
            %
            % Obtain shear values from disparities on viewpoint images. The 
            % conversion is based on the shearing formula presented by Tao 
            % et al:
            %
            %       Tao, Michael W., et al. "Depth from combining defocus 
            %       and correspondence using light-field cameras." 
            %       Proceedings of the IEEE ICCV. 2013.
            %
            % The shear value \alpha is given by:
            %
            %       \alpha = -1 / (dk_di - 1)
            %
            narginchk(1,1);
            
            shearValues = -1./(disparities - 1);
        end
    end    
end
