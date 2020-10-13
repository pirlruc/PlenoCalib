classdef Lightfield < lightfield.Lightfield
    %LIGHTFIELD
    %   Lightfield methods and properties with utilities regarding the
    %   decoding options used to obtain the lightfield.
    
    properties
        confidence            = []    % Confidence measurement associated with each item in lightfield
        decodeOptions         = []    % Decoding options used to obtain the lightfield
        lensletGridModel      = []    % Parameters of the grid model to decode the lightfield
        microlensCenters      = []    % Microlens centers on sensor according to lenslet grid model
        lightfieldMetadata    = []    % Lightfield metadata    
        rawLensletImage       = image.Image()   % Raw lenslet image instance
        lensletImage          = image.Image()   % Lenslet image after alignment
        rawWhiteImage         = image.Image()   % Raw white image instance
        whiteImage            = image.Image()   % White image after rescaling image values and alignment
        demosaickedLensletImage = image.Image()   % Lenslet image after rescaling image values and demosaicking
    end
    
    properties
        % The microlenses in the raw image have an hexagonal tilling. Thus,
        % the microlenses have two aligned grids. The microlenses in the
        % x-dimension are aligned.
        gridA       % Lightfield from microlenses aligned with first line 
        gridB       % Lightfield from microlenses aligned with second line
    end
    
    methods
        function self = Lightfield(varargin)
            %
            % Create decoded lightfield instance.
            %
            % INPUTS:
            %   01. lightfieldData - lightfield data values.
            %   02. confidence     - confidence measurements for each item
            %       in the lightfield data.
            %   03. decodeOptions  - decode options used to obtain the
            %       lightfield.
            %   04. lensletGridModel    - grid model data.
            %   05. microlensCenters    - microlens centers coordinates on
            %       image sensor.
            %   06. lightfieldMetadata  - lightfield metadata.
            %   07. rawLensletImageData - raw lenslet image data.
            %   08. demosaickedLensletImageData - lenslet image data after
            %       demosaicking and rescaling.
            %   09. lensletImageData  - lenslet image data after alignment.
            %   10. rawWhiteImageData - raw white image data.
            %   11. whiteImageData    - white image data after alignment.
            %
            narginchk(0,11);
            
            % Create superclass instance
            self = self@lightfield.Lightfield();
            
            if ~isempty(varargin)
                if nargin >= 11
                    self.whiteImage = varargin{11};
                end
                
                if nargin >= 10
                    self.rawWhiteImage = varargin{10};
                end
                
                if nargin >= 9
                    self.lensletImage = varargin{9};
                end
                
                if nargin >= 8
                    self.demosaickedLensletImage = varargin{8};
                end
                
                if nargin >= 7
                    self.rawLensletImage = varargin{7};
                end
                
                if nargin >= 6
                    self.lightfieldMetadata = varargin{6};
                end
                
                if nargin >= 5
                    self.microlensCenters = varargin{5};
                end
                
                if nargin >= 4
                    self.lensletGridModel = varargin{4};
                end
                
                if nargin >= 3
                    self.decodeOptions = varargin{3};
                end
                
                if nargin >= 2
                    self.confidence = varargin{2};
                end
                
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end

        function self = set.confidence(self,newConfidence)
            self.confidence = newConfidence;
        end
        
        function self = set.decodeOptions(self,newDecodeOptions)
            self.decodeOptions = newDecodeOptions;
        end
        
        function self = set.lensletGridModel(self,newLensletGridModel)
            self.lensletGridModel = newLensletGridModel;
        end
        
        function self = set.microlensCenters(self,newMicrolensCenters)
            self.microlensCenters = utils.enums.Classes.PIXEL().convert(newMicrolensCenters);
        end
        
        function self = set.lightfieldMetadata(self,newLightfieldMetadata)
            self.lightfieldMetadata = newLightfieldMetadata;
        end

        function self = set.rawLensletImage(self, newRawLensletImageData)
            self.rawLensletImage = utils.enums.Classes.IMAGE().convert(newRawLensletImageData);
        end
        
        function self = set.demosaickedLensletImage(self, newDemosaickedLensletImageData)
            self.demosaickedLensletImage = utils.enums.Classes.IMAGE().convert(newDemosaickedLensletImageData);
        end
        
        function self = set.lensletImage(self, newLensletImageData)
            self.lensletImage = utils.enums.Classes.IMAGE().convert(newLensletImageData);
        end
        
        function self = set.rawWhiteImage(self, newRawWhiteImageData)
            self.rawWhiteImage = utils.enums.Classes.IMAGE().convert(newRawWhiteImageData);
        end
        
        function self = set.whiteImage(self, newWhiteImageData)
            self.whiteImage = utils.enums.Classes.IMAGE().convert(newWhiteImageData);
        end
        
        function grid = get.gridA(self)
            grid = lightfield.lytro.Lightfield();
            if ~isempty(self.data)
                grid.confidence = self.confidence(:,:,:,:,1:2:end);
                grid.data       = self.data(:,:,:,:,1:2:end);
            end
        end
        
        function grid = get.gridB(self)
            grid = lightfield.lytro.Lightfield();
            if ~isempty(self.data)
                grid.confidence = self.confidence(:,:,:,:,2:2:end);
                grid.data       = self.data(:,:,:,:,2:2:end);
            end
        end
    end
    
    methods
        function self = correctColor(self,withWhiteBalancing)
            %
            % Applies a color correction matrix, balance vector, gamma, and
            % automatic white balancing for correcting the color of the 
            % lightfield.
            %
            % INPUTS:
            %   1. withWhiteBalancing - flag that indicates if an automatic
            %   white balancing should be performed for the lightfield.
            %   Default is false.
            %
            narginchk(1,2);
            
            if nargin <= 1
                withWhiteBalancing = false;
            end
            
            % If lightfield is empty, return empty list
            if isempty(self.data)
                self.data = [];
            
            else
                % This function deals with saturated input pixels by 
                % aggressively saturating output pixels and returns the
                % Dansereau's structure of lightfield. This function is
                % part of Dansereau's toolbox and needs that specific
                % format for lightfield.
                lightfieldData = LFColourCorrect( self.TOOLBOX_FORMAT.encode(self.data) ...
                                                , self.decodeOptions.ColourMatrix ...
                                                , self.decodeOptions.ColourBalance ...
                                                , self.decodeOptions.Gamma ...
                                                , withWhiteBalancing );

                % Convert structure back to the structure assumed 
                % internally
                self.data = self.TOOLBOX_FORMAT.decode(lightfieldData);
            end
        end
        
        function self = removeMicrolenses(self,numberMicrolenses)
            %
            % Remove microlenses from lightfield data.
            %
            % The decoding process of Dansereau gives poor results at the
            % borders of the viewpoint images. Thus, one might want to 
            % remove the viewpoint pixel borders (microlenses).
            %
            % INPUTS:
            %   1. numberMicrolenses - number of microlenses to remove from
            %   either direction of the viewpoint image.
            %
            narginchk(1,2);
            
            if nargin <= 1
                numberMicrolenses = 3;
            end
            
            % Concatenate lightfield and confidence measurement
            lightfieldData      = cat( lightfield.enums.Formats.INTERNAL_FORMAT().channelDimension ...
                                     , self.data, self.confidence );
            temporaryLightfield = lightfield.Lightfield(lightfieldData);

            % If lightfield is empty, return empty list
            if isempty(temporaryLightfield.data)
                self.data = [];
                self.confidence = [];
                
            % Otherwise, remove microlenses from lightfield. The
            % microlenses correspond to viewpoint pixels.
            else
                temporaryLightfield.data = temporaryLightfield.data(:,:,:, (1 + numberMicrolenses):(end - numberMicrolenses) ...
                                                                         , (1 + numberMicrolenses):(end - numberMicrolenses) );
                
                % Update lightfield
                self.data       = temporaryLightfield.data(1:self.size.numberChannels,:,:,:,:);
                
                % Update confidence data if there is data in the initial
                % instance.
                if ~isempty(self.confidence)
                    self.confidence = temporaryLightfield.data(end,:,:,:,:);
                end
            end
        end
        
        function self = removeViewpoints(self,numberViewpoints)
            %
            % Remove viewpoints from lightfield data.
            %
            % The decoding process of Dansereau obtains square images for
            % the microlenses. Nonetheless, the microlenses have a
            % hexagonal pattern. Thus, there are some dark pixels at the
            % borders of the microlenses images. This also results from 
            % vignetting. Thus, one might want to remove the microlens
            % pixel borders (viewpoints). 
            %
            % INPUTS:
            %   1. numberViewpoints - number of microlenses pixels to 
            %   remove from either direction of the microlens image.
            %
            narginchk(1,2);
            
            if nargin <= 1
                numberViewpoints = 1;
            end
            
            % Concatenate lightfield and confidence measurement
            lightfieldData      = cat( lightfield.enums.Formats.INTERNAL_FORMAT().channelDimension ...
                                     , self.data, self.confidence );
            temporaryLightfield = lightfield.Lightfield(lightfieldData);
            
            % If lightfield is empty, return empty list
            if isempty(temporaryLightfield.data)
                self.data = [];
                self.confidence = [];
            
            % Otherwise, remove viewpoints. This corresponds to pixels in
            % the microlens images.
            else
                temporaryLightfield.data = temporaryLightfield.data(:, (1 + numberViewpoints):(end - numberViewpoints) ...
                                                                     , (1 + numberViewpoints):(end - numberViewpoints),:,: );
                                                                 
                % Update lightfield
                self.data       = temporaryLightfield.data(1:self.size.numberChannels,:,:,:,:);
                
                % Update confidence data if there is data in the initial
                % instance.
                if ~isempty(self.confidence)
                    self.confidence = temporaryLightfield.data(end,:,:,:,:);
                end
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
            %   imresize help. The default interpolation method is linear.
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

            % Concatenate lightfield and confidence measurement
            lightfieldData      = cat( lightfield.enums.Formats.INTERNAL_FORMAT().channelDimension ...
                                     , self.data, self.confidence );
            temporaryLightfield = lightfield.Lightfield(lightfieldData);

            % Interpolate lightfield
            temporaryLightfield = temporaryLightfield.interpolate( newSize ...
                                                                 , interpolationMethod ...
                                                                 , imageType );

            % Update lightfield
            self.data       = temporaryLightfield.data(1:self.size.numberChannels,:,:,:,:);
            
            % Update confidence data if there is data in the initial
            % instance.
            if ~isempty(self.confidence)
                self.confidence = temporaryLightfield.data(end,:,:,:,:);
            end
        end
        
        function self = shear(self,dk_di,shearMethod,referenceViewpoint,memoryEfficient)
            %
            % Shear lightfield. This allows to obtain a new lightfield by
            % sampling the coordinates of the original lightfield.
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
            %   pixel.
            %   4. memoryEfficient  - shearing using memory efficient
            %   implementation.
            %
            narginchk(1,5);
            
            % If disparity slope is not provided, consider the virtual
            % focal plane to be the real world focal plane.
            if nargin <= 1
                dk_di = 0;
            end
            
            % Define default interpolation method
            if nargin <= 2
                shearMethod = lightfield.enums.ShearMethods.FREQUENCY;
            end
            
            % Define default reference viewpoint
            if nargin <= 3
                referenceViewpoint = [ round(self.size.numberPixels_i/2) ...
                                     , round(self.size.numberPixels_j/2) ...
                                     , nan, nan];
            end
            
            if nargin <= 4
                memoryEfficient = true;
            end
            
            % If data is not defined, do nothing
            if ~isempty(self.data)
                % Concatenate lightfield and confidence measurement
                lightfieldData      = cat( lightfield.enums.Formats.INTERNAL_FORMAT().channelDimension ...
                                         , self.data, self.confidence );
                temporaryLightfield = lightfield.Lightfield(lightfieldData);

                % Shear lightfield
                temporaryLightfield = temporaryLightfield.shear( dk_di ...
                                                               , shearMethod ...
                                                               , referenceViewpoint ...
                                                               , memoryEfficient );

                % Update lightfield
                self.data = temporaryLightfield.data(1:self.size.numberChannels,:,:,:,:);
                
                % Update confidence data if there is data in the initial
                % instance.
                if ~isempty(self.confidence)
                    self.confidence = temporaryLightfield.data(end,:,:,:,:);
                end
            end
        end
                
        function renderedImage = render(self,aperture)
            %
            % Obtain image from rendering the entire lightfield.
            %
            % INPUTS: 
            %   1. aperture - define viewpoints to be added to the rendered
            %   image.
            %
            narginchk(1,2);
            
            if nargin <= 1
                aperture = [1,self.size.numberPixels_i;1,self.size.numberPixels_j;nan,nan;nan,nan];
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

                % Weighted lightfield data
                confidenceWeights = repmat(self.confidence,3,1,1,1,1);
                lightfieldData    = self.data .* confidenceWeights;

                % Integrate aperture images to render image, setting to
                % zero the ligthfield coordinates with low precision
                % Additionally, consider that border pixels are not defined
                % for every viewpoint image. Thus, the border pixels must
                % be normalized by dividing by the weights of valid rays 
                % captured by the lightfield.
                weights        = nansum( nansum( confidenceWeights(:,minimum_i:maximum_i,minimum_j:maximum_j,:,:), 2), 3);
                invalidWeights = abs(weights) <= 10 * eps;
                
                imageData = nansum( nansum( lightfieldData(:,minimum_i:maximum_i,minimum_j:maximum_j,:,:), 2), 3);
                imageData = imageData ./ weights;
                imageData(invalidWeights) = 0;
                imageData = permute(imageData,lightfield.image.ViewpointImage.ENCODING_ORDER);

                % Obtain image instance
                renderedImage = image.Image(imageData);
            end
        end
    end
    
    methods (Static)
        function self = LightfieldFromLytroImage( cameraInstance, imageFilepath ...
                                                , decodingMethod ...
                                                , correctHexagonalSampling ...
                                                , squareMicrolenses ...
                                                , hotPixelCorrection ...
                                                , weightedDemosaicing ...
                                                , weightedInterpolation )
            % 
            % Create a decoded lightfiled instance from a lytro image using
            % the camera white images database path.
            %
            % Decode lightfield from raw plenoptic image in 2D. This method
            % allows to read all lytro files (.lfr, .lfp, *.raw) and it
            % will store the raw lenslet and white image. The white image
            % is used to correct vignetting.
            % 
            % INPUTS:
            %   1. imageFilepath - lytro image filepath
            %   2. correctHexagonalSampling - flag to indicate if hexagonal
            %   sampling must be corrected. Default is true.
            %   3. decodingMethod - decoding method. Default is Dansereau's
            %   decoding.
            %   4. squareMicrolenses - flag to indicate if microlenses
            %   resolution must be equal on the 2 dimensions. Default is
            %   true.
            %   5. hotPixelCoorection - flag to indicate if hot pixel
            %   correction using raw black images should be performed.
            %   Default is true.
            %   6. weightedDemosaicing - flag to indicate if white image
            %   guided demosaicing should be performed. Default is true.
            %   7. weightedInterpolation - flag to indicate if white image
            %   guided interpolation should be performed. Default is true.
            %
            narginchk(2,8);
            
            if nargin <= 2
                decodingMethod = lightfield.enums.DecodingMethods.DANSEREAU();
            end
            
            if nargin <= 3
                correctHexagonalSampling = true;
            end
            
            if nargin <= 4
                squareMicrolenses = true;
            end
            
            if nargin <= 5
                hotPixelCorrection = true;
            end

            if nargin <= 6
                weightedDemosaicing = true;
            end

            if nargin <= 7
                weightedInterpolation = true;
            end

            % Create instance
            self = lightfield.lytro.Lightfield();
            
            % Obtain white images database
            if ~exist(cameraInstance.whiteImagesFilepath,'file')
                % Assume that lytro archive folder path is in the root path
                % of white images and generate white images database
                [lytroArchivePath,~,~] = fileparts(cameraInstance.whiteImagesFilepath);
                camera.lytro.Camera.obtainWhiteImagesFromLytroArchive(lytroArchivePath);
            end

            % Define path for white images database
            decodeOptionsTemp = struct( 'WhiteImageDatabasePath', cameraInstance.whiteImagesFilepath ...
                                      , 'DecodingOutput', lightfield.enums.DecodingOutputs.LIGHTFIELD().symbol ...
                                      , 'ResampMethod', 'triangulation' ...
                                      , 'HotPixelCorrect', hotPixelCorrection ...
                                      , 'WeightedDemosaic', weightedDemosaicing ...
                                      , 'WeightedInterp', weightedInterpolation ...
                                      , 'DoDehex' , true ...
                                      , 'SquareST', squareMicrolenses );
                              
            if correctHexagonalSampling == false
                decodeOptionsTemp.ResampMethod = 'none';
            end
            
            % Dansereau's decoding method
            if decodingMethod == lightfield.enums.DecodingMethods.DANSEREAU()
                % Decode lightfield and obtain raw images from lenslet and
                % white images
                [ rawData, self.lightfieldMetadata, ~, self.lensletGridModel, self.decodeOptions ...
                , rawWhiteImage, rawLensletImage, demosaickedLensletImage ...
                , lensletImage , whiteImage ] = LFLytroDecodeImage(imageFilepath,decodeOptionsTemp);
            
            % KAIST decoding method
            else
                % Decode lightfield and obtain raw images from lenslet and
                % white images
                [ rawData, self.lightfieldMetadata, ~, self.lensletGridModel, self.decodeOptions ...
                , rawWhiteImage, rawLensletImage, demosaickedLensletImage ...
                , lensletImage , whiteImage ] = lightfield.lytro.utils.LFLytroDecodeImage_KAIST(imageFilepath,decodeOptionsTemp);
            end
            
            % Decode toolbox format to internal format
            rawData = self.TOOLBOX_FORMAT.decode(rawData);
            self.rawWhiteImage   = image.Image.MATLAB_FORMAT.decode(im2double(rawWhiteImage));
            self.rawLensletImage = image.Image.MATLAB_FORMAT.decode(im2double(rawLensletImage));
            self.lensletImage    = image.Image.MATLAB_FORMAT.decode(im2double(lensletImage));
            self.whiteImage      = image.Image.MATLAB_FORMAT.decode(im2double(whiteImage));
            self.demosaickedLensletImage = image.Image.MATLAB_FORMAT.decode(im2double(demosaickedLensletImage));
            
            % Obtain microlens centers
            self.microlensCenters = lightfield.lytro.Lightfield.obtainMicrolensesCentersFromLensletGridModel(self.lensletGridModel);
            
            % The lightfield data corresponds to the first 3 channels of
            % the raw data obtained. The data is converted to double during
            % the assignment.
            self.data = LFConvertToFloat( rawData(1:3,:,:,:,:), 'double' );
            
            % The confidence measurements correspond to the 4th channel of
            % the raw data obtained.
            self.confidence = LFConvertToFloat( rawData(4,:,:,:,:), 'double' );
        end
        
        function self = LightfieldFromDansereauFile(filepath)
            % 
            % Create a lightfiled instance from a decoded lytro image saved
            % on a Dansereau's file structure.
            %
            % INPUTS:
            %   1. filepath - decoded lytro image filepath (Dansereau's
            %   structure).
            %
            narginchk(1,1);
            
            % Create instance
            self = lightfield.lytro.Lightfield();

            % Load Dansereau's decode lightfield
            fileData = load(filepath);
            
            % Obtain information for lightfield
            self.decodeOptions      = fileData.DecodeOptions;
            self.lightfieldMetadata = fileData.LFMetadata;
            self.lensletGridModel   = fileData.LensletGridModel;
            
            % Obtain microlens centers
            self.microlensCenters = lightfield.lytro.Lightfield.obtainMicrolensesCentersFromLensletGridModel(self.lensletGridModel);

            % Obtain lightfield and decode toolbox format to internal
            % format
            rawData = fileData.LF;
            rawData = self.TOOLBOX_FORMAT.decode(rawData);
            
            % The lightfield data corresponds to the first 3 channels of
            % the raw data obtained. The data is converted to double during
            % the assignment.
            self.data = LFConvertToFloat( rawData(1:3,:,:,:,:), 'double' );
            
            % The confidence measurements correspond to the 4th channel of
            % the raw data obtained.
            self.confidence = LFConvertToFloat( rawData(4,:,:,:,:), 'double' );
        end
        
        function microlensesCenters = obtainMicrolensesCentersFromLensletGridModel(lensletGridModel)
            %
            % Obtain microlenses centers from lenslet grid model.
            %
            % INPUTS:
            %   1. lensletGridModel - lenslet grid model.
            %
            narginchk(1,1);
            
            % Obtain microlenses centers
            microlensesCenters   = lightfield.lytro.utils.LFBuildHexGrid( lensletGridModel );
            microlensesCenters_u = microlensesCenters(:,:,1);
            microlensesCenters_v = microlensesCenters(:,:,2);
            microlensesCenters   = image.Pixel( [microlensesCenters_u(:)';microlensesCenters_v(:)'],false );
        end
        
        function lenletGridModel = obtainLensletGridModelFromWhiteImage(whiteImagePath,radius)
            %
            % Obtain lenslet grid model from white image.
            %
            % INPUTS:
            %   1. whiteImagePath - filepath to white image file. It must
            %   be a file extension that can be read using imread.
            %   2. radius - radius considered for sampling. According to
            %   Bok, one should use 5 for 1st generation Lytro and 7 for
            %   Lytro Illum.
            %
            narginchk(1,2);
            
            if nargin <= 1
                radius = 5;
            end
            
            NO_DEBUG = false;
            
            % Define default grid model options
            gridModelOptions = struct( 'Orientation', 'horz' ...
                                     , 'FilterDiskRadiusMult', 1/3 ...
                                     , 'CropAmt', 25 ...
                                     , 'SkipStep', 250 ...
                                     , 'ApproxLensletSpacing', 2 * radius );

            % Obtain white image data and convert to gray level image
            whiteImage = image.Image.ImageFromFile(whiteImagePath);
            whiteImage = whiteImage.rgb2gray();
            
            % Transform to Matlab format
            whiteImage_matlab = whiteImage.MATLAB_FORMAT().encode(whiteImage.data);
            
            % Obtain lenslet grid model
            lenletGridModel = LFBuildLensletGridModel(whiteImage_matlab, gridModelOptions, NO_DEBUG);
        end
        
        function microlensesCenters = obtainMicrolensesCenterFromWhiteImage(whiteImagePath,radius)
            %
            % Obtain microlenses centers from white image. This code
            % obtains an irregular position for the microlenses centers.
            % For a regular position of the microlenses centers, please use
            % obtainLensletGridModelFromWhiteImages.
            %
            % INPUTS:
            %   1. whiteImagePath - filepath to white image file. It must
            %   be a file extension that can be read using imread.
            %   2. radius - radius considered for sampling. According to
            %   Bok, one should use 5 for 1st generation Lytro and 7 for
            %   Lytro Illum.
            %
            narginchk(1,2);
            
            if nargin <= 1
                radius = 5;
            end
            
            % Obtain microlenses centers using Bok's method
            microlensesCenters = image.Pixel( lightfield.lytro.utils.obtainMicrolensesCenters(radius,whiteImagePath) ...
                                            , false );
        end
    end
end
