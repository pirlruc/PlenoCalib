classdef SurfaceCameraCollection
    %SURFACECAMERACOLLECTION
    %   Surface camera collection methods and properties.

    properties
        data      = {}                          % Collection of surface camera images
        reference = lightfield.ImageRay()       % Lightfield viewpoint reference
    end
    
    properties (Dependent)
        size            % Surface camera collection size
        vectorize       % Vectorize surface camera collection
    end

    methods
        function self = SurfaceCameraCollection(varargin)
            %
            % Surface camera collection instance
            %
            % INPUTS:
            %   1. surfaceCameraImageData - surface camera image data.
            %   2. referenceViewpoint - reference viewpoint considered to
            %   obtain surface camera collection.
            %
            narginchk(0,2);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.data = varargin{1};
                end
                
                if nargin >= 2
                    self.reference = varargin{2};
                end
            end
        end
        
        function self = set.reference(self,newReference)
            self.reference = utils.enums.Classes.IMAGE_RAY().convert(newReference);
        end
        
        function imageSize = get.size(self)
            % If data is empty, consider the default size values.
            if isempty(self.data)
                imageSize = lightfield.image.ViewpointImageSize();
            else
                imageSize = lightfield.image.ViewpointImageSize([1,size(self.data)]);
            end
        end
        
        function vector = get.vectorize(self)
            % Vectorize surface camera collection. Surface camera 
            % collection has the size of an image with:
            %           channels x viewpoints x microlenses
            if isempty(self.data)
                vector = image.Image();
            else
                % Vectorize each surface camera image into image format
                vector = [self.data{:}];
                vector = [vector.vectorize];
                vector = [vector.data];
                
                % Obtain size to reshape final matrix
                numberChannels    = self.data{1,1}.size.numberChannels;
                numberViewpoints  = self.data{1,1}.size.numberPixels_u * self.data{1,1}.size.numberPixels_v;
                numberMicrolenses = self.size.numberMicrolenses_k * self.size.numberMicrolenses_l;

                % Reshape final matrix
                vector = image.Image( reshape(vector, [numberChannels,numberViewpoints,numberMicrolenses]) );
            end
        end
    end
    
    methods 
        function renderedImage = render(self)
            %
            % Obtain image from rendering the surface camera collection.
            %
            
            if isempty(self.data)
                renderedImage = image.Image();
            else
                % Render each surface camera image
                renderedImage = cellfun(@render,self.data,'UniformOutput',false);
                
                % Obtain size to reshape final matrix
                numberChannels = renderedImage{1,1}.size.numberChannels;
                renderedImageSize = self.size.data;
                renderedImageSize(image.enums.Formats.INTERNAL_FORMAT().channelDimension) = numberChannels;
                
                % Convert data to image format
                renderedImage = [renderedImage{:}];
                renderedImage = [renderedImage.data];
                renderedImage = image.Image( reshape(renderedImage,renderedImageSize) );
            end
        end
    end
    
    methods (Static)
        function self = SurfaceCameraCollectionFromLightfield( lightfieldInstance ...
                                                             , disparity_dk_di ...
                                                             , interpolationMethod ...
                                                             , referenceViewpoint ...
                                                             , switchSamplingMethods ...
                                                             , memoryEfficient )
            %
            % Create surface camera image instance from lightfield data.
            %
            % INPUTS:
            %   1. lightfieldInstance - lightfield instance.
            %   2. disparity_dk_di    - disparity on viewpoint images to 
            %   create surface camera image. This encodes the depth of 
            %   the center of projection of the surface camera.
            %   3. interpolationMethod - interpolation method to sample
            %   lightfield. Default is linear.
            %   4. referenceViewpoint - reference viewpoint considered to
            %   create the collection of surface camera images. The default 
            %   corresponds to the central viewpoint.
            %   5. switchSamplingMethods - flag that indicates to switch
            %   automatically between sampling methods for surface camera
            %   image. Default is false.
            %   6. memoryEfficient - flag that indicates if the sampling on
            %   the viewpoints (i,j) should consider only the filled
            %   values. Default is true.
            %
            narginchk(1,6);
            
            EXTRAPOLATION_VALUE = nan;
            
            % Convert input data to instances
            lightfieldInstance = utils.enums.Classes.LIGHTFIELD().convert(lightfieldInstance);

            % Define default disparity on viewpoint images.
            if nargin <= 1
                disparity_dk_di = 1;
            end
            
            % Define default interpolation method
            if nargin <= 2
                interpolationMethod = lightfield.enums.ShearMethods.LINEAR;
            end
            
            % Define image ray coordinates position
            if nargin <= 3
                referenceViewpoint = [ round(lightfieldInstance.size.numberPixels_i / 2) ...
                                     ; round(lightfieldInstance.size.numberPixels_j / 2) ...
                                     ; nan; nan ];
            end
            
            if nargin <= 4
                switchSamplingMethods = false;
            end
            
            if nargin <= 5
                memoryEfficient = true;
            end
            
            % Convert image ray coordinates to image ray
            referenceViewpoint = utils.enums.Classes.IMAGE_RAY().convert(referenceViewpoint);
            
            % Define threshold to exchange between methods for creating
            % surface camera image
            dk_di_threshold = +inf;
            if switchSamplingMethods == true
                dk_di_threshold = 1;
            end
            
            % Initialize surface camera collection
            self = lightfield.scams.SurfaceCameraCollection();
            self.data = cell( lightfieldInstance.size.numberMicrolenses_k ...
                            , lightfieldInstance.size.numberMicrolenses_l );
            self.reference = referenceViewpoint;
            
            % For each microlens of the lightfield, create a surface camera
            % image
            if abs(disparity_dk_di) <= dk_di_threshold
                % Shear lightfield
                lightfieldInstance = lightfieldInstance.shear( disparity_dk_di, interpolationMethod ...
                                                             , referenceViewpoint, memoryEfficient );
                                                         
                % Reorganize in microlens images
                microlensImages = permute( lightfieldInstance.data ...
                                         , lightfield.image.MicrolensImage.ENCODING_ORDER );
                
                for microlens_k = 1:lightfieldInstance.size.numberMicrolenses_k
                    for microlens_l = 1:lightfieldInstance.size.numberMicrolenses_l
                        % Create surface camera image instance from 
                        % lightfield
                        scam = lightfield.image.SurfaceCameraImage();
                        scam.data      = microlensImages(:,:,:,microlens_k,microlens_l);
                        scam.reference = [ referenceViewpoint.i; referenceViewpoint.j ...
                                         ; microlens_k; microlens_l ];

                        % Update collection of surface cameras
                        self.data{microlens_k,microlens_l} = scam;
                    end
                end
                
            % Surface camera collections have many coordinates without
            % information. It is a good idea to reduce the computer memory
            % requirements when using these objects.
            elseif memoryEfficient == true
                % Obtain lightfield coordinates
                [pixels_i,pixels_j,microlenses_k,microlenses_l] = ...
                                ndgrid( 1:lightfieldInstance.size.numberPixels_i ...
                                      , 1:lightfieldInstance.size.numberPixels_j ...
                                      , 1:lightfieldInstance.size.numberMicrolenses_k ...
                                      , 1:lightfieldInstance.size.numberMicrolenses_l );

                % Obtain new lightfield coordinates sampling on (i,j)
                [newMicrolenses_k,newMicrolenses_l,coordinates_k,coordinates_l] = ...
                                ndgrid( 1:lightfieldInstance.size.numberMicrolenses_k ...
                                      , 1:lightfieldInstance.size.numberMicrolenses_l ...
                                      , 1:lightfieldInstance.size.numberMicrolenses_k ...
                                      , 1:lightfieldInstance.size.numberMicrolenses_l );
                newPixels_i = referenceViewpoint.i + (newMicrolenses_k - coordinates_k)./ disparity_dk_di;
                newPixels_j = referenceViewpoint.j + (newMicrolenses_l - coordinates_l)./ disparity_dk_di;

                % Obtain maximum size of the surface camera images
                indices_i   = newPixels_i >= 0 ...
                            & newPixels_i <= lightfieldInstance.size.numberPixels_i + 1;
                indices_j   = newPixels_j >= 0 ...
                            & newPixels_j <= lightfieldInstance.size.numberPixels_j + 1;
                indices_ij  = indices_i & indices_j;
                maximumNumberMicrolenses_k = max(max(max(sum(indices_ij,1))));
                maximumNumberMicrolenses_l = max(max(max(sum(indices_ij,2))));
                scamInitialization         = nan( lightfieldInstance.size.numberChannels ...
                                                , maximumNumberMicrolenses_k ...
                                                , maximumNumberMicrolenses_l );
                
				% Create surface camera image for each microlens
                for microlens_k = 1:lightfieldInstance.size.numberMicrolenses_k
                    for microlens_l = 1:lightfieldInstance.size.numberMicrolenses_l
                        % Obtain number of samplings on (i,j)
                        numberMicrolenses_k = max(sum(indices_ij(:,:,microlens_k,microlens_l),1));
                        numberMicrolenses_l = max(sum(indices_ij(:,:,microlens_k,microlens_l),2));

                        % Sub-select coordinates for interpolation and 
                        % obtain original microlens structure
                        microlensCoordinates_i = newPixels_i(:,:,microlens_k,microlens_l);
                        microlensCoordinates_j = newPixels_j(:,:,microlens_k,microlens_l);
                        microlensCoordinates_k = newMicrolenses_k(:,:,microlens_k,microlens_l);
                        microlensCoordinates_l = newMicrolenses_l(:,:,microlens_k,microlens_l);

                        microlensCoordinates_i = reshape( microlensCoordinates_i(indices_ij(:,:,microlens_k,microlens_l)) ...
                                                        , [numberMicrolenses_k,numberMicrolenses_l] );
                        microlensCoordinates_j = reshape( microlensCoordinates_j(indices_ij(:,:,microlens_k,microlens_l)) ...
                                                        , [numberMicrolenses_k,numberMicrolenses_l] );
                        microlensCoordinates_k = reshape( microlensCoordinates_k(indices_ij(:,:,microlens_k,microlens_l)) ...
                                                        , [numberMicrolenses_k,numberMicrolenses_l] );
                        microlensCoordinates_l = reshape( microlensCoordinates_l(indices_ij(:,:,microlens_k,microlens_l)) ...
                                                        , [numberMicrolenses_k,numberMicrolenses_l] );

                        % Initialize surface camera image
                        scam = lightfield.image.SurfaceCameraImage();
                        scam.reference = [ referenceViewpoint.i,referenceViewpoint.j ...
                                         , microlens_k, microlens_l ];
                        scam.data = scamInitialization;
                        
						% Interpolate data for surface camera
                        for iChannel = 1:lightfieldInstance.size.numberChannels
                            colorLightfield = permute(lightfieldInstance.data(iChannel,:,:,:,:),[2,3,4,5,1]);
                            colorLightfield = interpn( pixels_i, pixels_j, microlenses_k, microlenses_l, colorLightfield ...
                                                     , microlensCoordinates_i, microlensCoordinates_j, microlensCoordinates_k, microlensCoordinates_l ...
                                                     , interpolationMethod.symbol, EXTRAPOLATION_VALUE );
                            scam.data(iChannel,1:numberMicrolenses_k,1:numberMicrolenses_l) = colorLightfield;
                        end

                        % Update collection of surface cameras
                        self.data{microlens_k,microlens_l} = scam;
                    end
                end

            % Obtain surface camera collection without any considerations
            % with computer memory requirements
            else
                for microlens_k = 1:lightfieldInstance.size.numberMicrolenses_k
                    for microlens_l = 1:lightfieldInstance.size.numberMicrolenses_l
                        % Obtain surface camera image coordinates
                        imageRay_ijkl = [referenceViewpoint.i;referenceViewpoint.j;microlens_k;microlens_l];

                        % Create surface camera image instance from lightfield
                        scam = lightfield.image.SurfaceCameraImage.SurfaceCameraImageFromLightfield( lightfieldInstance ...
                                                                                                   , disparity_dk_di ...
                                                                                                   , interpolationMethod ...
                                                                                                   , imageRay_ijkl ...
                                                                                                   , switchSamplingMethods ...
                                                                                                   , memoryEfficient );

                        % Update collection of surface cameras
                        self.data{microlens_k,microlens_l} = scam;
                    end
                end
            end
        end
    end
end