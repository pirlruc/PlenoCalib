classdef SurfaceCameraImage < abstract.TemplateImage
    %SURFACECAMERAIMAGE
    %   Surface camera image methods and properties.

    properties
        reference = lightfield.ImageRay()       % Lightfield reference coordinate
    end
    
    properties (Dependent)
        size            % Surface camera image size
        vectorize       % Vectorize surface camera image
    end

    methods
        function self = SurfaceCameraImage(varargin)
            %
            % Surface camera image instance
            %
            % INPUTS:
            %   1. surfaceCameraImageData - surface camera image data.
            %   2. referenceImageRay - lightfield coordinate used to obtain
            %   surface camera image.
            %
            narginchk(0,2);
            
            % Create super class instance
            self = self@abstract.TemplateImage();
            
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
            % If image data is empty, consider the default image size
            % values.
            if isempty(self.data)
                imageSize = image.ImageSize();
            else
                imageSize = image.ImageSize(size(self.data));
            end
        end
        
        function vector = get.vectorize(self)
            % Vectorize surface camera image. Surface camera image has the 
            % size of an image row with:
            %           channels x viewpoints x microlenses (1)
            if isempty(self.data)
                vector = [];
            else
                vector = image.Image( ...
                                reshape( self.data, [ self.size.numberChannels ...
                                                    , self.size.numberPixels_u * self.size.numberPixels_v ] ));
            end
        end
    end
    
    methods 
        function pixelImage = render(self)
            %
            % Obtain pixel value from rendering the surface camera image.
            %
            
            if isempty(self.data)
                pixelImage = image.Image();
            else
                % Count number of defined pixels in one channel
                numberPixels = ~isnan(self.data(1,:,:));
                numberPixels = sum(numberPixels(:));
                
                % Integrate to get pixel value. Do not consider NaN values
                pixelImage = image.Image(nansum(nansum(self.data,2),3) ./ numberPixels);
            end
        end
    end
    
    methods (Static)
        function self = SurfaceCameraImageFromLightfield( lightfieldInstance ...
                                                        , disparity_dk_di ...
                                                        , interpolationMethod ...
                                                        , imageRay_ijkl ...
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
            %   4. imageRay_ijkl      - image ray coordinates. This encodes
            %   the (x,y) position of the center of projection of the
            %   surface camera.
            %   5. switchSamplingMethods - flag that indicates to switch
            %   automatically between sampling methods for surface camera
            %   image. Default is false.
            %   6. memoryEfficient - flag that indicates if the sampling on
            %   the viewpoints (i,j) should consider only the filled
            %   values. Default is false.
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
                imageRay_ijkl = [ round(lightfieldInstance.size.numberPixels_i / 2) ...
                                ; round(lightfieldInstance.size.numberPixels_j / 2) ...
                                ; round(lightfieldInstance.size.numberMicrolenses_k / 2) ...
                                ; round(lightfieldInstance.size.numberMicrolenses_l / 2) ];
            end
            
            if nargin <= 4
                switchSamplingMethods = false;
            end
            
            if nargin <= 5
                memoryEfficient = false;
            end
            
            % Convert image ray coordinates to image ray
            imageRay_ijkl = utils.enums.Classes.IMAGE_RAY().convert(imageRay_ijkl);

            % Define threshold to exchange between methods for creating
            % surface camera image
            dk_di_threshold = +inf;
            if switchSamplingMethods == true
                dk_di_threshold = 1;
            end
            
            % Obtain lightfield coordinates
            [pixels_i,pixels_j,microlenses_k,microlenses_l] = ...
                            ndgrid( 1:lightfieldInstance.size.numberPixels_i ...
                                  , 1:lightfieldInstance.size.numberPixels_j ...
                                  , 1:lightfieldInstance.size.numberMicrolenses_k ...
                                  , 1:lightfieldInstance.size.numberMicrolenses_l );

            % Obtain new lightfield coordinates
            % Sampling on (i,j)
            if abs(disparity_dk_di) <= dk_di_threshold
                [newPixels_i,newPixels_j] = ...
                                ndgrid( 1:lightfieldInstance.size.numberPixels_i ...
                                      , 1:lightfieldInstance.size.numberPixels_j );
                newMicrolenses_k = imageRay_ijkl.k + disparity_dk_di .* (newPixels_i - imageRay_ijkl.i);
                newMicrolenses_l = imageRay_ijkl.l + disparity_dk_di .* (newPixels_j - imageRay_ijkl.j);

                % Obtain surface camera image data from lightfield
                surfaceCameraImageData = zeros( lightfieldInstance.size.numberChannels ...
                                              , lightfieldInstance.size.numberPixels_i ...
                                              , lightfieldInstance.size.numberPixels_j );
            % Sampling on (k,l)
            else
                [newMicrolenses_k,newMicrolenses_l] = ...
                                ndgrid( 1:lightfieldInstance.size.numberMicrolenses_k ...
                                      , 1:lightfieldInstance.size.numberMicrolenses_l );
                newPixels_i = imageRay_ijkl.i + (newMicrolenses_k - imageRay_ijkl.k) ./ disparity_dk_di;
                newPixels_j = imageRay_ijkl.j + (newMicrolenses_l - imageRay_ijkl.l) ./ disparity_dk_di;

                % Most of the coordinates obtained will not have values.
                % Thus, select only values that are going to be filled. 
                % Consider one pixel neighborhood for ensuring that all 
                % values are selected.
                if memoryEfficient == true
                    indices_i  = newPixels_i >= 0 ...
                               & newPixels_i <= lightfieldInstance.size.numberPixels_i + 1;
                    indices_j  = newPixels_j >= 0 ...
                               & newPixels_j <= lightfieldInstance.size.numberPixels_j + 1;
                    indices_ij = indices_i & indices_j;
                    numberMicrolenses_k = max(sum(indices_ij,1));
                    numberMicrolenses_l = max(sum(indices_ij,2));

                    % Sub-select coordinates for interpolation
                    newPixels_i = reshape( newPixels_i(indices_ij) ...
                                         , [numberMicrolenses_k,numberMicrolenses_l] );
                    newPixels_j = reshape( newPixels_j(indices_ij) ...
                                         , [numberMicrolenses_k,numberMicrolenses_l] );
                    newMicrolenses_k = reshape( newMicrolenses_k(indices_ij) ...
                                              , [numberMicrolenses_k,numberMicrolenses_l] );
                    newMicrolenses_l = reshape( newMicrolenses_l(indices_ij) ...
                                              , [numberMicrolenses_k,numberMicrolenses_l] );

                else
                    numberMicrolenses_k = lightfieldInstance.size.numberMicrolenses_k;
                    numberMicrolenses_l = lightfieldInstance.size.numberMicrolenses_l;
                end
                
                % Obtain surface camera image data from lightfield
                surfaceCameraImageData = nan( lightfieldInstance.size.numberChannels ...
                                            , numberMicrolenses_k, numberMicrolenses_l );
            end
            
            for iChannel = 1:lightfieldInstance.size.numberChannels
                colorLightfield = permute(lightfieldInstance.data(iChannel,:,:,:,:),[2,3,4,5,1]);
                colorLightfield = interpn( pixels_i, pixels_j, microlenses_k, microlenses_l, colorLightfield ...
                                         , newPixels_i, newPixels_j, newMicrolenses_k, newMicrolenses_l ...
                                         , interpolationMethod.symbol, EXTRAPOLATION_VALUE );
                surfaceCameraImageData(iChannel,:,:) = colorLightfield;
            end

            % Create image instance
            self = lightfield.image.SurfaceCameraImage();
            self.data      = surfaceCameraImageData;
            self.reference = imageRay_ijkl;
        end
    end
end