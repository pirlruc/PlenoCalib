classdef ViewpointProjectionMatrix < abstract.TemplatePlenopticProjectionMatrix
    %VIEWPOINTPROJECTIONMATRIX
    %   Viewpoint projection matrix camera instance
    
    properties (Dependent)
        incrementIntrinsicMatrix_i      % Increment intrinsic matrix data for the i-coordinate
        incrementIntrinsicMatrix_j      % Increment intrinsic matrix data for the j-coordinate
        incrementExtrinsicMatrix_i  	% Increment extrinsic matrix data for the i-coordinate
        incrementExtrinsicMatrix_j  	% Increment extrinsic matrix data for the j-coordinate
    end
    
    methods
        function self = ViewpointProjectionMatrix(varargin)
            %
            % Create viewpoint projection matrix instance.
            %
            % INPUTS:
            %   1. intrinsicMatrix - intrinsic matrix data.
            %   2. extrinsicMatrix - extrinsic matrix data.
            %   3. incrementIntrinsicMatrix_i - increment intrinsic matrix 
            %   data for i-coordinate.
            %   4. incrementIntrinsicMatrix_j - increment intrinsic matrix 
            %   data for j-coordinate.
            %   5. incrementExtrinsicMatrix_i - increment extrinsic matrix 
            %   data for i-coordinate.
            %   6. incrementExtrinsicMatrix_j - increment extrinsic matrix 
            %   data for j-coordinate.
            %
            narginchk(0,6);
            
            self = self@abstract.TemplatePlenopticProjectionMatrix;
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.intrinsicMatrix = varargin{1};
                end
                
                if nargin >= 2
                    self.extrinsicMatrix = varargin{2};
                end
                
                if nargin >= 3
                    self.incrementIntrinsicMatrix_x = varargin{3};
                end
                
                if nargin >= 4
                    self.incrementIntrinsicMatrix_y = varargin{4};
                end
                
                if nargin >= 5
                    self.incrementExtrinsicMatrix_x = varargin{5};
                end
                
                if nargin >= 6
                    self.incrementExtrinsicMatrix_y = varargin{6};
                end
            end
        end
        
        function increment = get.incrementIntrinsicMatrix_i(self)
            increment = self.incrementIntrinsicMatrix_x;
        end
        
        function increment = get.incrementIntrinsicMatrix_j(self)
            increment = self.incrementIntrinsicMatrix_y;
        end
        
        function increment = get.incrementExtrinsicMatrix_i(self)
            increment = self.incrementExtrinsicMatrix_x;
        end
        
        function increment = get.incrementExtrinsicMatrix_j(self)
            increment = self.incrementExtrinsicMatrix_y;
        end
    end
    
    methods
        function rays_ijkl = project(self,worldPoints,pixels_ij)
            %
            % Obtain projection of world points in image sensor.
            %
            % INPUTS: 
            %   1. worldPoints - points defined in the world coordinate
            %   system. The points should not be given in homogeneous 
            %   coordinates.
            %   2. pixels_ij - viewpoint camera indices to obtain the
            %   corresponding projections.
            % 
            narginchk(3,3);
            
            rays_ijkl = self.project@abstract.TemplatePlenopticProjectionMatrix( worldPoints ...
                                                                               , pixels_ij );
        end
        
        function rays_ijkl = reconstruct(self,rays_ijkl)
            %
            % Reconstruct world point from pixels in image sensor.
            %
            % INPUTS: 
            %   1. rays_ijkl - image ray coordinates for a given point.
            % 
            narginchk(2,2);
            
            CAMERA_TYPE = camera.models.plenoptic.enums.CameraTypes.VIEWPOINT(); 
            rays_ijkl   = self.reconstruct@abstract.TemplatePlenopticProjectionMatrix( rays_ijkl ...
                                                                                     , CAMERA_TYPE );
        end
        
        function projectionMatrices = obtainProjectionMatrices(self,pixels_ij)
            %
            % Obtain projection matrices for viewpoint cameras.
            %
            % INPUTS:
            %   1. pixels_ij - pixel indices to obtain projection matrices.
            %      Each pixel should be provided in different columns.
            %
            narginchk(1,2);
            
            % If pixels are not provided, obtain a regular sampling of
            % pixels
            if nargin <= 1
                [pixels_i,pixels_j] = meshgrid(1:4:10,1:4:10);
                pixels_ij = [pixels_i(:),pixels_j(:)]';
            end

            % Obtain projection matrices
            projectionMatrices = self.obtainProjectionMatrices@abstract.TemplatePlenopticProjectionMatrix(pixels_ij);
        end
    end
    
    methods (Static)
        function self = decode(parameters)
            %
            % Decode intrinsic and extrinsic parameters associated with
            % camera array.
            %
            % INPUTS:
            %   1. parameters - camera array intrinsic and extrinsic
            %   parameters.
            %
            narginchk(1,1);
            
            % Obtain template projection matrix
            template = abstract.TemplatePlenopticProjectionMatrix.decode(parameters);
            
            % Transform to viewpoint projection matrix
            self = camera.models.plenoptic.standard.ViewpointProjectionMatrix();
            self.intrinsicMatrix            = template.intrinsicMatrix;
            self.extrinsicMatrix            = template.extrinsicMatrix;
            self.incrementIntrinsicMatrix_x = template.incrementIntrinsicMatrix_x;
            self.incrementIntrinsicMatrix_y = template.incrementIntrinsicMatrix_y;
            self.incrementExtrinsicMatrix_x = template.incrementExtrinsicMatrix_x;
            self.incrementExtrinsicMatrix_y = template.incrementExtrinsicMatrix_y;
        end
        
        function projectionCenters = obtainProjectionCentersFromStandardPlenopticCamera(plenopticCamera,pixels_ij,warnings)
            %
            % Obtain projection centers for viewpoint cameras. In this
            % situation we are considering that (i,j) is fixed and (k,l)
            % can have different values.
            %
            % INPUTS:
            %   1. pixels_ij - pixel indices to obtain projection centers.
            %   Each pixel should be provided in different columns.
            %   2. warnings - flag to indicate if warnings should be
            %   displayed. Default is false.
            %
            narginchk(1,3);
            
            % If pixels are not provided, obtain a regular sampling of
            % pixels
            if nargin <= 1
                [pixels_i,pixels_j] = meshgrid(1:4:10,1:4:10);
                pixels_ij = [pixels_i(:),pixels_j(:)]';
            end
            pixels_ij = utils.enums.Classes.PIXEL().convert(pixels_ij);
            pixels_i  = pixels_ij.u;
            pixels_j  = pixels_ij.v;
            
            if nargin <= 2
                warnings = false;
            end
                    % Convert to plenoptic camera instance
            plenopticCamera = utils.enums.Classes.STANDARD_PLENOPTIC_CAMERA().convert(plenopticCamera);

            % Let us obtain the depth value from the derivatives ds/du and
            % dt/dv considering that they depend on k and l.
            z_su = -plenopticCamera.ds_dk * plenopticCamera.dk_du;
            z_tv = -plenopticCamera.dt_dl * plenopticCamera.dl_dv;
            
            % If the depth values for the projection center are not equal,
            % give warning and display the difference in depth values.
            if abs(z_su - z_tv) > eps && warnings == true
                warning( 'StandardPlenopticCamera:obtainViewpointsProjectionCenters' ...
                       , [ 'Depth values for projection center are different.\n' ...
                         , 'Depth (s,u) = %f\n', 'Depth (t,v) = %f\n', 'Difference  = %f' ] ...
                       , z_su,z_tv,abs(z_su - z_tv) );
            end
            
            % Consider that the depth of the center of projection is one of
            % the values obtained. We are assuming symmetry between (i,k)
            % and (j,l).
            % Obtain the remaining coordinates for the projection centers
            projectionCenters = math.Point( [ plenopticCamera.h_si .* pixels_i ...
                                            + plenopticCamera.h_s ...
                                            + z_su .* ( plenopticCamera.h_ui .* pixels_i ...
                                                      + plenopticCamera.h_u )...
                                            ; plenopticCamera.h_tj .* pixels_j ...
                                            + plenopticCamera.h_t ...
                                            + z_su .* ( plenopticCamera.h_vj .* pixels_j ...
                                                      + plenopticCamera.h_v ) ...
                                            ; z_su .* ones(1,pixels_ij.numberVectors) ], false);
        end
        
        function [baselineLengths, baselines] = obtainBaselinesFromStandardPlenopticCamera(plenopticCamera,pixels_ij)
            % 
            % Obtain baselines between different viewpoints centers of
            % projection.
            %
            % INPUTS:
            %   1. pixels_ij - pixel indices to obtain projection centers.
            %      Each pixel should be provided in different columns.
            %
            narginchk(1,2);
            
            % If pixels are not provided, obtain a regular sampling of
            % pixels
            if nargin <= 1
                [pixels_i,pixels_j] = meshgrid(1:4:10,1:4:10);
                pixels_ij = [pixels_i(:),pixels_j(:)]';
            end
            pixels_ij = utils.enums.Classes.PIXEL().convert(pixels_ij);
            delta_i   = pixels_ij.u;
            delta_j   = pixels_ij.v;
            
            % Convert to plenoptic camera instance
            plenopticCamera = utils.enums.Classes.STANDARD_PLENOPTIC_CAMERA().convert(plenopticCamera);

            % Obtain baselines in the x- and y-axis
            baselines_x = delta_i * (plenopticCamera.h_si - plenopticCamera.h_ui * plenopticCamera.h_sk / plenopticCamera.h_uk);
            baselines_y = delta_j * (plenopticCamera.h_tj - plenopticCamera.h_vj * plenopticCamera.h_sk / plenopticCamera.h_uk);
            baselines   = [baselines_x; baselines_y];
            
            % Since there is an x- and y-component for the baseline, the
            % baseline is given by the norm
            baselineLengths = sqrt(baselines_x.^2 + baselines_y.^2);
        end
        
        function self = ViewpointProjectionMatrixFromStandardPlenopticCamera(plenopticCamera)
            %
            % Obtain parameters to describe a camera array from a standard
            % plenoptic camera. In this situation we are considering that 
            % (i,j) is fixed and (k,l) can have different values.
            %
            narginchk(1,1);
            
            % Convert to plenoptic camera instance
            plenopticCamera = utils.enums.Classes.STANDARD_PLENOPTIC_CAMERA().convert(plenopticCamera);
            
            % Obtain projection matrix for viewpoint i = j = 0 and for i =
            % j = 1.
            pixels_ij = image.Pixel([0,1;0,1],false);
            
            % Obtain the projection centers of the viewpoint cameras
            projectionCenters = camera.models.plenoptic.standard.ViewpointProjectionMatrix.obtainProjectionCentersFromStandardPlenopticCamera( ...
                                        plenopticCamera, pixels_ij );
            
            % Define intrinsic parameters
            % Remember that in pixel notation u = i and v = j
            pixelSize      = [ 1 / abs(plenopticCamera.h_uk), 1 / abs(plenopticCamera.h_vl) ];
            focalLength    = 1;
            skewFactor     = -plenopticCamera.h_ul / (plenopticCamera.h_uk * plenopticCamera.h_vl);
            principalPoint = image.Pixel ( - [ ( plenopticCamera.h_ui  .* pixels_ij.u ...
                                               + plenopticCamera.h_u ) ./ plenopticCamera.h_uk ...
                                             + (plenopticCamera.h_ul / plenopticCamera.h_uk) ...
                                             * ( plenopticCamera.h_vj  .* pixels_ij.v ...
                                               + plenopticCamera.h_v ) ./ plenopticCamera.h_vl ...
                                             ; ( plenopticCamera.h_vj  .* pixels_ij.v ...
                                               + plenopticCamera.h_v ) ./ plenopticCamera.h_vl ] );
                    
            % Define the extrinsic parameters. The rotation matrix
            % corresponds to the rotation matrix transforming the world to
            % the camera coordinate system. The translation corresponds to
            % the difference between the translation of the world to camera
            % coordinate system and the projection center of the viewpoint
            rotation    = plenopticCamera.extrinsic.rotationMatrix;
            translation = plenopticCamera.extrinsic.translationVector;

            % Correct extrinsics and projection centers if needed
            rotation.data(1,:)          = sign(plenopticCamera.h_uk) .* rotation.data(1,:);
            rotation.data(2,:)          = sign(plenopticCamera.h_vl) .* rotation.data(2,:);
            translation.data(1,:)       = sign(plenopticCamera.h_uk) .* translation.data(1,:);
            translation.data(2,:)       = sign(plenopticCamera.h_vl) .* translation.data(2,:);
            projectionCenters.data(1,:) = sign(plenopticCamera.h_uk) .* projectionCenters.data(1,:);
            projectionCenters.data(2,:) = sign(plenopticCamera.h_vl) .* projectionCenters.data(2,:);
            
            % Projection matrix for i = j = 0
            % Obtain intrinsic matrix
            intrinsic_ij_0 = camera.models.pinhole.IntrinsicMatrix( pixelSize ...
                                                                  , principalPoint.obtainVectors(1) ...
                                                                  , focalLength ...
                                                                  , skewFactor );
                
            % Obtain the extrinsic matrix
            extrinsic_ij_0 = camera.models.ExtrinsicMatrix( rotation ...
                                                          , translation.data - projectionCenters.obtainVectors(1).data );
                                                      
            % Projection matrix for i = j = 1
            % Obtain intrinsic matrix
            intrinsic_ij_1 = camera.models.pinhole.IntrinsicMatrix( pixelSize ...
                                                                  , principalPoint.obtainVectors(2) ...
                                                                  , focalLength ...
                                                                  , skewFactor );
                
            % Obtain the extrinsic matrix
            extrinsic_ij_1 = camera.models.ExtrinsicMatrix( rotation ...
                                                          , translation.data - projectionCenters.obtainVectors(2).data );
                
            % Obtain increment intrinsic and extrinsic matrix
            incrementIntrinsic = intrinsic_ij_1.intrinsicMatrix - intrinsic_ij_0.intrinsicMatrix;
            incrementExtrinsic = extrinsic_ij_1.extrinsicMatrix - extrinsic_ij_0.extrinsicMatrix;
            
            % Let us isolate each of the components
            incrementIntrinsic_i = zeros(3,3);
            incrementIntrinsic_i(1,:) = incrementIntrinsic(1,:);
            incrementIntrinsic_j = zeros(3,3);
            incrementIntrinsic_j(2,:) = incrementIntrinsic(2,:);
            incrementExtrinsic_i = zeros(3,4);
            incrementExtrinsic_i(1,:) = incrementExtrinsic(1,:);
            incrementExtrinsic_j = zeros(3,4);
            incrementExtrinsic_j(2,:) = incrementExtrinsic(2,:);
            
            % Create viewpoint projection matrix
            self = camera.models.plenoptic.standard.ViewpointProjectionMatrix();
            self.intrinsicMatrix  = intrinsic_ij_0.intrinsicMatrix;
            self.extrinsicMatrix  = extrinsic_ij_0.extrinsicMatrix;
            self.incrementIntrinsicMatrix_x = incrementIntrinsic_i;
            self.incrementIntrinsicMatrix_y = incrementIntrinsic_j;
            self.incrementExtrinsicMatrix_x = incrementExtrinsic_i;
            self.incrementExtrinsicMatrix_y = incrementExtrinsic_j;
        end
        
        function arrays = ViewpointProjectionMatrixFromHomographies(homographies)
            %
            % Estimate intrinsic and extrinsic parameters that characterize 
            % the camera array from homographies. This method obtains a 
            % collection of viewpoint projection matrices corresponding to
            % the different homographies.
            %
            % INPUTS:
            %   1. homographies - collection of camera array homographies.
            %
            narginchk(1,1);

            % This should be the number of viewpoints, but we can use just
            % these camera indices
            cameras_ij = [ 0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9,0 ,10,0 ,11,0 ,12,0 ,13,0 ,14,0 ,15,0 ,20 ...
                         ; 1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9,0,10,0 ,11,0 ,12,0 ,13,0 ,14,0 ,15,0 ,20,0 ];

            % Perform rectangular sampling estimation
			intrinsics  = abstract.TemplatePlenopticProjectionMatrix.IntrinsicsFromHomographies( ...
								homographies, cameras_ij );
			projections = abstract.TemplatePlenopticProjectionMatrix.ExtrinsicMatrixFromHomographies( ...
								intrinsics, homographies );
            
            
            % Correct class for estimated projection matrices
            arrays = [];
            for iViewpoint = 1:length(projections)
                % Obtain projection matrix for pose
                poseProjection = projections(iViewpoint);
                
                % Define camera array
                cameraArray = camera.models.plenoptic.standard.ViewpointProjectionMatrix();
                cameraArray.intrinsicMatrix = poseProjection.intrinsicMatrix;
                cameraArray.extrinsicMatrix = poseProjection.extrinsicMatrix;
                cameraArray.incrementIntrinsicMatrix_x = poseProjection.incrementIntrinsicMatrix_x;
                cameraArray.incrementIntrinsicMatrix_y = poseProjection.incrementIntrinsicMatrix_y;
                cameraArray.incrementExtrinsicMatrix_x = poseProjection.incrementExtrinsicMatrix_x;
                cameraArray.incrementExtrinsicMatrix_y = poseProjection.incrementExtrinsicMatrix_y;
                
                % Update collection of arrays
                arrays = cat(2,arrays,cameraArray);
            end
        end
    end
end
