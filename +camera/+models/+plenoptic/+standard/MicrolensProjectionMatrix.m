classdef MicrolensProjectionMatrix < abstract.TemplatePlenopticProjectionMatrix
    %MICROLENSPROJECTIONMATRIX
    %   Microlens projection matrix camera instance
    
    properties (Dependent)
        incrementIntrinsicMatrix_k  	% Increment intrinsic matrix data for k-coordinate
        incrementIntrinsicMatrix_l  	% Increment intrinsic matrix data for l-coordinate
        incrementExtrinsicMatrix_k  	% Increment extrinsic matrix data for k-coordinate
        incrementExtrinsicMatrix_l  	% Increment extrinsic matrix data for l-coordinate
    end
    
    methods
        function self = MicrolensProjectionMatrix(varargin)
            %
            % Create microlens projection matrix instance.
            %
            % INPUTS:
            %   1. intrinsicMatrix - intrinsic matrix data.
            %   2. extrinsicMatrix - extrinsic matrix data.
            %   3. incrementIntrinsicMatrix_k - increment intrinsic matrix 
            %   data for k-coordinate.
            %   4. incrementIntrinsicMatrix_l - increment intrinsic matrix 
            %   data for l-coordinate.
            %   5. incrementExtrinsicMatrix_k - increment extrinsic matrix 
            %   data for k-coordinate.
            %   6. incrementExtrinsicMatrix_l - increment extrinsic matrix 
            %   data for l-coordinate.
            %
            narginchk(0,6);
            
            self = self@abstract.TemplatePlenopticProjectionMatrix();
            
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
        
        function increment = get.incrementIntrinsicMatrix_k(self)
            increment = self.incrementIntrinsicMatrix_x;
        end
        
        function increment = get.incrementIntrinsicMatrix_l(self)
            increment = self.incrementIntrinsicMatrix_y;
        end
        
        function increment = get.incrementExtrinsicMatrix_k(self)
            increment = self.incrementExtrinsicMatrix_x;
        end
        
        function increment = get.incrementExtrinsicMatrix_l(self)
            increment = self.incrementExtrinsicMatrix_y;
        end
    end
    
    methods
        function rays_ijkl = project(self,worldPoints,microlenses_kl)
            %
            % Obtain projection of world points in image sensor.
            %
            % INPUTS: 
            %   1. worldPoints - points defined in the world coordinate
            %   system. The points should not be given in homogeneous 
            %   coordinates.
            %   2. microlenses_kl - microlenses indices to obtain the
            %   corresponding projections.
            % 
            narginchk(3,3);
            
            % Transform world points and camera indices to vector class 
            worldPoints   = utils.enums.Classes.POINT().convert(worldPoints);

            % The processing is made assuming viewpoint cameras.
            rays_ijkl = self.project@abstract.TemplatePlenopticProjectionMatrix( worldPoints ...
                                                                               , microlenses_kl );
                                                                     
            % Thus, for considering the microlenses cameras one should
            % switch the coordinates
            for iPoint = 1:worldPoints.numberVectors
                % Obtain point rays
                pointRays_ijkl = rays_ijkl(iPoint);
                rays_ijkl(iPoint).data = [ pointRays_ijkl.k ...
                                         ; pointRays_ijkl.l ...
                                         ; pointRays_ijkl.i ...
                                         ; pointRays_ijkl.j ];
            end
        end
        
        function rays_ijkl = reconstruct(self,rays_ijkl)
            %
            % Reconstruct world point from pixels in image sensor.
            %
            % INPUTS: 
            %   1. rays_ijkl - image ray coordinates for a given point.
            % 
            narginchk(2,2);
            
            CAMERA_TYPE = camera.models.plenoptic.enums.CameraTypes.MICROLENS(); 
            rays_ijkl   = self.reconstruct@abstract.TemplatePlenopticProjectionMatrix( rays_ijkl ...
                                                                                     , CAMERA_TYPE );
        end
        
        function projectionMatrices = obtainProjectionMatrices(self,microlenses_kl)
            %
            % Obtain projection matrices for microlens cameras.
            %
            % INPUTS:
            %   1. microlenses_kl - microlenses indices to obtain 
            %      projection matrices. Each microlens should be provided 
            %      in different columns.
            %
            narginchk(1,2);
            
            if nargin <= 1
                [microlenses_k,microlenses_l] = meshgrid(1:20:100,1:20:100);
                microlenses_kl = [microlenses_k(:),microlenses_l(:)]';
            end
            
            % Obtain projection matrices
            projectionMatrices = self.obtainProjectionMatrices@abstract.TemplatePlenopticProjectionMatrix(microlenses_kl);
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
            
            % Transform to microlens projection matrix
            self = camera.models.plenoptic.standard.MicrolensProjectionMatrix();
            self.intrinsicMatrix            = template.intrinsicMatrix;
            self.extrinsicMatrix            = template.extrinsicMatrix;
            self.incrementIntrinsicMatrix_x = template.incrementIntrinsicMatrix_x;
            self.incrementIntrinsicMatrix_y = template.incrementIntrinsicMatrix_y;
            self.incrementExtrinsicMatrix_x = template.incrementExtrinsicMatrix_x;
            self.incrementExtrinsicMatrix_y = template.incrementExtrinsicMatrix_y;
        end
        
        function projectionCenters = obtainProjectionCentersFromStandardPlenopticCamera(plenopticCamera,microlenses_kl,warnings)
            %
            % Obtain projection centers for microlens cameras. In this
            % situation we are considering that (k,l) is fixed and (i,j)
            % can have different values.
            %
            % INPUTS:
            %   1. microlenses_kl - microlenses indices to obtain
            %   projection centers. Each microlens should be provided in 
            %   different columns.
            %   2. warnings - flag to indicate if warnings should be
            %   displayed. Default is false.
            %
            narginchk(1,3);
            
            % If microlenses are not provided, obtain a regular sampling of
            % microlenses
            if nargin <= 1
                [microlenses_k,microlenses_l] = meshgrid(1:20:100,1:20:100);
                microlenses_kl = [microlenses_k(:),microlenses_l(:)]';
            end
            microlenses_kl = utils.enums.Classes.VECTOR().convert(microlenses_kl);
            
            if nargin <= 2
                warnings = false;
            end
            
            % If microlenses are provided in rows instead of columns, 
            % transpose the microlenses
            if microlenses_kl.numberComponents ~= 2
                microlenses_kl.data = microlenses_kl.data';
            end
            microlenses_k = microlenses_kl.obtainComponents(1).data;
            microlenses_l = microlenses_kl.obtainComponents(2).data;
            
            % Convert to plenoptic camera instance
            plenopticCamera = utils.enums.Classes.STANDARD_PLENOPTIC_CAMERA().convert(plenopticCamera);

            % Let us obtain the depth value from the derivatives ds/du and
            % dt/dv considering that they depend on i and j.
            z_su = -plenopticCamera.ds_di * plenopticCamera.di_du;
            z_tv = -plenopticCamera.dt_dj * plenopticCamera.dj_dv;
            
            % If the depth values for the projection center are not equal,
            % give warning and display the difference in depth values.
            if abs(z_su - z_tv) > eps && warnings == true
                warning( 'StandardPlenopticCamera:obtainMicrolensesProjectionCenters' ...
                       , [ 'Depth values for projection center are different.\n' ...
                         , 'Depth (s,u) = %f\n', 'Depth (t,v) = %f\n', 'Difference  = %f' ] ...
                       , z_su,z_tv,abs(z_su - z_tv) );
            end
            
            % Consider that the depth of the center of projection is one of
            % the values obtained. We are assuming symmetry between (i,k)
            % and (j,l).
            % Obtain the remaining coordinates for the projection centers
            projectionCenters = math.Point( [ plenopticCamera.h_sk .* microlenses_k ...
                                            + plenopticCamera.h_sl .* microlenses_l ...
                                            + plenopticCamera.h_s ...
                                            + z_su .* ( plenopticCamera.h_uk .* microlenses_k ...
                                                      + plenopticCamera.h_ul .* microlenses_l ...
                                                      + plenopticCamera.h_u )...
                                            ; plenopticCamera.h_tl .* microlenses_l ...
                                            + plenopticCamera.h_t ...
                                            + z_su .* ( plenopticCamera.h_vl .* microlenses_l ...
                                                      + plenopticCamera.h_v ) ...
                                            ; z_su .* ones(1,microlenses_kl.numberVectors) ], false);
        end
        
        function [baselineLengths, baselines] = obtainBaselinesFromStandardPlenopticCamera(plenopticCamera,microlenses_kl)
            % 
            % Obtain baselines between different microlenses centers of
            % projection.
            %
            % INPUTS:
            %   1. microlenses_kl - microlenses indices to obtain the
            %      baselines. Each microlens should be provided in 
            %      different columns.
            %
            narginchk(1,2);
            
            % If microlenses are not provided, obtain a regular sampling of
            % microlenses
            if nargin <= 1
                [microlenses_k,microlenses_l] = meshgrid(1:20:100,1:20:100);
                microlenses_kl = [microlenses_k(:),microlenses_l(:)]';
            end
            microlenses_kl = utils.enums.Classes.VECTOR().convert(microlenses_kl);
            
            % If microlenses are provided in rows instead of columns, 
            % transpose the microlenses
            if microlenses_kl.numberComponents ~= 2
                microlenses_kl.data = microlenses_kl.data';
            end
            delta_k = microlenses_kl.obtainComponents(1).data;
            delta_l = microlenses_kl.obtainComponents(2).data;
            
            % Convert to plenoptic camera instance
            plenopticCamera = utils.enums.Classes.STANDARD_PLENOPTIC_CAMERA().convert(plenopticCamera);

            % Obtain baselines in the x- and y-axis
            baselines_x = delta_k .* (plenopticCamera.h_sk - plenopticCamera.h_uk * plenopticCamera.h_si / plenopticCamera.h_ui) ...
                        + delta_l .* (plenopticCamera.h_sl - plenopticCamera.h_ul * plenopticCamera.h_si / plenopticCamera.h_ui);
            baselines_y = delta_l .* (plenopticCamera.h_tl - plenopticCamera.h_vl * plenopticCamera.h_si / plenopticCamera.h_ui);
            baselines   = [baselines_x; baselines_y];
            
            % Since there is an x- and y-component for the baseline, the
            % baseline is given by the norm
            baselineLengths = sqrt(baselines_x.^2 + baselines_y.^2);
        end
        
        function self = MicrolensProjectionMatrixFromStandardPlenopticCamera(plenopticCamera)
            %
            % Obtain parameters to describe a camera array from a standard
            % plenoptic camera. In this situation we are considering that 
            % (k,l) is fixed and (i,j) can have different values.
            %
            % INPUTS:
            %   1. plenopticCamera - plenoptic camera data.
            %
            narginchk(1,1);
            
            % Convert to plenoptic camera instance
            plenopticCamera = utils.enums.Classes.STANDARD_PLENOPTIC_CAMERA().convert(plenopticCamera);
            
            % Obtain projection matrix for microlens k = l = 0, for 
            % k = 1 and l = 0, and for k = 0 and l = 1.
            microlenses_kl = math.Vector([0,1,0;0,0,1],false);
            microlenses_k  = microlenses_kl.obtainComponents(1).data;
            microlenses_l  = microlenses_kl.obtainComponents(2).data;
            
            % Obtain the projection centers of the microlens cameras
            projectionCenters = camera.models.plenoptic.standard.MicrolensProjectionMatrix.obtainProjectionCentersFromStandardPlenopticCamera( ...
                                        plenopticCamera, microlenses_kl );
            
            % Define intrinsic parameters
            pixelSize      = [ 1 / abs(plenopticCamera.h_ui), 1 / abs(plenopticCamera.h_vj) ];
            focalLength    = 1;
            skewFactor     = 0;
            principalPoint = image.Pixel( - [ ( plenopticCamera.h_uk .* microlenses_k ...
                                              + plenopticCamera.h_ul .* microlenses_l ...
                                              + plenopticCamera.h_u) ./ plenopticCamera.h_ui ...
                                            ; ( plenopticCamera.h_vl .* microlenses_l ...
                                              + plenopticCamera.h_v) ./ plenopticCamera.h_vj ], false );
                    
            % Define the extrinsic parameters. The rotation matrix
            % corresponds to the rotation matrix transforming the world to
            % the camera coordinate system. The translation corresponds to
            % the difference between the translation of the world to camera
            % coordinate system and the projection center of the microlens
            rotation    = plenopticCamera.extrinsic.rotationMatrix;
            translation = plenopticCamera.extrinsic.translationVector;
                        
            % Correct extrinsics and projection centers if needed
            rotation.data(1,:)          = sign(plenopticCamera.h_ui) .* rotation.data(1,:);
            rotation.data(2,:)          = sign(plenopticCamera.h_vj) .* rotation.data(2,:);
            translation.data(1,:)       = sign(plenopticCamera.h_ui) .* translation.data(1,:);
            translation.data(2,:)       = sign(plenopticCamera.h_vj) .* translation.data(2,:);
            projectionCenters.data(1,:) = sign(plenopticCamera.h_ui) .* projectionCenters.data(1,:);
            projectionCenters.data(2,:) = sign(plenopticCamera.h_vj) .* projectionCenters.data(2,:);

            % Projection matrix for k = l = 0
            % Obtain intrinsic matrix
            intrinsic_0 = camera.models.pinhole.IntrinsicMatrix( pixelSize ...
                                                               , principalPoint.obtainVectors(1) ...
                                                               , focalLength ...
                                                               , skewFactor );
                
            % Obtain the extrinsic matrix
            extrinsic_0 = camera.models.ExtrinsicMatrix( rotation ...
                                                       , translation.data - projectionCenters.obtainVectors(1).data );
                
            % Projection matrix for k = 1 and l = 0
            % Obtain intrinsic matrix
            intrinsic_k = camera.models.pinhole.IntrinsicMatrix( pixelSize ...
                                                               , principalPoint.obtainVectors(2) ...
                                                               , focalLength ...
                                                               , skewFactor );
                
            % Obtain the extrinsic matrix
            extrinsic_k = camera.models.ExtrinsicMatrix( rotation ...
                                                       , translation.data - projectionCenters.obtainVectors(2).data );
                
            % Projection matrix for k = 0 and l = 1
            % Obtain intrinsic matrix
            intrinsic_l = camera.models.pinhole.IntrinsicMatrix( pixelSize ...
                                                               , principalPoint.obtainVectors(3) ...
                                                               , focalLength ...
                                                               , skewFactor );
                
            % Obtain the extrinsic matrix
            extrinsic_l = camera.models.ExtrinsicMatrix( rotation ...
                                                       , translation.data - projectionCenters.obtainVectors(3).data );
                
            % Create microlens pinhole camera
            self = camera.models.plenoptic.standard.MicrolensProjectionMatrix();
            self.intrinsicMatrix = intrinsic_0.intrinsicMatrix; 
            self.extrinsicMatrix = extrinsic_0.extrinsicMatrix;
            self.incrementIntrinsicMatrix_x = intrinsic_k.intrinsicMatrix ...
                                            - intrinsic_0.intrinsicMatrix;
            self.incrementIntrinsicMatrix_y = intrinsic_l.intrinsicMatrix ...
                                            - intrinsic_0.intrinsicMatrix;
            self.incrementExtrinsicMatrix_x = extrinsic_k.extrinsicMatrix ...
                                            - extrinsic_0.extrinsicMatrix;
            self.incrementExtrinsicMatrix_y = extrinsic_l.extrinsicMatrix ...
                                            - extrinsic_0.extrinsicMatrix;
        end
    end
end
