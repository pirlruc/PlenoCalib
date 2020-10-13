classdef Camera
    %CAMERA
    %   Standard plenoptic camera model. The intrinsic matrix for
    %   this class should be provided with a lightfield parameterized with
    %   one point and one direction.
    %
    %   Remember that the indices (i,j,k,l) are one-based indices.
    
    properties
        intrinsicMatrix  = camera.models.plenoptic.standard.IntrinsicMatrix().intrinsicMatrixDirection
                % Intrinsic matrix to convert rays from image to object space
        distortion       = camera.models.plenoptic.Distortion()
                % Distortion instance to remove distortion from rays in the object space
        extrinsic        = camera.models.ExtrinsicMatrix()
                % Extrinsic parameters for camera
        lightfieldSize   = lightfield.LightfieldSize([3,10,10,380,380])   
                % Lightfield size
    end
    
    properties (Dependent)
        % Intrinsic matrix entries
        h_si                    % Entry (1,1) of intrinsic matrix
        h_sk                    % Entry (1,3) of intrinsic matrix
        h_sl                    % Entry (1,4) of intrinsic matrix
        h_s                     % Entry (1,5) of intrinsic matrix
        h_tj                    % Entry (2,2) of intrinsic matrix
        h_tl                    % Entry (2,4) of intrinsic matrix
        h_t                     % Entry (2,5) of intrinsic matrix
        h_ui                    % Entry (3,1) of intrinsic matrix
        h_uk                    % Entry (3,3) of intrinsic matrix
        h_ul                    % Entry (3,4) of intrinsic matrix
        h_u                     % Entry (3,5) of intrinsic matrix
        h_vj                    % Entry (4,2) of intrinsic matrix
        h_vl                    % Entry (4,4) of intrinsic matrix
        h_v                     % Entry (4,5) of intrinsic matrix
        
        % Orientations
        du_di                   % Change in direction u for change in pixel i
        du_dk                   % Change in direction u for change in microlens k
        dv_dj                   % Change in direction v for change in pixel j
        dv_dl                   % Change in direction v for change in microlens l

        % Projection centers
        ds_di                   % Change in position s for change in pixel i
        ds_dk                   % Change in position s for change in microlens k
        dt_dj                   % Change in position t for change in pixel j
        dt_dl                   % Change in position t for change in microlens l
        di_du                   % Change in pixel i for change in direction u
        dk_du                   % Change in microlens k for change in direction u
        dj_dv                   % Change in pixel j for change in direction v
        dl_dv                   % Change in microlens l for change in direction v

        % Matrix blocks
        H_stij                  % H(1:2,1:2) entries
        H_stkl                  % H(1:2,3:4) entries
        H_uvij                  % H(3:4,1:2) entries
        H_uvkl                  % H(3:4,3:4) entries
        h_st                    % H(1:2,5) entries
        h_uv                    % H(3:4,5) entries
        
        % Image space projection
        singularities           % Derivative singularities
    end
    
    methods
        function self = Camera(varargin)
            %
            % Create standard plenoptic camera instance.
            % 
            % INPUTS: 
            %   1. intrinsicMatrix  - intrinsic matrix data.
            %   2. lightfieldSize   - lightfield size data.
            %   3. distortion       - distortion parameters and center of
            %   distortion.
            %   4. extrinsic        - extrinsic matrix data.
            %
            narginchk(0,4);

            if ~isempty(varargin)
                if nargin >= 4
                    self.extrinsic = varargin{4};
                end
                
                if nargin >= 3
                    self.distortion = varargin{3};
                end
                
                if nargin >= 2
                    self.lightfieldSize = varargin{2};
                end
                
                if nargin >= 1
                    self.intrinsicMatrix = varargin{1};
                end
            end
        end
        
        function self = set.intrinsicMatrix(self,newIntrinsicMatrixData)
            self.intrinsicMatrix = newIntrinsicMatrixData;
        end
        
        function self = set.lightfieldSize(self,newLightfieldSize)
            self.lightfieldSize = utils.enums.Classes.LIGHTFIELD_SIZE().convert(newLightfieldSize);
        end
        
        function self = set.distortion(self,newDistortion)
            self.distortion = newDistortion;
        end
        
        function self = set.extrinsic(self,newExtrinsic)
            self.extrinsic = newExtrinsic;
        end
        
        function h_si = get.h_si(self)
            h_si = self.intrinsicMatrix(1,1);
        end
        
        function h_sk = get.h_sk(self)
            h_sk = self.intrinsicMatrix(1,3);
        end
        
        function h_sl = get.h_sl(self)
            h_sl = self.intrinsicMatrix(1,4);
        end
        
        function h_s = get.h_s(self)
            h_s = self.intrinsicMatrix(1,5);
        end
        
        function h_tj = get.h_tj(self)
            h_tj = self.intrinsicMatrix(2,2);
        end
        
        function h_tl = get.h_tl(self)
            h_tl = self.intrinsicMatrix(2,4);
        end
        
        function h_t = get.h_t(self)
            h_t = self.intrinsicMatrix(2,5);
        end
        
        function h_ui = get.h_ui(self)
            h_ui = self.intrinsicMatrix(3,1);
        end
        
        function h_uk = get.h_uk(self)
            h_uk = self.intrinsicMatrix(3,3);
        end
        
        function h_ul = get.h_ul(self)
            h_ul = self.intrinsicMatrix(3,4);
        end
        
        function h_u = get.h_u(self)
            h_u = self.intrinsicMatrix(3,5);
        end
        
        function h_vj = get.h_vj(self)
            h_vj = self.intrinsicMatrix(4,2);
        end
        
        function h_vl = get.h_vl(self)
            h_vl = self.intrinsicMatrix(4,4);
        end
        
        function h_v = get.h_v(self)
            h_v = self.intrinsicMatrix(4,5);
        end
        
        function du_di = get.du_di(self)
            du_di = self.h_ui;
        end
        
        function du_dk = get.du_dk(self)
            du_dk = self.h_uk;
        end
        
        function dv_dj = get.dv_dj(self)
            dv_dj = self.h_vj;
        end
        
        function dv_dl = get.dv_dl(self)
            dv_dl = self.h_vl;
        end
        
        function ds_di = get.ds_di(self)
            ds_di = self.h_si;
        end
        
        function ds_dk = get.ds_dk(self)
            ds_dk = self.h_sk;
        end
        
        function dt_dj = get.dt_dj(self)
            dt_dj = self.h_tj;
        end
        
        function dt_dl = get.dt_dl(self)
            dt_dl = self.h_tl;
        end
        
        function di_du = get.di_du(self)
            di_du = 1/self.h_ui;
        end
        
        function dk_du = get.dk_du(self)
            dk_du = 1/self.h_uk;
        end
        
        function dj_dv = get.dj_dv(self)
            dj_dv = 1/self.h_vj;
        end
        
        function dl_dv = get.dl_dv(self)
            dl_dv = 1/self.h_vl;
        end

        function matrix = get.H_stij(self)
            matrix = self.intrinsicMatrix(1:2,1:2);
        end
        
        function matrix = get.H_stkl(self)
            matrix = self.intrinsicMatrix(1:2,3:4);
        end
        
        function matrix = get.H_uvij(self)
            matrix = self.intrinsicMatrix(3:4,1:2);
        end
        
        function matrix = get.H_uvkl(self)
            matrix = self.intrinsicMatrix(3:4,3:4);
        end
        
        function matrix = get.h_st(self)
            matrix = self.intrinsicMatrix(1:2,5);
        end
        
        function matrix = get.h_uv(self)
            matrix = self.intrinsicMatrix(3:4,5);
        end
        
        function singularities = get.singularities(self)
            % Remember that singularities occur when the denominator is
            % zero, so the zero of di_dk gives the singularity of dk_di.

            % Corresponds to the microlenses projection centers - points in
            % focus by the main lens
            singularities.ik = -self.h_si / self.h_ui;
            singularities.jl = -self.h_tj / self.h_vj;

            % Corresponds to the viewpoints projection centers
            singularities.ki = -self.h_sk / self.h_uk;
            singularities.lj = -self.h_tl / self.h_vl;
        end
    end
    
    methods
        function undistortedObjectRays_stuv = obtainObjectRaysFromImageRays( self ...
                                                                           , distortedImageRays_ijkl ...
                                                                           , removeDistortion )
            %
            % Obtains object space rays from image space rays.
            %
            % INPUTS:
            %   1. distortedImageRays_ijkl - distorted image space rays in 
            %      pixels and microlenses indices. Each ray should be 
            %      provided in different columns.
            %   2. removeDistortion - flag to indicate if distortion should 
            %      be removed from the image rays. Default is true.
            %
            narginchk(1,3);
            
            % If image rays are not provided, obtain a regular sampling of
            % rays
            if nargin <= 1
                [pixel_i,pixel_j,microlens_k,microlens_l] = ndgrid(1:4:10,1:4:10,1:20:100,1:20:100);
                distortedImageRays_ijkl = [pixel_i(:),pixel_j(:),microlens_k(:),microlens_l(:)]';
            end
            
            % If flag for distortion is not provided, consider that it is
            % true.
            if nargin <= 2
                removeDistortion = true;
            end
            
            distortedImageRays_ijkl = utils.enums.Classes.IMAGE_RAY().convert(distortedImageRays_ijkl);

            % Represent image rays with homogeneous coordinates
            distortedImageRays_ijkl = distortedImageRays_ijkl.setHomogeneousCoordinates();
                
            % Obtain rays in object space
            distortedObjectRays_stuv = lightfield.ObjectRay(   self.intrinsicMatrix ...
                                                             * distortedImageRays_ijkl.data ...
                                                           , true );
            
            % Remove homogeneous coordinates
            distortedObjectRays_stuv = distortedObjectRays_stuv.removeHomogeneousCoordinates();

            % Remove distortion from object rays
            if removeDistortion == true
                undistortedObjectRays_stuv = self.distortion.removeDistortion( distortedObjectRays_stuv );
            else
                undistortedObjectRays_stuv = distortedObjectRays_stuv;
            end
        end
        
        function projectionRays = projectUsingCameraIndices( self, worldPoint ...
                                                           , imageRayCoordinates )
            %
            % Obtain projection rays for a world point in specific image 
            % ray coordinates. This part does not consider distortion.
            %
            % INPUTS:
            %   1. worldPoint - point in the world coordinate system. The 
            %   point should not be given in homogeneous coordinates.
            %   2. imageRayCoordinates - image ray coordinates considered
            %   for projection.
            %
            narginchk(3,3);
            
            % Convert to object classes
            imageRayCoordinates = utils.enums.Classes.IMAGE_RAY().convert(imageRayCoordinates);
            worldPoint          = utils.enums.Classes.POINT().convert(worldPoint);
            
            % Define projection method
            if all(~isnan(imageRayCoordinates.i))
                projection_kl = true;
            else
                projection_kl = false;    
            end

            % Obtain microlens coordinates
            if projection_kl == true
                pixels_ij = [imageRayCoordinates.i;imageRayCoordinates.j];
                
                % Create viewpoint projection matrix
                viewpointProjection = camera.models.plenoptic.standard.ViewpointProjectionMatrix.ViewpointProjectionMatrixFromStandardPlenopticCamera(self);

                % Obtain the microlenses coordinates
                projectionRays = viewpointProjection.project(worldPoint,pixels_ij);

                % If there is distortion, display warning
                if any(self.distortion.encode ~= 0)
                    warning( 'StandardPlenopticCamera:distortionNotImplemented' ...
                           , 'Projection with distortion correction is not yet implemented' );
                end
                
            % Obtain pixel coordinates
            else
                microlenses_kl = [imageRayCoordinates.k;imageRayCoordinates.l];
            
                % Create microlens projection matrix
                microlensProjection = camera.models.plenoptic.standard.MicrolensProjectionMatrix.MicrolensProjectionMatrixFromStandardPlenopticCamera(self);

                % Obtain the microlenses coordinates
                projectionRays = microlensProjection.project(worldPoint,microlenses_kl);
                
                % If there is distortion, display warning
                if any(self.distortion.encode ~= 0)
                    warning( 'StandardPlenopticCamera:distortionNotImplemented' ...
                           , 'Projection with distortion correction is not yet implemented' );
                end
            end
        end

        function projectionRays = project( self, worldPoint ...
                                         , pixelsToNearestInteger ...
                                         , microlensesToNearestInteger ...
                                         , addDistortion ...
                                         , projectionType )
            %
            % Obtain projection rays for a world point in image rays using
            % rounding.
            %
            % INPUTS:
            %   1. worldPoint - point in the world coordinate system. The 
            %   point should not be given in homogeneous coordinates.
            %   2. pixelsToNearestInteger - flag that indicates if pixels
            %   must be rounded to the nearest integer. The default is 
            %   true.
            %   3. microlensesToNearestInteger - flag that indicates if 
            %   microlenses must be rounded to the nearest integer. The 
            %   default is true.
            %   4. addDistortion - flag to indicate if projection should
            %   consider distortion. Default is true.
            %   5. projectionType - flag to indicate if projection should
            %   be made on viewpoint or microlens images or a mixture.
            %   Default is a mixture projection since this allows to obtain
            %   more points in the lightfield.
            %
            % We are assuming that there is some symmetry regarding the 
            % pair of coordinates (i,k) and (j,l). If this symmetry does 
            % not occur display a warning.
            narginchk(1,6);
            
            % If no point is provided, generate random point.
            if nargin <= 1
                worldPoint = rand(3,1);
            end
            
            % If pixels round flag is not provided, assume that the pixels 
            % should be returned rounded to the nearest integer.
            if nargin <= 2
                pixelsToNearestInteger = true;
            end
            
            % If microlenses round flag is not provided, assume that the 
            % microlenses should be returned rounded to the nearest integer
            if nargin <= 3
                microlensesToNearestInteger = true;
            end
            
            % If flag for distortion is not provided, consider that it is
            % true.
            if nargin <= 4
                addDistortion = true;
            end
            
            if nargin <= 5
                projectionType = camera.lytro.enums.ProjectionTypes.MIXTURE();
            end
            
            % Number iterations considered to add distortion in the object
            % space
            NUMBER_ITERATIONS_CONVERGENCE = 10;
            
            % Define world points as points
            worldPoint = utils.enums.Classes.POINT().convert(worldPoint);
            
            % Obtain points in the camera coordinate system
            cameraPoint = self.extrinsic.obtainCameraPoints(worldPoint);
            
            % Obtain projection rays for camera point
            % Obtain derivative di/dk
            di_dk = - (self.h_sk + cameraPoint.z * self.h_uk) ...
                    / (self.h_si + cameraPoint.z * self.h_ui);
                
            % The change of microlens produces a change of more than one
            % pixel. Since more than one pixel can be filled for the same
            % microlens, we should project discretizing the pixels
            if projectionType == camera.lytro.enums.ProjectionTypes.VIEWPOINTS() ...
            || (  projectionType == camera.lytro.enums.ProjectionTypes.MIXTURE() ...
               && abs(di_dk) > 1 )
                pixels_i = 1:self.lightfieldSize.numberPixels_i;
                
                % Obtain the microlenses coordinates
                microlenses_k = - ( pixels_i .* (self.h_si + cameraPoint.z * self.h_ui) + self.h_s + cameraPoint.z * self.h_u - cameraPoint.x ) ...
                                / ( self.h_sk + cameraPoint.z * self.h_uk );
            
            % The change of microlens produces a change of less than one
            % pixel. Since the same pixel can be filled for different 
            % microlenses, we should project discretizing the microlenses
            else
                microlenses_k = 1:self.lightfieldSize.numberMicrolenses_k;
            
                % Obtain the pixel coordinates
                pixels_i = - ( microlenses_k .* (self.h_sk + cameraPoint.z * self.h_uk) + self.h_s + cameraPoint.z * self.h_u - cameraPoint.x ) ...
                           / ( self.h_si + cameraPoint.z * self.h_ui );
            end
            
            % Obtain derivative dj/dl
            dj_dl = - (self.h_tl + cameraPoint.z * self.h_vl) ...
                    / (self.h_tj + cameraPoint.z * self.h_vj);

            % The change of microlens produces a change of more than one
            % pixel. Since more than one pixel can be filled for the same
            % microlens, we should project discretizing the pixels
            if projectionType == camera.lytro.enums.ProjectionTypes.VIEWPOINTS() ...
            || (  projectionType == camera.lytro.enums.ProjectionTypes.MIXTURE() ...
               && abs(dj_dl) > 1 )
                pixels_j = 1:self.lightfieldSize.numberPixels_j;

                % Obtain the microlenses coordinates
                microlenses_l = - ( pixels_j .* (self.h_tj + cameraPoint.z * self.h_vj) + self.h_t + cameraPoint.z * self.h_v - cameraPoint.y ) ...
                                / ( self.h_tl + cameraPoint.z * self.h_vl );
                
            % The change of microlens produces a change of less than one
            % pixel. Since the same pixel can be filled for different 
            % microlenses, we should project discretizing the microlenses
            else
                microlenses_l = 1:self.lightfieldSize.numberMicrolenses_l;

                % Obtain the pixel coordinates
                pixels_j = - ( microlenses_l .* (self.h_tl + cameraPoint.z * self.h_vl) + self.h_t + cameraPoint.z * self.h_v - cameraPoint.y ) ...
                           / ( self.h_tj + cameraPoint.z * self.h_vj );
            end
            
            if addDistortion == true
                % Obtain combination of pixels and microlenses
                [pixels_i,pixels_j]           = meshgrid( pixels_i, pixels_j );
                [microlenses_k,microlenses_l] = meshgrid( microlenses_k, microlenses_l );

                % Obtain undistorted image rays in sensor
                undistortedImageRays_ijkl = lightfield.ImageRay( [ pixels_i(:)' ...
                                                                 ; pixels_j(:)' ...
                                                                 ; microlenses_k(:)' ...
                                                                 ; microlenses_l(:)' ], false );
                                                             
                % Obtain undistorted image rays in object space
                undistortedObjectRays_stuv = self.obtainObjectRaysFromImageRays( undistortedImageRays_ijkl ...
                                                                               , false );

                % Add distortion to image rays in object space
                distortedObjectRays_stuv = self.distortion.addDistortion( undistortedObjectRays_stuv ...
                                                                        , NUMBER_ITERATIONS_CONVERGENCE );
                
                % Obtain distorted image rays in sensor plane
                distortedObjectRays_stuv = distortedObjectRays_stuv.setHomogeneousCoordinates();
                distortedImageRays_ijkl  = lightfield.ImageRay( mldivide(self.intrinsicMatrix,distortedObjectRays_stuv.data), true);
                distortedImageRays_ijkl  = distortedImageRays_ijkl.removeHomogeneousCoordinates();
                
                % Obtain valid microlenses and pixels. The valid values for the
                % microlenses are 0.5 <= k < M + 0.5, assuming that we have
                % M microlenses. The valid values for the pixels are 
                % 0.5 <= i < N + 0.5, assuming that we have N pixels.
                validIndices_ijkl = distortedImageRays_ijkl.k >= 0.5 & distortedImageRays_ijkl.k < self.lightfieldSize.numberMicrolenses_k + 0.5 ...
                                  & distortedImageRays_ijkl.i >= 0.5 & distortedImageRays_ijkl.i < self.lightfieldSize.numberPixels_i + 0.5 ...
                                  & distortedImageRays_ijkl.l >= 0.5 & distortedImageRays_ijkl.l < self.lightfieldSize.numberMicrolenses_l + 0.5 ...
                                  & distortedImageRays_ijkl.j >= 0.5 & distortedImageRays_ijkl.j < self.lightfieldSize.numberPixels_j + 0.5;

                % Obtain valid lightfield image rays
                distortedImageRays_ijkl = distortedImageRays_ijkl.obtainVectors(validIndices_ijkl);
                pixels_i      = distortedImageRays_ijkl.i;
                pixels_j      = distortedImageRays_ijkl.j;
                microlenses_k = distortedImageRays_ijkl.k;
                microlenses_l = distortedImageRays_ijkl.l;
            else
                % Obtain valid microlenses and pixels. The valid values for the
                % microlenses are 0.5 <= k < M + 0.5, assuming that we have
                % M microlenses. The valid values for the pixels are 
                % 0.5 <= i < N + 0.5, assuming that we have N pixels.
                validIndices_ik = microlenses_k >= 0.5 & microlenses_k < self.lightfieldSize.numberMicrolenses_k + 0.5 ...
                                & pixels_i      >= 0.5 & pixels_i      < self.lightfieldSize.numberPixels_i + 0.5;
                validIndices_jl = microlenses_l >= 0.5 & microlenses_l < self.lightfieldSize.numberMicrolenses_l + 0.5 ...
                                & pixels_j      >= 0.5 & pixels_j      < self.lightfieldSize.numberPixels_j + 0.5;

                % Obtain valid pixels and microlenses
                [pixels_i,pixels_j]           = meshgrid( pixels_i(validIndices_ik), pixels_j(validIndices_jl));
                [microlenses_k,microlenses_l] = meshgrid( microlenses_k(validIndices_ik), microlenses_l(validIndices_jl) );
            end
            
            % Obtain projection rays
            projectionRays = lightfield.ImageRay([ pixels_i(:)'; pixels_j(:)'; microlenses_k(:)'; microlenses_l(:)' ],false);
            
            % Round projection ray coordinates
            % Round pixels and microlenses
            if     pixelsToNearestInteger == true  && microlensesToNearestInteger == true
                projectionRays.data = round(projectionRays.data);
                projectionRays.data = unique(projectionRays.data','rows')';

            % Round pixels
            elseif pixelsToNearestInteger == true  && microlensesToNearestInteger == false
                projectionRays.data = [ round([projectionRays.i;projectionRays.j]); [projectionRays.k;projectionRays.l] ];
                projectionRays.data = unique(projectionRays.data','rows')';
            
            % Round microlenses
            elseif pixelsToNearestInteger == false && microlensesToNearestInteger == true
                projectionRays.data = [ [projectionRays.i;projectionRays.j]; round([projectionRays.k;projectionRays.l]) ];
                projectionRays.data = unique(projectionRays.data','rows')';
            end
        end
        
       function worldPoint = reconstructUsingObjectRays(self, distortedObjectRays_stuv ...
                                                             , removeDistortion )
            %
            % Reconstruct camera point from correspondences on object rays.
            %
            % INPUTS:
            %   1. distortedObjectRays_stuv - distorted object rays 
            %      correspondences.
            %   2. removeDistortion - flag to indicate if directions should
            %      be corrected by removing distortion to rays in the 
            %      object space. Default is true.
            %
            narginchk(2,3);
            
            if nargin <= 2
                removeDistortion = true;
            end
            
            distortedObjectRays_stuv = utils.enums.Classes.OBJECT_RAY().convert(distortedObjectRays_stuv);
            
            % Remove distortion from object rays
            if removeDistortion == true
                undistortedObjectRays_stuv = self.distortion.removeDistortion( distortedObjectRays_stuv );
            else
                undistortedObjectRays_stuv = distortedObjectRays_stuv;
            end
            
            % Initialize observation matrix and vector
            observationMatrix = zeros(2 * undistortedObjectRays_stuv.numberVectors,3);
            observationVector = zeros(2 * undistortedObjectRays_stuv.numberVectors,1);
            
            % Define observation matrix
            % x-coordinate
            observationMatrix(1:undistortedObjectRays_stuv.numberVectors,1) = ones(undistortedObjectRays_stuv.numberVectors,1);
            observationMatrix(1:undistortedObjectRays_stuv.numberVectors,3) = -undistortedObjectRays_stuv.u;
            observationVector(1:undistortedObjectRays_stuv.numberVectors,1) =  undistortedObjectRays_stuv.s;
            
            % y-coordinate
            observationMatrix(undistortedObjectRays_stuv.numberVectors + 1:end,2) = ones(undistortedObjectRays_stuv.numberVectors,1);
            observationMatrix(undistortedObjectRays_stuv.numberVectors + 1:end,3) = -undistortedObjectRays_stuv.v;
            observationVector(undistortedObjectRays_stuv.numberVectors + 1:end,1) =  undistortedObjectRays_stuv.t;
            
            % Obtain camera point
            cameraPoint = mldivide(observationMatrix,observationVector);
            
            % Obtain world point
            worldPoint  = self.extrinsic.obtainWorldPoints(cameraPoint);
        end
        
        function [cameraPoint,rmse,sse,residuals] = reconstruct( self, distortedImageRays_ijkl ...
                                                                     , removeDistortion )
            %
            % Reconstruct camera point from correspondences on image rays.
            %
            % INPUTS:
            %   1. distortedImageRays_ijkl - distorted image rays 
            %      correspondences.
            %   2. removeDistortion - flag to indicate if directions should
            %      be corrected by removing distortion to rays in the 
            %      object space. Default is true.
            %
            narginchk(2,3);
            
            % If flag for distortion is not provided, consider that it is
            % true.
            if nargin <= 2
                removeDistortion = true;
            end
            
            distortedImageRays_ijkl = utils.enums.Classes.IMAGE_RAY().convert(distortedImageRays_ijkl);
            
            % Obtain object space rays. I can remove validations since the
            % inputs are mandatory.
            undistortedObjectRays_stuv = self.obtainObjectRaysFromImageRays( distortedImageRays_ijkl ...
                                                                           , removeDistortion );
                        
            % Camera point estimation
            if nargout <= 1
                cameraPoint = self.reconstructFromObjectRays(undistortedObjectRays_stuv,false);
            else
                [cameraPoint,rmse,sse,residuals] = self.reconstructFromObjectRays(undistortedObjectRays_stuv,false);
            end
        end
        
        function [cameraPoint,rmse,sse,residuals] = reconstructFromObjectRays( self, distortedObjectRays_stuv ...
                                                                                   , removeDistortion )
            %
            % Reconstruct camera point from correspondences on object rays.
            %
            % INPUTS:
            %   1. distortedObjectRays_stuv - distorted object rays 
            %      correspondences.
            %   2. removeDistortion - flag to indicate if directions should
            %      be corrected by removing distortion to rays in the 
            %      object space. Default is true.
            %
            narginchk(2,3);
            
            % If flag for distortion is not provided, consider that it is
            % true.
            if nargin <= 2
                removeDistortion = true;
            end
            
            distortedObjectRays_stuv = utils.enums.Classes.OBJECT_RAY().convert(distortedObjectRays_stuv);
            
            % Remove distortion from object rays
            if removeDistortion == true
                undistortedObjectRays_stuv = self.distortion.removeDistortion( distortedObjectRays_stuv );
            else
                undistortedObjectRays_stuv = distortedObjectRays_stuv;
            end
            
            % Camera point estimation
            % Obtain position and direction
            position_st  = [undistortedObjectRays_stuv.s;undistortedObjectRays_stuv.t];
            direction_uv = [undistortedObjectRays_stuv.u;undistortedObjectRays_stuv.v];

            % Construct vector and matrix for least squares
            st_vector = position_st(:);
            uv_matrix = [repmat(eye(2,2),distortedObjectRays_stuv.numberVectors,1), -direction_uv(:)];

            % Solve least squares problem
            cameraPoint = math.Point(mldivide(uv_matrix,st_vector),false);
            
            % Obtain residuals and errors for least squares problem
            if nargout > 1
                residuals = abs(uv_matrix * cameraPoint.data - st_vector);
                sse  = sum(residuals.^2);               % Sum of squared errors
                rmse = sqrt(sse / length(residuals));   % Root mean squared error
            end
        end
        
        function [cameraPoint,rmse,sse,residuals] = reconstructByImposingLines(self, distortedImageRays_ijkl ...
                                                                                   , removeDistortion ...
                                                                                   , fitEachEpi )
            %
            % Reconstruct camera point from correspondences on image rays.
            % Since the image ray coordinates can be integers due to the
            % discretization performed during the lightfield acquisition,
            % we fit a line to the pairs (i,k) and (j,l) to obtain a higher
            % precision for the image ray coordinates. This reconstruction
            % uses point reconstruction after correcting the observations
            % to fit a line in each pair (i,k) and (j,l).
            %
            % INPUTS:
            %   1. imageRays_ijkl - distorted image rays correspondences. 
            %      Each correspondence should be given in different 
            %      columns.
            %   2. removeDistortion - flag to indicate if directions should
            %      be corrected by removing distortion to rays in the 
            %      object space. Default is true.
            %   3. fitEachEpi     - fit each line to an EPI like image.
            %      For example, this means that if true the fitting for 
            %      (i,k) is performed for each different value of (j,l). If
            %      false, only one fitting is performed. The default is
            %      false.
            %
            narginchk(2,4);
            
            % If distortion flag is not provided, assume default.
            if nargin <= 2
                removeDistortion = true;
            end
            
            % If fit each epi flag is not provided, assume default.
            if nargin <= 3
                fitEachEpi = false;
            end

            % Define image rays as image ray
            distortedImageRays_ijkl = utils.enums.Classes.IMAGE_RAY().convert(distortedImageRays_ijkl);

            % Remove distortion from image rays
            if removeDistortion == true
                % Obtain undistorted image rays in object space
                undistortedObjectRays_stuv = self.obtainObjectRaysFromImageRays( distortedImageRays_ijkl ...
                                                                               , removeDistortion );
                                                           
                % Obtain undistorted image rays in sensor
                undistortedObjectRays_stuv = undistortedObjectRays_stuv.setHomogeneousCoordinates();
                undistortedImageRays_ijkl  = lightfield.ImageRay( mldivide(self.intrinsicMatrix,undistortedObjectRays_stuv.data), true);
                undistortedImageRays_ijkl  = undistortedImageRays_ijkl.removeHomogeneousCoordinates();
            else
                undistortedImageRays_ijkl = distortedImageRays_ijkl;
            end
            
            %
            % Fit unique line for every observation of (i,k) and (j,l)
            %
            if fitEachEpi == false
                % Remove repeated observations to avoid memory problems in
                % singular value decomposition
                epiRays_ik = [undistortedImageRays_ijkl.i',undistortedImageRays_ijkl.k'];
                epiRays_ik = unique(epiRays_ik,'rows');
                epiRays_jl = [undistortedImageRays_ijkl.j',undistortedImageRays_ijkl.l'];
                epiRays_jl = unique(epiRays_jl,'rows');
                
                % Obtain line parameters for (i,k) and (j,l)
                lineFitting_ik = math.LineFitting(epiRays_ik(:,2),epiRays_ik(:,1));
                lineFitting_jl = math.LineFitting(epiRays_jl(:,2),epiRays_jl(:,1));
                lineParameters_ik = lineFitting_ik.fit();
                lineParameters_jl = lineFitting_jl.fit();

                % The estimation is done considering:
                %       a_hat * y + m_hat * x + b_hat = 0
                % 
                % Thus, we have to change the parameters obtained. Similarly to
                % what we have done in the project method, here we evaluate the
                % derivative di/dk to determine if we rewrite the parameters as
                % a funtion of k(i) or i(k).

                % Obtain derivative values
                di_dk = lineParameters_ik.m / lineParameters_ik.a;
                dj_dl = lineParameters_jl.m / lineParameters_jl.a;

                % If di/dk is greater than one, we write k(i). So we have to
                % transform the parameters to have something like:
                %       x = (1/m_hat) * (-a_hat * y - b1_hat).
                if abs(di_dk) > 1
                    lineParameters_ik = math.LineFitting.convertLineParameters_x(lineParameters_ik);

                    % Obtain new image ray coordinates
                    i = undistortedImageRays_ijkl.i;
                    k = lineParameters_ik.a .* undistortedImageRays_ijkl.i + lineParameters_ik.b;

                % If di/dk is less than one, we write i(k). So we have to
                % transform the parameters to have something like:
                %       y = (1/a_hat) * (-m_hat * x - b1_hat).
                else
                    lineParameters_ik = math.LineFitting.convertLineParameters_y(lineParameters_ik);

                    % Obtain new image ray coordinates
                    i = lineParameters_ik.m .* undistortedImageRays_ijkl.k + lineParameters_ik.b;
                    k = undistortedImageRays_ijkl.k;
                end

                % If dj/dl is greater than one, we write l(j). So we have to
                % transform the parameters to have something like:
                %       x = (1/m_hat) * (-a_hat * y - b1_hat).
                if abs(dj_dl) > 1
                    lineParameters_jl = math.LineFitting.convertLineParameters_x(lineParameters_jl);

                    % Obtain new image ray coordinates
                    j = undistortedImageRays_ijkl.j;
                    l = lineParameters_jl.a .* undistortedImageRays_ijkl.j + lineParameters_jl.b;

                % If dj/dl is less than one, we write j(l). So we have to
                % transform the parameters to have something like:
                %       y = (1/a_hat) * (-m_hat * x - b1_hat).
                else
                    lineParameters_jl = math.LineFitting.convertLineParameters_y(lineParameters_jl);

                    % Obtain new image ray coordinates
                    j = lineParameters_jl.m .* undistortedImageRays_ijkl.l + lineParameters_jl.b;
                    l = undistortedImageRays_ijkl.l;
                end

                % Obtain new image ray coordinates
                newImageRays_ijkl = [i;j;k;l];
                
            %
            % Fit line to each EPI. Apply this process to (i,k) and (j,l)
            %
            else
                % Round projection rays to obtain indices
                indices_ijkl = lightfield.ImageRay(round(undistortedImageRays_ijkl.data),false);

                % Obtain line parameters for (i,k). 
                % The line parameters should be obtained in the EPI like
                % images. For example, this means that the fitting for (i,k) is
                % performed for each combination of (j,l) values.
                newImageRays_ijkl = [];
                for pixel_j = 1:self.lightfieldSize.numberPixels_j
                    for microlens_l = 1:self.lightfieldSize.numberMicrolenses_l
                        % Obtain indices to select (i,k) observations for (j,l)
                        validPixels         = indices_ijkl.j == pixel_j;
                        validMicrolenses    = indices_ijkl.l == microlens_l;
                        validImageRays_ijkl = validPixels & validMicrolenses;

                        % If there is no projections for a particular
                        % coordinate pair (j,l) skip fitting
                        if any(validImageRays_ijkl) == false
                            continue
                        end
                        
                        % Obtain coordinates (j,l)
                        j = undistortedImageRays_ijkl.j(validImageRays_ijkl);
                        l = undistortedImageRays_ijkl.l(validImageRays_ijkl);

                        % Obtain line parameters for (i,k)
                        lineFitting_ik    = math.LineFitting( undistortedImageRays_ijkl.k(validImageRays_ijkl) ...
                                                            , undistortedImageRays_ijkl.i(validImageRays_ijkl) );
                        lineParameters_ik = lineFitting_ik.fit();

                        % Correct the projection rays by using the line
                        % parameters.
                        %
                        % The estimation is done considering:
                        %       a_hat * y + m_hat * x + b_hat = 0
                        % 
                        % Thus, we have to change the parameters obtained. Similarly to
                        % what we have done in the project method, here we evaluate the
                        % derivative di/dk to determine if we rewrite the parameters as
                        % a funtion of k(i) or i(k).

                        % Obtain derivative values
                        di_dk = lineParameters_ik.m / lineParameters_ik.a;

                        % If di/dk is greater than one, we write k(i). So we have to
                        % transform the parameters to have something like:
                        %       x = (1/m_hat) * (-a_hat * y - b1_hat).
                        if abs(di_dk) > 1
                            lineParameters_ik = math.LineFitting.convertLineParameters_x(lineParameters_ik);

                            % Obtain new image ray coordinates
                            i = undistortedImageRays_ijkl.i(validImageRays_ijkl);
                            k = lineParameters_ik.a .* i + lineParameters_ik.b;

                        % If di/dk is less than one, we write i(k). So we have to
                        % transform the parameters to have something like:
                        %       y = (1/a_hat) * (-m_hat * x - b1_hat).
                        else
                            lineParameters_ik = math.LineFitting.convertLineParameters_y(lineParameters_ik);

                            % Obtain new image ray coordinates
                            k = undistortedImageRays_ijkl.k(validImageRays_ijkl);
                            i = lineParameters_ik.m .* k + lineParameters_ik.b;
                        end

                        % Add new image rays
                        newImageRays_ijkl = cat(2,newImageRays_ijkl,[i;j;k;l]);
                    end
                end

                % Update image rays
                distortedImageRays_ijkl.data = newImageRays_ijkl;

                % Obtain line parameters for (j,l) applying the same strategy
                % applied to (i,k).
                newImageRays_ijkl = [];
                for pixel_i = 1:self.lightfieldSize.numberPixels_i
                    for microlens_k = 1:self.lightfieldSize.numberMicrolenses_k
                        % Obtain indices to select (j,l) observations for (i,k)
                        validPixels         = indices_ijkl.i == pixel_i;
                        validMicrolenses    = indices_ijkl.k == microlens_k;
                        validImageRays_ijkl = validPixels & validMicrolenses;

                        % If there is no projections for a particular
                        % coordinate pair (i,k) skip fitting
                        if any(validImageRays_ijkl) == false
                            continue
                        end
                        
                        % Obtain coordinates (i,k)
                        i = undistortedImageRays_ijkl.i(validImageRays_ijkl);
                        k = undistortedImageRays_ijkl.k(validImageRays_ijkl);

                        % Obtain line parameters for (j,l)
                        lineFitting_jl    = math.LineFitting( undistortedImageRays_ijkl.l(validImageRays_ijkl) ...
                                                            , undistortedImageRays_ijkl.j(validImageRays_ijkl) );
                        lineParameters_jl = lineFitting_jl.fit();

                        % Correct the projection rays by using the line
                        % parameters applying the same strategy of (i,k).
                        % Obtain derivative values
                        dj_dl = lineParameters_jl.m / lineParameters_jl.a;

                        % If di/dk is greater than one, we write k(i). So we have to
                        % transform the parameters to have something like:
                        %       x = (1/m_hat) * (-a_hat * y - b1_hat).
                        if abs(dj_dl) > 1
                            lineParameters_jl = math.LineFitting.convertLineParameters_x(lineParameters_jl);

                            % Obtain new image ray coordinates
                            j = undistortedImageRays_ijkl.j(validImageRays_ijkl);
                            l = lineParameters_jl.a .* j + lineParameters_jl.b;

                        % If di/dk is less than one, we write i(k). So we have to
                        % transform the parameters to have something like:
                        %       y = (1/a_hat) * (-m_hat * x - b1_hat).
                        else
                            lineParameters_jl = math.LineFitting.convertLineParameters_y(lineParameters_jl);

                            % Obtain new image ray coordinates
                            l = imageRays_ijkll(validImageRays_ijkl);
                            j = lineParameters_jl.m .* l + lineParameters_jl.b;
                        end

                        % Add new image rays
                        newImageRays_ijkl = cat(2,newImageRays_ijkl,[i;j;k;l]);
                    end
                end
            end
            
            % Reconstruct points and obtain reconstruction errors.
            % Reconstruction will not remove distortion since it was
            % already removed during the analysis.
            REMOVE_DISTORTION = false;
            if nargout == 1
                cameraPoint = self.reconstruct(newImageRays_ijkl,REMOVE_DISTORTION);
            else
                [cameraPoint, rmse, sse, residuals] = self.reconstruct( newImageRays_ijkl ...
                                                                      , REMOVE_DISTORTION );
            end
        end
        
        function [cameraPoint,rmse,sse,residuals] = reconstructUsingLineParameters( self ...
                                                                                  , distortedImageRays_ijkl ...
                                                                                  , removeDistortion ...
                                                                                  , fitEachEpi)
            %
            % Reconstruct camera point from correspondences on image rays.
            % Since the image ray coordinates can be integers due to the
            % discretization performed during the lightfield acquisition,
            % we fit a line to the pairs (i,k) and (j,l) to obtain a higher
            % precision for the image ray coordinates. This method
            % reconstructs the camera point by looking into the relation
            % between the camera point and the parameters of the lines
            % fitted to the pairs (i,k) and (j,l).
            %
            % INPUTS:
            %   1. distortedImageRays_ijkl - distorted image rays 
            %      correspondences. Each correspondence should be given in 
            %      different columns.
            %   2. removeDistortion - flag to indicate if directions should
            %      be corrected by remove distortion to rays in the object
            %      space. Default is true.
            %   3. fitEachEpi     - fit each line to an EPI like image.
            %      For example, this means that if true the fitting for 
            %      (i,k) is performed for each different value of (j,l). If
            %      false, only one fitting is performed. The default is
            %      false.
            %
            narginchk(2,4);
            
            NUMBER_LINE_PARAMETERS = 3;

            % If flag is not provided, assume default.
            if nargin <= 2
                removeDistortion = true;
            end
            
            % If flag is not provided, assume default.
            if nargin <= 3
                fitEachEpi = false;
            end
            
            % Define image rays as image ray
            distortedImageRays_ijkl = utils.enums.Classes.IMAGE_RAY().convert(distortedImageRays_ijkl);

            % Remove distortion from image rays
            if removeDistortion == true
                % Obtain undistorted image rays in object space
                undistortedObjectRays_stuv = self.obtainObjectRaysFromImageRays( distortedImageRays_ijkl ...
                                                                               , removeDistortion );
                                                           
                % Obtain undistorted image rays in sensor
                undistortedObjectRays_stuv = undistortedObjectRays_stuv.setHomogeneousCoordinates();
                undistortedImageRays_ijkl  = lightfield.ImageRay( mldivide(self.intrinsicMatrix,undistortedObjectRays_stuv.data), true);
                undistortedImageRays_ijkl  = undistortedImageRays_ijkl.removeHomogeneousCoordinates();
            else
                undistortedImageRays_ijkl = distortedImageRays_ijkl;
            end
            
            %
            % Fit unique line for every observation of (i,k) and (j,l)
            %
            if fitEachEpi == false
                % Remove repeated observations to avoid memory problems in
                % singular value decomposition
                epiRays_ik = [undistortedImageRays_ijkl.i',undistortedImageRays_ijkl.k'];
                epiRays_ik = unique(epiRays_ik,'rows');
                epiRays_jl = [undistortedImageRays_ijkl.j',undistortedImageRays_ijkl.l'];
                epiRays_jl = unique(epiRays_jl,'rows');
                
                % Obtain line parameters for (i,k) and (j,l)
                lineFitting_ik = math.LineFitting(epiRays_ik(:,2),epiRays_ik(:,1));
                lineFitting_jl = math.LineFitting(epiRays_jl(:,2),epiRays_jl(:,1));
                lineParameters_ik = lineFitting_ik.fit();
                lineParameters_jl = lineFitting_jl.fit();

                % The estimation is done considering:
                %       ai_hat * i + mi_hat * k + bi_hat = 0
                %       aj_hat * j + mj_hat * l + bj_hat = 0
                % 
                % Remember that the line parameters are obtained
                % considering that the norm is equal to one. Thus, these 
                % parameters are related with the scene point m = [x,y,z]
                % by the following equations:
                %       lambda_ik * ai_hat = h_si + z * h_ui
                %       lambda_jl * aj_hat = h_tj + z * h_vj
                %       lambda_ik * mi_hat = h_sk + z * h_uk
                %       lambda_jl * mj_hat = h_tl + z * h_vl
                %       lambda_ik * bi_hat = h_s + z * h_u - x
                %       lambda_jl * bj_hat = h_t + z * h_v - y

                % Obtain the observation matrix
                observationMatrix = [  0,  0, self.h_ui, lineParameters_ik.a, 0 ...
                                    ;  0,  0, self.h_vj, 0, lineParameters_jl.a ...
                                    ;  0,  0, self.h_uk, lineParameters_ik.m, 0 ...
                                    ;  0,  0, self.h_vl, 0, lineParameters_jl.m ...
                                    ; -1,  0, self.h_u , lineParameters_ik.b, 0 ...
                                    ;  0, -1, self.h_v , 0, lineParameters_jl.b ];

                % Obtain the observation vector
                observationVector = [ -self.h_si ...
                                    ; -self.h_tj ...
                                    ; -self.h_sk ...
                                    ; -self.h_tl ...
                                    ; -self.h_s ...
                                    ; -self.h_t ];
                                
            %
            % Fit line to each EPI. Apply this process to (i,k) and (j,l)
            %
            else
                % Round projection rays to obtain indices
                indices_ijkl = lightfield.ImageRay(round(undistortedImageRays_ijkl.data),false);

                % The estimation is done considering:
                %       ai_hat * i + mi_hat * k + bi_hat = 0
                %       aj_hat * j + mj_hat * l + bj_hat = 0
                % 
                % Remember that the line parameters are obtained
                % considering that the norm is equal to one. Thus, these 
                % parameters are related with the scene point m = [x,y,z]
                % by the following equations:
                %       lambda_ik * ai_hat = h_si + z * h_ui
                %       lambda_jl * aj_hat = h_tj + z * h_vj
                %       lambda_ik * mi_hat = h_sk + z * h_uk
                %       lambda_jl * mj_hat = h_tl + z * h_vl
                %       lambda_ik * bi_hat = h_s + z * h_u - x
                %       lambda_jl * bj_hat = h_t + z * h_v - y

                % Define the observation matrices
                observationMatrixBlock_ik = [  0,  0, self.h_ui ...
                                            ;  0,  0, self.h_uk ...
                                            ; -1,  0, self.h_u  ];
                observationMatrixBlock_jl = [  0,  0, self.h_vj ...
                                            ;  0,  0, self.h_vl ...
                                            ;  0, -1, self.h_v  ];
                                        
                % Define the observation vectors
                observationVector_ik = [ -self.h_si ...
                                       ; -self.h_sk ...
                                       ; -self.h_s ];

                observationVector_jl = [ -self.h_tj ...
                                       ; -self.h_tl ...
                                       ; -self.h_t ];

                % Obtain line parameters for (i,k). 
                % The line parameters should be obtained in the EPI like
                % images. For example, this means that the fitting for (i,k) is
                % performed for each combination of (j,l) values.
                observationMatrix = [];
                observationVector = [];
                lineParameters    = [];
                fittedLines = 0;
                for pixel_j = 1:self.lightfieldSize.numberPixels_j
                    for microlens_l = 1:self.lightfieldSize.numberMicrolenses_l
                        % Obtain indices to select (i,k) observations for (j,l)
                        validPixels         = indices_ijkl.j == pixel_j;
                        validMicrolenses    = indices_ijkl.l == microlens_l;
                        validImageRays_ijkl = validPixels & validMicrolenses;

                        % If there is no projections for a particular
                        % coordinate pair (j,l) skip fitting
                        if any(validImageRays_ijkl) == false
                            continue
                        end
                        
                        % Obtain line parameters for (i,k)
                        fittedLines = fittedLines + 1;
                        lineFitting_ik    = math.LineFitting( undistortedImageRays_ijkl.k(validImageRays_ijkl) ...
                                                            , undistortedImageRays_ijkl.i(validImageRays_ijkl) );
                        lineParameters_ik = lineFitting_ik.fit();

                        % Obtain the line parameters
                        lineParameters = cat(1,lineParameters, [ lineParameters_ik.a ...
                                                               ; lineParameters_ik.m ...
                                                               ; lineParameters_ik.b ] );

                        % Add observations to observation vector
                        observationVector = cat(1,observationVector,observationVector_ik);
                        
                        % Add block to observation matrix
                        observationMatrix = cat(1,observationMatrix,observationMatrixBlock_ik);
                    end
                end

                % Obtain line parameters for (j,l) applying the same strategy
                % for (i,k).
                for pixel_i = 1:self.lightfieldSize.numberPixels_i
                    for microlens_k = 1:self.lightfieldSize.numberMicrolenses_k
                        % Obtain indices to select (j,l) observations for (i,k)
                        validPixels         = indices_ijkl.i == pixel_i;
                        validMicrolenses    = indices_ijkl.k == microlens_k;
                        validImageRays_ijkl = validPixels & validMicrolenses;

                        % If there is no projections for a particular
                        % coordinate pair (i,k) skip fitting
                        if any(validImageRays_ijkl) == false
                            continue
                        end
                        
                        % Obtain line parameters for (j,l)
                        fittedLines = fittedLines + 1;
                        lineFitting_jl    = math.LineFitting( undistortedImageRays_ijkl.l(validImageRays_ijkl) ...
                                                            , undistortedImageRays_ijkl.j(validImageRays_ijkl) );
                        lineParameters_jl = lineFitting_jl.fit();

                        % Obtain the line parameters
                        lineParameters = cat(1,lineParameters, [ lineParameters_jl.a ...
                                                               ; lineParameters_jl.m ...
                                                               ; lineParameters_jl.b ] );

                        % Add observations to observation vector
                        observationVector = cat(1,observationVector,observationVector_jl);
                        
                        % Add block to observation matrix
                        observationMatrix = cat(1,observationMatrix,observationMatrixBlock_jl);
                    end
                end
                
                % Add line parameters to observation matrix
                numberEquations = size(observationMatrix,1);
                numberVariables = size(observationMatrix,2);
                
                % Initialize new observation matrix
                observationMatrix = [observationMatrix,zeros(numberEquations,fittedLines)];
                parameterIndex = 0;
                % Set line parameters in observation matrix
                for iParameter = 1:fittedLines
                    observationMatrix( parameterIndex + 1: parameterIndex + NUMBER_LINE_PARAMETERS ...
                                     , numberVariables + iParameter ) = ...
                            lineParameters(parameterIndex + 1: parameterIndex + NUMBER_LINE_PARAMETERS);
                    parameterIndex = parameterIndex + NUMBER_LINE_PARAMETERS;
                end
            end
        
            % Solve least squares problem
            leastSquaresSolution = mldivide(observationMatrix,observationVector);
            
            % The least squares solution gives the point and the scale
            % factors used in the fitting of the lines. We are only
            % interested on the point.
            cameraPoint = math.Point(leastSquaresSolution(1:3));
            
            % Obtain residuals and errors for least squares problem
            if nargout > 1
                residuals = abs(observationMatrix * leastSquaresSolution - observationVector);
                sse  = sum(residuals.^2);               % Sum of squared errors
                rmse = sqrt(sse / length(residuals));   % Root mean squared error
            end
        end
        
        function cameraPoint = reconstructFromObjectRaysUsingLineParameters( self ...
                                                                           , distortedObjectRays_stuv ...
                                                                           , removeDistortion )
            %
            % Reconstruct camera point from correspondences on object rays.
            % This method reconstructs the camera point by looking into the 
            % relation between the camera point and the parameters of the 
            % lines fitted to the pairs (s,u) and (t,v).
            %
            % INPUTS:
            %   1. distortedObjectRays_stuv - distorted object rays 
            %      correspondences. Each correspondence should be given in 
            %      different columns.
            %   2. removeDistortion - flag to indicate if directions should
            %      be corrected by remove distortion to rays in the object
            %      space. Default is true.
            %
            narginchk(2,3);
            
            % If flag is not provided, assume default.
            if nargin <= 2
                removeDistortion = true;
            end
            
            % Define image rays as image ray
            distortedObjectRays_stuv = utils.enums.Classes.OBJECT_RAY().convert(distortedObjectRays_stuv);

            % Remove distortion from object rays
            if removeDistortion == true
                undistortedObjectRays_stuv = self.distortion.removeDistortion( distortedObjectRays_stuv );
            else
                undistortedObjectRays_stuv = distortedObjectRays_stuv;
            end
            
            % Remove repeated observations to avoid memory problems in
            % singular value decomposition
            epiRays_su = [undistortedObjectRays_stuv.s',undistortedObjectRays_stuv.u'];
            epiRays_su = unique(epiRays_su,'rows');
            epiRays_tv = [undistortedObjectRays_stuv.t',undistortedObjectRays_stuv.v'];
            epiRays_tv = unique(epiRays_tv,'rows');

            % Obtain line parameters for (s,u) and (t,v)
            lineFitting_su = math.LineFitting(epiRays_su(:,2),epiRays_su(:,1));
            lineFitting_tv = math.LineFitting(epiRays_tv(:,2),epiRays_tv(:,1));
            lineParameters_su = lineFitting_su.fit();
            lineParameters_tv = lineFitting_tv.fit();

            % The estimation is done considering:
            %       as_hat * s + ms_hat * u + bs_hat = 0
            %       at_hat * t + mt_hat * v + bt_hat = 0
            % 
            % Remember that the line parameters are obtained
            % considering that the norm is equal to one. Thus, these 
            % parameters are related with the scene point m = [x,y,z]
            % by the following equations:
            %       lambda_su = as_hat
            %       lambda_tv = at_hat
            %       z = ms_hat / lambda_su
            %       z = mt_hat / lambda_tv
            %       x = -bs_hat / lambda_su
            %       y = -bt_hat / lambda_tv

            % Obtain the scale factor of the line equations
            lambda_su = lineParameters_su.a;
            lambda_tv = lineParameters_tv.a;
            
            % Obtain (x,y,z) coordinates
            point_z = mean( [ lineParameters_su.m / lambda_su ...
                            ; lineParameters_tv.m / lambda_tv ] );
            point_x = -lineParameters_su.b / lambda_su;
            point_y = -lineParameters_tv.b / lambda_tv;

            % Define camera point
            cameraPoint = math.Point([point_x;point_y;point_z],false);
        end
        
         function pointsInFocus = obtainPointsInFocusByMicrolenses(self,microlenses_kl)
            %
            % Obtain points in focus by a given set of microlens cameras.
            % The points in focus by the microlens camera corresponds to a
            % point in the focal plane of the main lens.
            %
            % INPUTS:
            %   1. microlenses_kl - microlenses indices to obtain
            %      points in focus. Each microlens should be provided in 
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
            
            % Obtain points in focus by the microlenses
            pointsInFocus = camera.models.plenoptic.standard.MicrolensProjectionMatrix.obtainProjectionCentersFromStandardPlenopticCamera( self ...
                                                                                                                                         , microlenses_kl );
        end

        function lfield = rectifyLightfield(self,lfield)
            %
            % Rectify lightfield by removing distortion.
            %
            % INPUTS:
            %   1. lightfieldData  - lightfield data.
            %  
            narginchk(2,2);
            
            % Convert lightfield and intrinsic matrix
            lfield  = utils.enums.Classes.LIGHTFIELD().convert(lfield);

            % Obtain indices for lightfield
            pixels_i      = (1:lfield.size.numberPixels_i) - 0.5;
            pixels_j      = (1:lfield.size.numberPixels_j) - 0.5;
            microlenses_k = (1:lfield.size.numberMicrolenses_k) - 0.5;
            microlenses_l = (1:lfield.size.numberMicrolenses_l) - 0.5;
            [pixels_i,pixels_j,microlenses_k,microlenses_l] = ndgrid( pixels_i ...
                                                                    , pixels_j ...
                                                                    , microlenses_k ...
                                                                    , microlenses_l );

            % Obtain image rays
            undistortedImageRays_ijkl  = lightfield.ImageRay([pixels_i(:),pixels_j(:),microlenses_k(:),microlenses_l(:)], false);

            % Obtain object rays
            DO_NOT_REMOVE_DISTORTION   = false;
            undistortedObjectRays_stuv = self.obtainObjectRaysFromImageRays( undistortedImageRays_ijkl ...
                                                                           , DO_NOT_REMOVE_DISTORTION );
            
            % Add distortion to rays
            distortedObjectRays_stuv   = self.distortion.addDistortion(undistortedObjectRays_stuv,2);
            
            % Obtain image rays from distorted object rays
            distortedObjectRays_stuv   = distortedObjectRays_stuv.setHomogeneousCoordinates();
            distortedImageRays_ijkl    = lightfield.ImageRay(mldivide(self.intrinsicMatrix,distortedObjectRays_stuv.data),true);
            distortedImageRays_ijkl    = distortedImageRays_ijkl.removeHomogeneousCoordinates();

            % Remap intensity values
            for iChannel = 1:lfield.size.numberChannels
                % Obtain channel
                lightfieldChannel = permute(lfield.data(iChannel,:,:,:,:),[2,3,4,5,1]);
                
                % Interpolate information of lightfield
                lightfieldChannel = interpn( pixels_i ...
                                           , pixels_j ...
                                           , microlenses_k ...
                                           , microlenses_l ...
                                           , lightfieldChannel ...
                                           , distortedImageRays_ijkl.i ...
                                           , distortedImageRays_ijkl.j ...
                                           , distortedImageRays_ijkl.k ...
                                           , distortedImageRays_ijkl.l ...
                                           , 'cubic', nan );

                lightfieldChannel = reshape( lightfieldChannel ...
                                           , lfield.size.numberPixels_i ...
                                           , lfield.size.numberPixels_j ...
                                           , lfield.size.numberMicrolenses_k ...
                                           , lfield.size.numberMicrolenses_l );
                                       
                % Update lightfield data for channel
                lfield.data(iChannel,:,:,:,:) = lightfieldChannel;
            end

            % Clip interpolation result
            lfield.data(isnan(lfield.data)) = 0;
            lfield.data = max(0,min(1,lfield.data));
        end
    end
    
    methods (Static)
        function self = CameraFromCalibrationFile(calibrationFilepath)
            %
            % Create standard plenoptic camera instance from calibration
            % file.
            %
            % INPUTS:
            %   1. calibrationFilepath - path to calibration file.
            %
            narginchk(1,1);
            
            % Obtain calibration parameters used to define the standard
            % plenoptic camera
            [intrinsic,~,distortion,lightfieldSize] = ...
                camera.lytro.Camera.obtainCalibrationParametersFromCalibrationFile(calibrationFilepath);
            
            % Create standard plenoptic camera instance
            self = camera.models.plenoptic.standard.Camera( intrinsic ...
                                                          , lightfieldSize ...
                                                          , distortion );
        end
        
        function plenopticCameras = CameraFromMicrolensProjectionMatrix( microlensCameraArrays ...
                                                                       , lightfieldSize ...
                                                                       , cameraCoordinateOriginRay )
            % 
            % Create standard plenoptic camera model instance using
            % microlens camera array projection matrices.
            %
            % INPUTS:
            %   1. microlensCameraArrays - microlens camera arrays 
            %   projection matrices.
            %   2. lightfieldSize   - lightfield size. 
            %   3. cameraCoordinateOriginRay - ray constraint for the
            %   origin of the camera coordinate system.
            %
            narginchk(2,3);
            
            % Transform to lightfield size
            lightfieldSize = utils.enums.Classes.LIGHTFIELD_SIZE().convert(lightfieldSize);

            if nargin <= 2
                % Obtain central ray
                cameraCoordinateOriginRay = lightfield.ImageRay([ lightfieldSize.numberPixels_i - 1 ...
                                                                , lightfieldSize.numberPixels_j - 1 ...
                                                                , lightfieldSize.numberMicrolenses_k - 1 ...
                                                                , lightfieldSize.numberMicrolenses_l - 1 ] ./ 2 + 1,false);
            end

            % Transform to image ray
            cameraCoordinateOriginRay = utils.enums.Classes.IMAGE_RAY().convert(cameraCoordinateOriginRay);

            % For each camera array projection matrix, create intrinsic 
            % matrix for plenoptic camera
            numberMicrolensArrays = length(microlensCameraArrays);
			h_sk = zeros(1,numberMicrolensArrays);
			h_sl = zeros(1,numberMicrolensArrays);
			h_s  = zeros(1,numberMicrolensArrays);
			h_tl = zeros(1,numberMicrolensArrays);
			h_t  = zeros(1,numberMicrolensArrays);
			h_ui = zeros(1,numberMicrolensArrays);
			h_uk = zeros(1,numberMicrolensArrays);
			h_ul = zeros(1,numberMicrolensArrays);
			h_u  = zeros(1,numberMicrolensArrays);
			h_vj = zeros(1,numberMicrolensArrays);
			h_vl = zeros(1,numberMicrolensArrays);
			h_v  = zeros(1,numberMicrolensArrays);
			translationVectors = zeros(3,numberMicrolensArrays);
			for iMicrolensArray = 1:numberMicrolensArrays
				% Obtain camera array projection matrix
				projectionMatrix = microlensCameraArrays(iMicrolensArray);

				% Obtain information regarding directions (u,v)
				h_ui(iMicrolensArray) = 1 / projectionMatrix.intrinsicMatrix(1,1);
				h_vj(iMicrolensArray) = 1 / projectionMatrix.intrinsicMatrix(2,2);
				h_u(iMicrolensArray)  = -projectionMatrix.intrinsicMatrix(1,3) * h_ui(iMicrolensArray);
				h_v(iMicrolensArray)  = -projectionMatrix.intrinsicMatrix(2,3) * h_vj(iMicrolensArray);
				h_uk(iMicrolensArray) = -projectionMatrix.incrementIntrinsicMatrix_k(1,3) * h_ui(iMicrolensArray);
				h_ul(iMicrolensArray) = -projectionMatrix.incrementIntrinsicMatrix_l(1,3) * h_ui(iMicrolensArray);
				h_vl(iMicrolensArray) = -projectionMatrix.incrementIntrinsicMatrix_l(2,3) * h_vj(iMicrolensArray);

				% Obtain information regarding positions (s,t)
				h_sk(iMicrolensArray) = -projectionMatrix.incrementExtrinsicMatrix_k(1,4);
				h_sl(iMicrolensArray) = -projectionMatrix.incrementExtrinsicMatrix_l(1,4);
				h_tl(iMicrolensArray) = -projectionMatrix.incrementExtrinsicMatrix_l(2,4);
				h_s(iMicrolensArray)  = -h_sk(iMicrolensArray) * cameraCoordinateOriginRay.k ...
									  -  h_sl(iMicrolensArray) * cameraCoordinateOriginRay.l;
				h_t(iMicrolensArray)  = -h_tl(iMicrolensArray) * cameraCoordinateOriginRay.l;

				% Obtain translation components
				translationVectors(:,iMicrolensArray) = [ projectionMatrix.extrinsicMatrix(1,4) + h_s(iMicrolensArray) ...
														; projectionMatrix.extrinsicMatrix(2,4) + h_t(iMicrolensArray) ...
														; projectionMatrix.extrinsicMatrix(3,4) ];
			end
			h_si = 0;
			h_tj = 0;
			h_sk = median(h_sk);
			h_sl = median(h_sl);
			h_s  = median(h_s);
			h_tl = median(h_tl);
			h_t  = median(h_t);
			h_ui = median(h_ui);
			h_uk = median(h_uk);
			h_ul = median(h_ul);
			h_u  = median(h_u);
			h_vj = median(h_vj);
			h_vl = median(h_vl);
			h_v  = median(h_v);
            
            % Obtain intrinsic matrix
            intrinsicMatrix  = [ h_si 0 h_sk h_sl h_s ...
                               ; 0 h_tj    0 h_tl h_t ...
                               ; h_ui 0 h_uk h_ul h_u ...
                               ; 0 h_vj    0 h_vl h_v ...
                               ; 0 0 0 0 1 ];

            plenopticCameras = []; 
            for iMicrolensArray = 1:numberMicrolensArrays
                % Obtain camera array projection matrix
                projectionMatrix  = microlensCameraArrays(iMicrolensArray);
                extrinsicMatrix   = camera.models.ExtrinsicMatrix.ExtrinsicMatrixFromMatrix(projectionMatrix.extrinsicMatrix);

                % Redefine translation vector
                extrinsicMatrix.translationVector = translationVectors(:,iMicrolensArray);
                
                % Create standard plenoptic camera instance
                posePlenopticCamera = camera.models.plenoptic.standard.Camera();
                posePlenopticCamera.intrinsicMatrix = intrinsicMatrix;
                posePlenopticCamera.extrinsic       = extrinsicMatrix;
                posePlenopticCamera.lightfieldSize  = lightfieldSize;
                
                % Update cameras
                plenopticCameras = cat(2,plenopticCameras,posePlenopticCamera);
            end
        end
        
        function plenopticCameras = CameraFromViewpointProjectionMatrix( viewpointCameraArrays ...
                                                                       , lightfieldSize ...
                                                                       , cameraCoordinateOriginRay )
            %
            % Create standard plenoptic camera model instance using
            % viewpoint camera array projection matrices.
            %
            % INPUTS:
            %   1. viewpointCameraArrays - viewpoint camera arrays 
            %   projection matrices.
            %   2. lightfieldSize   - lightfield size. 
            %   3. cameraCoordinateOriginRay - ray constraint for the
            %   origin of the camera coordinate system.
            %
            narginchk(2,3);
            
            % Transform to lightfield size
            lightfieldSize = utils.enums.Classes.LIGHTFIELD_SIZE().convert(lightfieldSize);

            if nargin <= 2
                % Obtain central ray
                cameraCoordinateOriginRay = lightfield.ImageRay([ lightfieldSize.numberPixels_i - 1 ...
                                                                , lightfieldSize.numberPixels_j - 1 ...
                                                                , lightfieldSize.numberMicrolenses_k - 1 ...
                                                                , lightfieldSize.numberMicrolenses_l - 1 ] ./ 2 + 1,false);
            end
            
            % Transform to image ray
            cameraCoordinateOriginRay = utils.enums.Classes.IMAGE_RAY().convert(cameraCoordinateOriginRay);

            % For each camera array projection matrix, create intrinsic 
            % matrix for plenoptic camera
            numberViewpointArrays = length(viewpointCameraArrays);
			h_si = zeros(1,numberViewpointArrays);
			h_s  = zeros(1,numberViewpointArrays);
			h_tj = zeros(1,numberViewpointArrays);
			h_t  = zeros(1,numberViewpointArrays);
			h_ui = zeros(1,numberViewpointArrays);
			h_uk = zeros(1,numberViewpointArrays);
			h_ul = zeros(1,numberViewpointArrays);
			h_u  = zeros(1,numberViewpointArrays);
			h_vj = zeros(1,numberViewpointArrays);
			h_vl = zeros(1,numberViewpointArrays);
			h_v  = zeros(1,numberViewpointArrays);
			translationVectors = zeros(3,numberViewpointArrays);
			for iViewpointArray = 1:numberViewpointArrays
				% Obtain camera array projection matrix
				projectionMatrix = viewpointCameraArrays(iViewpointArray);

				% Obtain information regarding directions (u,v)
				h_uk(iViewpointArray) = 1 / projectionMatrix.intrinsicMatrix(1,1);
				h_vl(iViewpointArray) = 1 / projectionMatrix.intrinsicMatrix(2,2);
				h_ul(iViewpointArray) = -projectionMatrix.intrinsicMatrix(1,2) * h_uk(iViewpointArray) * h_vl(iViewpointArray);
				h_v(iViewpointArray)  = -projectionMatrix.intrinsicMatrix(2,3) * h_vl(iViewpointArray);
				h_u(iViewpointArray)  = -projectionMatrix.intrinsicMatrix(1,3) * h_uk(iViewpointArray) ...
									  -  projectionMatrix.intrinsicMatrix(1,2) * h_v(iViewpointArray) * h_uk(iViewpointArray);
				h_ui(iViewpointArray) = -projectionMatrix.incrementIntrinsicMatrix_i(1,3) * h_uk(iViewpointArray);
				h_vj(iViewpointArray) = -projectionMatrix.incrementIntrinsicMatrix_j(2,3) * h_vl(iViewpointArray);

				% Obtain information regarding positions (s,t)
				h_si(iViewpointArray) = -projectionMatrix.incrementExtrinsicMatrix_i(1,4);
				h_tj(iViewpointArray) = -projectionMatrix.incrementExtrinsicMatrix_j(2,4);
				h_s(iViewpointArray)  = -h_si(iViewpointArray) * cameraCoordinateOriginRay.i;
				h_t(iViewpointArray)  = -h_tj(iViewpointArray) * cameraCoordinateOriginRay.j;

				% Obtain translation components
				translationVectors(:,iViewpointArray) = [ projectionMatrix.extrinsicMatrix(1,4) + h_s(iViewpointArray) ...
														; projectionMatrix.extrinsicMatrix(2,4) + h_t(iViewpointArray) ...
														; projectionMatrix.extrinsicMatrix(3,4) ];
			end
			h_sk = 0;
			h_sl = 0;
			h_tl = 0;
			h_si = median(h_si);
			h_tj = median(h_tj);
			h_s  = median(h_s);
			h_t  = median(h_t);
			h_ui = median(h_ui);
			h_uk = median(h_uk);
			h_ul = median(h_ul);
			h_u  = median(h_u);
			h_vl = median(h_vl);
			h_vj = median(h_vj);
			h_v  = median(h_v);                
            
            % Obtain intrinsic matrix
            intrinsicMatrix  = [ h_si 0 h_sk h_sl h_s ...
                               ; 0 h_tj    0 h_tl h_t ...
                               ; h_ui 0 h_uk h_ul h_u ...
                               ; 0 h_vj    0 h_vl h_v ...
                               ; 0 0 0 0 1 ];
            
            plenopticCameras = []; 
            for iViewpointArray = 1:numberViewpointArrays
                % Obtain camera array projection matrix
                projectionMatrix  = viewpointCameraArrays(iViewpointArray);
                
                % Obtain extrinsic parameters
                extrinsicMatrix = camera.models.ExtrinsicMatrix.ExtrinsicMatrixFromMatrix(projectionMatrix.extrinsicMatrix);
                extrinsicMatrix.translationVector = translationVectors(:,iViewpointArray);
                
                % Create standard plenoptic camera instance
                posePlenopticCamera = camera.models.plenoptic.standard.Camera();
                posePlenopticCamera.intrinsicMatrix = intrinsicMatrix;
                posePlenopticCamera.extrinsic       = extrinsicMatrix;
                posePlenopticCamera.lightfieldSize  = lightfieldSize;
                
                % Update cameras
                plenopticCameras = cat(2,plenopticCameras,posePlenopticCamera);
            end            
        end

        function plenopticCameras = CameraFromViewpointPinhole( viewpointCameraArrays ...
                                                              , lightfieldSize ...
                                                              , cameraCoordinateOriginRay )
            % 
            % Create standard plenoptic camera model instance using
            % viewpoint camera array pinhole cameras. This method
            % considers that the sampling can only be rectangular.
            %
            % INPUTS:
            %   1. viewpointCameraArrays - viewpoint pinhole camera arrays.
            %   2. lightfieldSize   - lightfield size. 
            %   3. cameraCoordinateOriginRay - ray constraint for the
            %   origin of the camera coordinate system.
            %
            narginchk(2,3);
            
            % Transform to lightfield size
            lightfieldSize = utils.enums.Classes.LIGHTFIELD_SIZE().convert(lightfieldSize);

            if nargin <= 2
                % Obtain central ray
                cameraCoordinateOriginRay = lightfield.ImageRay([ lightfieldSize.numberPixels_i - 1 ...
                                                                , lightfieldSize.numberPixels_j - 1 ...
                                                                , lightfieldSize.numberMicrolenses_k - 1 ...
                                                                , lightfieldSize.numberMicrolenses_l - 1 ] ./ 2 + 1,false);
            end
            
            % Obtain number of viewpoint camera arrays
            numberViewpointArrays = length(viewpointCameraArrays);
            
            % Obtain coordinates of projection matrix used to represent
            % viewpoint camera array
            centralCoordinates = lightfield.ImageRay([0;0;0;0],false);

            % Step 1 - obtain parameters related with directions (u,v)
            % For each viewpoint camera array, obtain observation matrix
            % and vector
            iLine = 0;
            observationMatrix = zeros(12 * numberViewpointArrays,6);
            observationVector = zeros(12 * numberViewpointArrays,1);
            for iViewpointArray = 1:numberViewpointArrays
                % Obtain projection matrix and increment matrix
                projectionMatrix = viewpointCameraArrays(iViewpointArray).projectionMatrix;
                incrementMatrix  = viewpointCameraArrays(iViewpointArray).incrementMatrix_i ...
                                 + viewpointCameraArrays(iViewpointArray).incrementMatrix_j;

                % Obtain rotation matrix 
                rotationMatrix   = projectionMatrix.extrinsicMatrix.rotationMatrix;
            
                % Projection matrix components
                observationMatrix(iLine + 1,1:3) = [ rotationMatrix.data(3,1) * centralCoordinates.i ...
                                                   , projectionMatrix.projectionMatrix(1,1) ...
                                                   , rotationMatrix.data(3,1) ];
                observationMatrix(iLine + 2,4:6) = [ rotationMatrix.data(3,1) * centralCoordinates.j ...
                                                   , projectionMatrix.projectionMatrix(2,1) ...
                                                   , rotationMatrix.data(3,1) ];
                observationMatrix(iLine + 3,1:3) = [ rotationMatrix.data(3,2) * centralCoordinates.i ...
                                                   , projectionMatrix.projectionMatrix(1,2) ...
                                                   , rotationMatrix.data(3,2) ];
                observationMatrix(iLine + 4,4:6) = [ rotationMatrix.data(3,2) * centralCoordinates.j ...
                                                   , projectionMatrix.projectionMatrix(2,2) ...
                                                   , rotationMatrix.data(3,2) ];
                observationMatrix(iLine + 5,1:3) = [ rotationMatrix.data(3,3) * centralCoordinates.i ...
                                                   , projectionMatrix.projectionMatrix(1,3) ...
                                                   , rotationMatrix.data(3,3) ];
                observationMatrix(iLine + 6,4:6) = [ rotationMatrix.data(3,3) * centralCoordinates.j ...
                                                   , projectionMatrix.projectionMatrix(2,3) ...
                                                   , rotationMatrix.data(3,3) ];
                observationVector(iLine + 1) = rotationMatrix.data(1,1);
                observationVector(iLine + 2) = rotationMatrix.data(2,1);
                observationVector(iLine + 3) = rotationMatrix.data(1,2);
                observationVector(iLine + 4) = rotationMatrix.data(2,2);
                observationVector(iLine + 5) = rotationMatrix.data(1,3);
                observationVector(iLine + 6) = rotationMatrix.data(2,3);

                % Increments matrix components
                observationMatrix(iLine + 07,1:2) = [rotationMatrix.data(3,1),incrementMatrix(1,1)];
                observationMatrix(iLine + 08,4:5) = [rotationMatrix.data(3,1),incrementMatrix(2,1)];
                observationMatrix(iLine + 09,1:2) = [rotationMatrix.data(3,2),incrementMatrix(1,2)];
                observationMatrix(iLine + 10,4:5) = [rotationMatrix.data(3,2),incrementMatrix(2,2)];
                observationMatrix(iLine + 11,1:2) = [rotationMatrix.data(3,3),incrementMatrix(1,3)];
                observationMatrix(iLine + 12,4:5) = [rotationMatrix.data(3,3),incrementMatrix(2,3)];
                
                % Update line
                iLine = iLine + 12;
            end
                
            % Decode parameters of intrinsic matrix H
            parameters_uv = mldivide(observationMatrix,observationVector);
            h_ui = parameters_uv(1);
            h_uk = parameters_uv(2);
            h_u  = parameters_uv(3);
            h_vj = parameters_uv(4);
            h_vl = parameters_uv(5);
            h_v  = parameters_uv(6);

            % Step 2 - obtain parameters related with positions (s,t) and
            % translation vector
            numberTranslationColumns = 3 * numberViewpointArrays;
			observationMatrix = zeros(6 * numberViewpointArrays + 3,numberTranslationColumns + 6);
			observationVector = zeros(6 * numberViewpointArrays + 3,1);
            
            % For each viewpoint camera array, obtain observation matrix
            % and vector
            iLine   = 0;
            iColumn = 0;
            for iViewpointArray = 1:numberViewpointArrays
                % Obtain projection matrix and increment matrix
                projectionMatrix = viewpointCameraArrays(iViewpointArray).projectionMatrix;
                incrementMatrix  = viewpointCameraArrays(iViewpointArray).incrementMatrix_i ...
                                 + viewpointCameraArrays(iViewpointArray).incrementMatrix_j;

                % Obtain projection center
                projectionCenter = camera.models.pinhole.Pinhole([projectionMatrix.intrinsicMatrix,zeros(3,1)]).projectionCenter;

                % Projection matrix components
                observationMatrix(iLine + 1,iColumn + 1) = -1;
                observationMatrix(iLine + 1,iColumn + 3) = h_u + h_ui * centralCoordinates.i;
                observationMatrix(iLine + 1,numberTranslationColumns + 1) = centralCoordinates.i;
                observationMatrix(iLine + 1,numberTranslationColumns + 3) = 1;
                observationMatrix(iLine + 2,iColumn + 2) = -1;
                observationMatrix(iLine + 2,iColumn + 3) = h_v + h_vj * centralCoordinates.j;
                observationMatrix(iLine + 2,numberTranslationColumns + 4) = centralCoordinates.j;
                observationMatrix(iLine + 2,numberTranslationColumns + 6) = 1;
                observationMatrix(iLine + 3,iColumn + 3)    = h_uk;
                observationMatrix(iLine + 3,numberTranslationColumns + 2) = 1;

                observationVector(iLine + 1) = -h_uk * projectionMatrix.projectionMatrix(1,4);
                observationVector(iLine + 2) = -h_vl * projectionMatrix.projectionMatrix(2,4);
                observationVector(iLine + 3) =  h_uk * projectionMatrix.projectionMatrix(3,4);
            
                % Increments matrix components
                observationMatrix(iLine + 4,iColumn + 3) = h_ui;
                observationMatrix(iLine + 4,numberTranslationColumns + 1) = 1;
                observationMatrix(iLine + 5,iColumn + 3) = h_vj;
                observationMatrix(iLine + 5,numberTranslationColumns + 4) = 1;

                observationVector(iLine + 4) = -h_uk * incrementMatrix(1,4);
                observationVector(iLine + 5) = -h_vl * incrementMatrix(2,4);
            
                % Plane (s,t) position containing viewpoint projection centers
				observationMatrix(iLine + 6,numberTranslationColumns + 2) = 1;
				observationVector(iLine + 6) = 0;
            
                % Update line and column
                iLine   = iLine + 6;
                iColumn = iColumn + 3;
            end
            
            % Central ray constraint
            observationMatrix(iLine + 1,numberTranslationColumns + (1:3)) = [cameraCoordinateOriginRay.i, cameraCoordinateOriginRay.k, 1];
            observationMatrix(iLine + 2,numberTranslationColumns + (4:6)) = [cameraCoordinateOriginRay.j, cameraCoordinateOriginRay.l, 1];

            % Viewpoint pinholes
			observationMatrix(iLine + 3,numberTranslationColumns + 2) =  h_vl;
			observationMatrix(iLine + 3,numberTranslationColumns + 5) = -h_uk;

            % Decode parameters translation vector and intrinsic matrix H
            parameters_st      = mldivide(observationMatrix,observationVector);
            translationVectors = parameters_st(1:numberTranslationColumns);
            translationVectors = reshape(translationVectors,3,numberViewpointArrays);
            h_si = parameters_st(numberTranslationColumns + 1);
            h_sk = parameters_st(numberTranslationColumns + 2);
            h_s  = parameters_st(numberTranslationColumns + 3);
            h_tj = parameters_st(numberTranslationColumns + 4);
            h_tl = parameters_st(numberTranslationColumns + 5);
            h_t  = parameters_st(numberTranslationColumns + 6);

            % Obtain intrinsic matrix H
            intrinsicMatrix = [ h_si 0 h_sk 0 h_s ...
                              ; 0 h_tj 0 h_tl h_t ...
                              ; h_ui 0 h_uk 0 h_u ...
                              ; 0 h_vj 0 h_vl h_v ...
                              ; 0 0 0 0 1 ];

            % There is also ambiguity on the solutions obtained for the
            % intrinsic matrix. Test if the solution is correct by ensuring
            % that the baseline distances are positive.
            intrinsic = camera.models.plenoptic.IntrinsicMatrixDirection(intrinsicMatrix);
            
            % Obtain extrinsic matrix for each viewpoint array
            plenopticCameras = [];
            for iViewpointArray = 1:numberViewpointArrays
                % Obtain projection matrix
                projectionMatrix = viewpointCameraArrays(iViewpointArray).projectionMatrix;
                
                % Obtain extrinsic matrix
                extrinsicMatrix  = projectionMatrix.extrinsicMatrix;
                extrinsicMatrix.translationVector = translationVectors(:,iViewpointArray);
            
                % Create standard plenoptic camera instance
                posePlenopticCamera = camera.models.plenoptic.standard.Camera();
                posePlenopticCamera.intrinsicMatrix = intrinsic.baselineCorrectedIntrinsicMatrix.data;
                posePlenopticCamera.extrinsic       = extrinsicMatrix;
                posePlenopticCamera.lightfieldSize  = lightfieldSize;
                
                % Update cameras and line
                plenopticCameras = cat(2,plenopticCameras,posePlenopticCamera);
            end
        end
    end
end
