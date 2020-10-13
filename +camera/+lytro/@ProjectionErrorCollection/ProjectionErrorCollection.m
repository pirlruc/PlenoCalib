classdef ProjectionErrorCollection
    %PROJECTIONERRORCOLLECTION
    %   Utility to evaluate projection error for lytro camera calibration.
    
    properties
        data = []           % Projection error instances collection
        plenoptic   = camera.models.plenoptic.standard.Camera()
                            % Plenoptic camera instance. This camera does
                            % not have information regarding the extrinsic
                            % matrix.
    end
    
    properties (Dependent)
        numberImages                % Number of plenoptic images used for reconstruction
        numberPoints                % Number of points considered for reconstruction
        groundTruthWorldPoints      % List all world ground truth points in collection
        groundTruthCameraPoints     % List all camera ground truth points in collection
        estimatedCameraPoints       % List all estimated points in collection
        groundTruthRays_ijkl        % Ground truth rays in image space
        projectionRays_ijkl         % Projection rays in image space
        projectionRays_stuv         % Projection rays in object space
        reconstructionError         % 3D reconstruction error
        projectionError             % Projection error on image plane. RMSE for each point
        pointToRayDistanceByPoint   % RMSE ray-reprojection error for each point
        pointToRayDistance          % RMSE ray-reprojection error for each ray
        dansereauPointToRayDistance % RMSE dansereau ray-reprojection error for all points
    end
    
    methods
        function self = ProjectionErrorCollection(varargin)
            %
            % Create projection error collection instance.
            %
            % INPUTS:
            %   1. projectionErrors - projection error instances to be
            %      executed.
            %   2. plenoptic  - plenoptic camera instance.
            %
            narginchk(0,2);

            if ~isempty(varargin)
                if nargin >= 2
                    self.plenoptic = varargin{2};
                end
                
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
        
        function self = set.plenoptic(self,newPlenoptic)
            self.plenoptic = utils.enums.Classes.STANDARD_PLENOPTIC_CAMERA().convert(newPlenoptic);
            
            % Since this is a collection of images with different extrinsic 
            % parameters, remove extrinsic parameters
            self.plenoptic.extrinsic = camera.models.ExtrinsicMatrix();
        end
        
        function number = get.numberImages(self)
            number = length(self.data);
        end
        
        function number = get.numberPoints(self)
            % If there is no projection error information, there are no
            % points
            if isempty(self.data)
                number = 0;
            else
                number = sum([self.data.numberPoints]);
            end
        end
        
        function data = get.groundTruthWorldPoints(self)
            % If there is no projection error information, provide empty 
            % list
            if isempty(self.data)
                data = math.Point();
            else
                points = [self.data.groundTruthWorldPoints];
                data   = math.Point([points.data],false);
            end
        end
        
        function data = get.groundTruthCameraPoints(self)
            % If there is no projection error information, provide empty 
            % list
            if isempty(self.data)
                data = math.Point();
            else
                points = [self.data.groundTruthCameraPoints];
                data   = math.Point([points.data],false);
            end
        end
        
        function data = get.estimatedCameraPoints(self)
            % If there is no projection error information, provide empty 
            % list
            if isempty(self.data)
                data = math.Point();
            else
                estimatedPoints = [self.data.estimatedCameraPoints];
                data            = math.Point([estimatedPoints.data],false);
            end
        end
        
        function data = get.groundTruthRays_ijkl(self)
            % If there is no projection error information, provide empty 
            % list
            if isempty(self.data)
                data = [];
            else
                data = [self.data.correspondences];
            end
        end
        
        function data = get.projectionRays_ijkl(self)
            % If there is no projection error information, provide empty 
            % list
            if isempty(self.data)
                data = [];
            else
                % Remove entries with no points
                indices = [self.data.numberPoints] > 0;

                % Concatenate projection rays
                data = [self.data(indices).projectionRays_ijkl];
            end
        end
        
        function data = get.projectionRays_stuv(self)
            % If there is no projection error information, provide empty 
            % list
            if isempty(self.data)
                data = [];
            else
                % Remove entries with no points
                indices = [self.data.numberPoints] > 0;

                % Concatenate projection rays
                data = [self.data(indices).projectionRays_stuv];
            end
        end
        
        function reconstruction = get.reconstructionError(self)
            % Obtain reconstruction error considering the l2-norm to
            % compute the distance between the estimated and ground truth
            % points in the camera coordinate system.
            %
            % This gives an error associated to each point.
            %
            reconstruction = math.Point( self.groundTruthCameraPoints.data - self.estimatedCameraPoints.data ...
                                       , false );
            reconstruction = reconstruction.norm.data;
        end
        
        function projection = get.projectionError(self)
            % For each ground truth point, project world point into image
            % sensor and compute the distance to the detected interest
            % points. The error for each point corresponds to the RMSE of
            % the projections.
            
            % Obtain camera type
            cameraType = self.data(1).cameraType;
            
            if self.numberPoints == 0
                projection = [];
            else
                % Obtain projection rays in image space
                estimatedRays_ijkl = self.projectionRays_ijkl;
                groundTruthRays    = self.groundTruthRays_ijkl;
                
                projection = nan(1,self.numberPoints);
                for iPoint = 1:self.numberPoints
                    % Obtain estimated and ground truth lightfield
                    % coordinates
                    pointEstimatedRays_ijkl   = estimatedRays_ijkl(iPoint);
                    pointGroundTruthRays_ijkl = groundTruthRays(iPoint);

                    % If number of vectors is equal compare points directly
                    matchIndividualRays = true;
                    if pointGroundTruthRays_ijkl.numberVectors == pointEstimatedRays_ijkl.numberVectors
                        if cameraType == camera.models.plenoptic.enums.CameraTypes.VIEWPOINT
                            indices_i = pointGroundTruthRays_ijkl.i == pointEstimatedRays_ijkl.i;
                            indices_j = pointGroundTruthRays_ijkl.j == pointEstimatedRays_ijkl.j;
                            if all(indices_i & indices_j) == true
                                projectionDistances = lightfield.ImageRay( pointGroundTruthRays_ijkl.data - pointEstimatedRays_ijkl.data ...
                                                                         , false );
                                matchIndividualRays = false;
                            end
                        else % Microlens projection
                            indices_k = pointGroundTruthRays_ijkl.k == pointEstimatedRays_ijkl.k;
                            indices_l = pointGroundTruthRays_ijkl.l == pointEstimatedRays_ijkl.l;
                            if all(indices_k & indices_l) == true
                                projectionDistances = lightfield.ImageRay( pointGroundTruthRays_ijkl.data - pointEstimatedRays_ijkl.data ...
                                                                         , false );
                                matchIndividualRays = false;
                            end
                        end
                    end
                            
                    % If number of vectors is different, match each image
                    % ray coordinates and obtain distance. Assume that
                    % features are on viewpoint images
                    if matchIndividualRays == true
                        projectionDistances = nan(4,pointGroundTruthRays_ijkl.numberVectors);
                        for iRay = 1:pointGroundTruthRays_ijkl.numberVectors
                            % Obtain image ray
                            imageRay = pointGroundTruthRays_ijkl.obtainVectors(iRay);
                            
                            % Match ground truth ray
                            if cameraType == camera.models.plenoptic.enums.CameraTypes.VIEWPOINT
                                indices_i     = imageRay.i == pointEstimatedRays_ijkl.i;
                                indices_j     = imageRay.j == pointEstimatedRays_ijkl.j;
                                cameraIndices = indices_i & indices_j;
                            else
                                indices_k     = imageRay.k == pointEstimatedRays_ijkl.k;
                                indices_l     = imageRay.l == pointEstimatedRays_ijkl.l;
                                cameraIndices = indices_k & indices_l;
                            end
                            
                            % If there is a match compute projection
                            % distance
                            if any(cameraIndices)
                                % Obtain ground truth ray
                                groundTruthRay = pointEstimatedRays_ijkl.obtainVectors(cameraIndices);
                                
                                % Compute projection distance. Due to
                                % distortion, there may be more than one
                                % ray per viewpoint. Therefore, consider a
                                % mean error
                                projectionDistances(:,iRay) = mean(imageRay.data - groundTruthRay.data,2);
                            end
                        end
                        projectionDistances = lightfield.ImageRay(projectionDistances,false);
                    end
                    
                    % Remove nan entries
                    indices             = ~isnan(projectionDistances.i);
                    projectionDistances = projectionDistances.obtainVectors(indices);
                    
                    % Obtain mean projection distance for point using 
                    % l2-norm
                    projectionDistances = projectionDistances.norm;
                    projectionDistances = sum(projectionDistances.data) ...
                                        / projectionDistances.numberVectors;
                    projection(iPoint)  = projectionDistances;
                end
            end
        end
        
        function rayReprojection = get.pointToRayDistanceByPoint(self)
            % For each ground truth point, transform point from the world
            % to camera coordinate system and project the point to the ray.
            % This measures the distance from 3D points to 3D rays.
            %
            % This gives an error associated to each point.
            %
            
            if self.numberPoints == 0
                rayReprojection = [];
            else
                % Obtain camera ground truth points
                groundTruthPoints  = self.groundTruthCameraPoints;

                % Obtain rays in object space
                estimatedRays_stuv = self.projectionRays_stuv;
  
                rayReprojection = nan(1,self.numberPoints);
                for iPoint = 1:self.numberPoints
                    % Obtain ground truth point and estimated rays in
                    % object space
                    cameraGroundTruthPoint  = groundTruthPoints.obtainVectors(iPoint);
                    pointEstimatedRays_stuv = estimatedRays_stuv(iPoint);
                    
                    % Obtain (s,t)-plane point and (u,v)-plane point. The
                    % (u,v)-plane point corresponds to a direction since 
                    % the plane is at distance equal to 1 and (u,v) are 
                    % local coordinates.
                    pointObjectRays_st = math.Point( [ pointEstimatedRays_stuv.s;pointEstimatedRays_stuv.t ...
                                                     ; zeros(1,pointEstimatedRays_stuv.numberVectors) ], false );
                    pointObjectRays_uv = math.Point( [ pointEstimatedRays_stuv.u;pointEstimatedRays_stuv.v ...
                                                     ; ones(1,pointEstimatedRays_stuv.numberVectors) ], false );

                    % Normalize ray directions
                    pointObjectRays_uv = pointObjectRays_uv.data ./ pointObjectRays_uv.norm.data;

                    % Obtain direction defined between (s,t)-plane point 
                    % and camera point and project to camera rays
                    pointRayDirection  = repmat(cameraGroundTruthPoint.data,1,pointObjectRays_st.numberVectors) ...
                                       - pointObjectRays_st.data;
                    pointRayProjection = repmat(dot(pointRayDirection,pointObjectRays_uv),3,1) ...
                                      .* pointObjectRays_uv;

                    % Obtain distance between the camera point ray and the
                    % projected camera ray
                    distanceVector = math.Vector(pointRayDirection - pointRayProjection,false);

                    % Store l2-norm for distance vector
                    rayReprojection(iPoint) = sum(distanceVector.norm.data) / distanceVector.numberVectors;
                end
            end
        end
        
        function rayReprojection = get.pointToRayDistance(self)
            % For each ground truth point, transform point from the world
            % to camera coordinate system and project the point to the ray.
            % This measures the distance from 3D points to 3D rays.
            %
            % This gives an error associated to each ray.
            %
            
            if self.numberPoints == 0
                rayReprojection = [];
            else
                % Obtain camera ground truth points
                groundTruthPoints  = self.groundTruthCameraPoints;

                % Obtain rays in object space
                estimatedRays_stuv = self.projectionRays_stuv;
  
                rayReprojection = [];
                for iPoint = 1:self.numberPoints
                    % Obtain ground truth point and estimated rays in
                    % object space
                    cameraGroundTruthPoint  = groundTruthPoints.obtainVectors(iPoint);
                    pointEstimatedRays_stuv = estimatedRays_stuv(iPoint);
                    
                    % Obtain (s,t)-plane point and (u,v)-plane point. The
                    % (u,v)-plane point corresponds to a direction since 
                    % the plane is at distance equal to 1 and (u,v) are 
                    % local coordinates.
                    pointObjectRays_st = math.Point( [ pointEstimatedRays_stuv.s;pointEstimatedRays_stuv.t ...
                                                     ; zeros(1,pointEstimatedRays_stuv.numberVectors) ], false );
                    pointObjectRays_uv = math.Point( [ pointEstimatedRays_stuv.u;pointEstimatedRays_stuv.v ...
                                                     ; ones(1,pointEstimatedRays_stuv.numberVectors) ], false );

                    % Normalize ray directions
                    pointObjectRays_uv = pointObjectRays_uv.data ./ pointObjectRays_uv.norm.data;

                    % Obtain direction defined between (s,t)-plane point 
                    % and camera point and project to camera rays
                    pointRayDirection  = repmat(cameraGroundTruthPoint.data,1,pointObjectRays_st.numberVectors) ...
                                       - pointObjectRays_st.data;
                    pointRayProjection = repmat(dot(pointRayDirection,pointObjectRays_uv),3,1) ...
                                      .* pointObjectRays_uv;

                    % Obtain distance between the camera point ray and the
                    % projected camera ray
                    distanceVector = math.Vector(pointRayDirection - pointRayProjection,false);

                    % Store l2-norm for distance vector
                    rayReprojection = cat(2,rayReprojection,distanceVector.norm.data);
                end
            end
        end
        
        function rayReprojection = get.dansereauPointToRayDistance(self)
            % For each ground truth point, transform point from the world
            % to camera coordinate system and project the point to the ray.
            % This measures the distance from 3D points to 3D rays.
            %
            % This gives an unique error that corresponds to the mean value
            % of all points.
            %
            
            if self.numberPoints == 0
                rayReprojection = [];
            else
                % Obtain camera ground truth points
                groundTruthPoints  = self.groundTruthCameraPoints;

                % Obtain rays in object space
                estimatedRays_stuv = self.projectionRays_stuv;

                rayReprojection = [];
                for iPoint = 1:self.numberPoints
                    % Obtain ground truth point and estimated rays in
                    % object space
                    cameraGroundTruthPoint  = groundTruthPoints.obtainVectors(iPoint);
                    pointEstimatedRays_stuv = estimatedRays_stuv(iPoint);

                    % Obtain (s,t)-plane point and (u,v)-plane point. The
                    % (u,v)-plane point corresponds to a direction since 
                    % the plane is at distance equal to 1 and (u,v) are 
                    % local coordinates.
                    pointObjectRays_st = math.Point( [ pointEstimatedRays_stuv.s;pointEstimatedRays_stuv.t ...
                                                     ; zeros(1,pointEstimatedRays_stuv.numberVectors) ], false );
                    pointObjectRays_uv = math.Point( [ pointEstimatedRays_stuv.u;pointEstimatedRays_stuv.v ...
                                                     ; ones(1,pointEstimatedRays_stuv.numberVectors) ], false );

                    % Normalize ray directions
                    pointObjectRays_uv = pointObjectRays_uv.data ./ pointObjectRays_uv.norm.data;

                    % Obtain direction defined between (s,t)-plane point 
                    % and camera point and project to camera rays
                    pointRayDirection  = repmat(cameraGroundTruthPoint.data,1,pointObjectRays_st.numberVectors) ...
                                       - pointObjectRays_st.data;
                    pointRayProjection = repmat(dot(pointRayDirection,pointObjectRays_uv),3,1) ...
                                      .* pointObjectRays_uv;

                    % Obtain distance between the camera point ray and the
                    % projected camera ray
                    distanceVector = math.Vector(pointRayDirection - pointRayProjection,false);

                    % Store l2-norm for distance vector
                    rayReprojection = cat(2,rayReprojection,distanceVector.norm.data);
                end
                
                % Obtain RMSE for all points
                rayReprojection = sqrt(mean(rayReprojection.^2));
            end
        end
    end
    
    methods (Static)
		function self = ProjectionErrorCollectionFromCalibrationFile( calibrationFilepath,reconstructionMethod )
            %
            % Obtain projection error collection from calibration file.
            %
            % INPUTS:
            %   1. calibrationFilepath  - filepath to calibration file.
            %   2. reconstructionMethod - reconstruction method.
            %
            narginchk(1,2);
            
            % If reconstruction method is not provided, assume point
            % reconstruction.
            if nargin <= 1
                reconstructionMethod = camera.models.plenoptic.enums.ReconstructionTypes.POINT_RECONSTRUCTION();
            end
            
            % Obtain calibration folder path
            calibrationFolderPath = fileparts(calibrationFilepath);
            
            % Obtain calibration parameters
            [intrinsic,extrinsic,distortion,lightfieldSize] = ...
                    camera.lytro.Camera.obtainCalibrationParametersFromCalibrationFile( calibrationFilepath );

            % Create plenoptic camera instance with intrinsic matrix and
            % lightfield size.
            plenopticCamera = camera.models.plenoptic.standard.Camera( intrinsic ...
                                                                     , lightfieldSize ...
                                                                     , distortion );
            
            % Obtain corner points information 
            cornersFilepath = [calibrationFolderPath filesep 'CheckerboardCorners.mat'];
            [imageCorners, worldCornersGroundTruth] = camera.lytro.Camera.obtainCornersFromCornersFile(cornersFilepath);
            
            % Obtain general information
            numberImages = size(imageCorners,3);
            
            projections     = [];
            posePlenopticCamera = plenopticCamera;
            for iImage = 1:numberImages
                % Obtain extrinsic parameters
                rotationVector      = extrinsic.rotationVector(iImage,:);
                translationVector   = extrinsic.translationVector(iImage,:);
                extrinsicParameters = [translationVector,rotationVector];
                
                % Update plenoptic camera information with extrinsic matrix
                % parameters
                posePlenopticCamera.extrinsic = camera.models.ExtrinsicMatrix.ExtrinsicMatrixFromEncodedParameters(extrinsicParameters);
                
                % Obtain correspondences from calibration pattern points in
                % viewpoint images. Remember that each pose is given in the
                % last dimension
                numberPoints    = worldCornersGroundTruth.numberVectors;
                correspondences = camera.lytro.Camera.obtainCorrespondencesFromCorners( ...
                                            imageCorners(:,:,iImage) ...
                                          , numberPoints );
                
                % Create reconstruction instance and add it to the
                % collection.
                projections = cat( 1, projections ...
                                    , camera.lytro.ProjectionError( posePlenopticCamera ...
                                                                  , correspondences ...
                                                                  , worldCornersGroundTruth ...
                                                                  , reconstructionMethod ...
                                                                  , camera.models.plenoptic.enums.CameraTypes.VIEWPOINT ));
            end
            
            % Create point reconstruction collection instance
            self = camera.lytro.ProjectionErrorCollection( projections ...
                                                         , plenopticCamera );
        end
    end
end
