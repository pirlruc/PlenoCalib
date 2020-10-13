classdef ProjectionError < camera.lytro.PointReconstruction
    %PROJECTIONERROR
    %   Utility to evaluate projection error for lytro camera calibration.
    
    properties
        cameraType = camera.models.plenoptic.enums.CameraTypes.VIEWPOINT()    % Camera type considered for projections
    end
    
    properties (Dependent)
        projectionRays_ijkl         % Projection rays in image space
        projectionRays_stuv         % Projection rays in object space
        reconstructionError         % 3D reconstruction error
        projectionError             % Projection error on image plane. RMSE for each point
        pointToRayDistanceByPoint   % RMSE ray-reprojection error for each point
        pointToRayDistance          % RMSE ray-reprojection error for each ray
        dansereauPointToRayDistance % RMSE dansereau ray-reprojection error for all points
    end
    
    methods
        function self = ProjectionError(varargin)
            %
            % Create projection error collection instance. In this object, 
            % each item in correspondences corresponds to a different 
            % point.
            %
            % INPUTS:
            %   1. plenoptic   - plenoptic camera instance.
            %   2. correspondences   - correspondences identified in 
            %   ligthfield.
            %   3. groundTruthWorldPoints - ground truth points in the
            %   world coordinate system.
            %   4. reconstructionMethod   - reconstruction method.
            %   5. cameraType - camera type to compute projections.
            %
            narginchk(0,5);

            % Create super class instance
            self = self@camera.lytro.PointReconstruction();
            
            if ~isempty(varargin)
                if nargin >= 5
                    self.cameraType = varargin{5};
                end
                
                if nargin >= 4
                    self.reconstructionMethod = varargin{4};
                else
                    self.reconstructionMethod = camera.models.plenoptic.enums.ReconstructionTypes.LINE_PARAMETER_RECONSTRUCTION();
                end
                
                if nargin >= 3
                    self.groundTruthWorldPoints = varargin{3};
                end
                
                if nargin >= 2
                    self.correspondences = varargin{2};
                end
                
                if nargin >= 1
                    self.plenoptic = varargin{1};
                end
            end
        end
        
        function self = set.cameraType(self,newCameraType)
            self.cameraType = newCameraType;
        end
        
        function rays_ijkl = get.projectionRays_ijkl(self)
            % For each ground truth point, project world point into image
            % sensor.
            if self.cameraType == camera.models.plenoptic.enums.CameraTypes.VIEWPOINT
                PROJECTION_TYPE                = camera.lytro.enums.ProjectionTypes.VIEWPOINTS();
                PIXELS_TO_NEAREST_INTEGER      = true;
                MICROLENSES_TO_NEAREST_INTEGER = false;
                ADD_DISTORTION                 = true;
            else
                PROJECTION_TYPE                = camera.lytro.enums.ProjectionTypes.MICROLENSES();
                PIXELS_TO_NEAREST_INTEGER      = false;
                MICROLENSES_TO_NEAREST_INTEGER = true;
                ADD_DISTORTION                 = true;
            end
            
            if self.numberPoints == 0
                rays_ijkl = lightfield.ImageRay();
            else
                rays_ijkl = [];
                for iPoint = 1:self.numberPoints
                    % Obtain world ground truth point
                    groundTruthPoint = self.groundTruthWorldPoints.obtainVectors(iPoint);

                    % Project ground truth point
                    pointRays_ijkl   = self.plenoptic.project( groundTruthPoint ...
                                                             , PIXELS_TO_NEAREST_INTEGER ...
                                                             , MICROLENSES_TO_NEAREST_INTEGER ...
                                                             , ADD_DISTORTION ...
                                                             , PROJECTION_TYPE );
                                                          
                    % Add rays in image space
                    rays_ijkl = cat(1,rays_ijkl,pointRays_ijkl);
                end
            end
        end
        
        function rays_stuv = get.projectionRays_stuv(self)
            % For each ground truth point, transform point from the world
            % to camera coordinate system and project the point to the ray.
            REMOVE_DISTORTION = true;
            
            if self.numberPoints == 0
                rays_stuv = lightfield.ObjectRay();
            else
                rays_stuv = [];
                for iPoint = 1:self.numberPoints
                    % Obtain image ray coordinates for correspondences and 
                    % transform to (s,t,u,v) coordinates.
                    pointImageRays_ijkl = self.correspondences(iPoint);
                    pointRays_stuv      = self.plenoptic.obtainObjectRaysFromImageRays( pointImageRays_ijkl ...
                                                                                      , REMOVE_DISTORTION );

                    % Add rays in image space
                    rays_stuv = cat(2,rays_stuv,pointRays_stuv);
                end
            end
        end
        
        function reconstruction = get.reconstructionError(self)
            % Obtain reconstruction error considering the l2-norm to
            % compute the distance between the estimated and ground truth
            % points in the camera coordinate system.
            reconstruction = math.Point( self.groundTruthCameraPoints.data - self.estimatedCameraPoints.data ...
                                       , false );
            reconstruction = reconstruction.norm.data;
        end
        
        function projection = get.projectionError(self)
            % For each ground truth point, project world point into image
            % sensor and compute the distance to the detected interest
            % points. The error for each point corresponds to the RMSE of
            % the projections.
            
            if self.numberPoints == 0
                projection = [];
            else
                % Obtain projection rays in image space
                estimatedRays_ijkl = self.projectionRays_ijkl;
                
                projection = nan(1,self.numberPoints);
                for iPoint = 1:self.numberPoints
                    % Obtain estimated and ground truth lightfield
                    % coordinates
                    pointEstimatedRays_ijkl   = estimatedRays_ijkl(iPoint);
                    pointGroundTruthRays_ijkl = self.correspondences(iPoint);

                    % If number of vectors is equal compare points directly
                    matchIndividualRays = true;
                    if pointGroundTruthRays_ijkl.numberVectors == pointEstimatedRays_ijkl.numberVectors
                        % Sanity check, just to ensure that coordinates are the same   
                        if self.cameraType == camera.models.plenoptic.enums.CameraTypes.VIEWPOINT
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
                            if self.cameraType == camera.models.plenoptic.enums.CameraTypes.VIEWPOINT
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
            % This measures the distance from 3D points to 3D rays. Gives
            % the RMSE for each point.
            
            if self.numberPoints == 0
                rayReprojection = [];
            else
                % Obtain camera ground truth points
                groundTruthCameraPoints = self.groundTruthCameraPoints;

                % Obtain rays in object space
                estimatedRays_stuv      = self.projectionRays_stuv;
                
                rayReprojection = nan(1,self.numberPoints);
                for iPoint = 1:self.numberPoints
                    % Obtain ground truth point and estimated rays in
                    % object space
                    cameraGroundTruthPoint  = groundTruthCameraPoints.obtainVectors(iPoint);
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
            % This measures the distance from 3D points to 3D rays. Gives
            % the RMSE for each ray.
            
            if self.numberPoints == 0
                rayReprojection = [];
            else
                % Obtain camera ground truth points
                groundTruthCameraPoints = self.groundTruthCameraPoints;

                % Obtain rays in object space
                estimatedRays_stuv      = self.projectionRays_stuv;
                
                rayReprojection = [];
                for iPoint = 1:self.numberPoints
                    % Obtain ground truth point and estimated rays in
                    % object space
                    cameraGroundTruthPoint  = groundTruthCameraPoints.obtainVectors(iPoint);
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
            
            if self.numberPoints == 0
                rayReprojection = [];
            else
                % Obtain camera ground truth points
                groundTruthCameraPoints = self.groundTruthCameraPoints;

                % Obtain rays in object space
                estimatedRays_stuv      = self.projectionRays_stuv;

                rayReprojection = [];
                for iPoint = 1:self.numberPoints
                    % Obtain ground truth point and estimated rays in
                    % object space
                    cameraGroundTruthPoint  = groundTruthCameraPoints.obtainVectors(iPoint);
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
end
