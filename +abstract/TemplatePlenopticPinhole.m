classdef TemplatePlenopticPinhole
    %TEMPLATEPLENOPTICPINHOLE
    %   Template plenoptic pinhole camera instance
    
    properties
        projectionMatrix  = camera.models.pinhole.Pinhole()     % Pinhole projection matrix
        incrementMatrix_x = zeros(3,4)                          % Increment matrix data for x-coordinate
        incrementMatrix_y = zeros(3,4)                          % Increment matrix data for y-coordinate
    end
    
    properties (Dependent)
        baseline            % Baseline between consecutive (0,0) and (1,1) camera
    end
    
    methods
        function self = TemplatePlenopticPinhole(varargin)
            %
            % Create template plenoptic pinhole instance.
            %
            % INPUTS:
            %   1. projectionMatrix  - projection matrix data
            %   2. incrementMatrix_x - increment matrix data for
            %   x-coordinate.
            %   3. incrementMatrix_y - increment matrix data for
            %   y-coordinate.
            %
            narginchk(0,3);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.projectionMatrix = varargin{1};
                end
                
                if nargin >= 2
                    self.incrementMatrix_x = varargin{2};
                end
                
                if nargin >= 3
                    self.incrementMatrix_y = varargin{3};
                end
            end
        end
        
        function self = set.projectionMatrix(self,newProjectionMatrix)
            self.projectionMatrix = utils.enums.Classes.PINHOLE().convert(newProjectionMatrix);
        end
        
        function baseline = get.baseline(self)
            % Obtain projection center for (0,0)
            projectionCenter_0 = self.projectionMatrix.projectionCenter;
            
            % Obtain projection matrix for (1,1)
            CAMERA_INDICES     = [1;1];
            projectionMatrix_1 = self.obtainProjectionMatrices(CAMERA_INDICES);
            projectionCenter_1 = projectionMatrix_1.projectionCenter;
            
            % Obtain baseline
            baseline = math.Point( projectionCenter_1.data ...
                                 - projectionCenter_0.data, false );
        end
    end
    
    methods
        function projectionMatrices = obtainProjectionMatrices(self,cameraIndices)
            %
            % Obtain projection matrices for cameras.
            %
            % INPUTS:
            %   1. cameraIndices - camera indices to obtain projection 
            %      matrices. Each camera should be provided in different
            %      columns.
            %
            narginchk(1,2);
            
            if nargin <= 1
                [camera_x,camera_y] = meshgrid(1:4:10,1:4:10);
                cameraIndices = [camera_x(:),camera_y(:)]';
            end
            cameraIndices = utils.enums.Classes.PIXEL().convert(cameraIndices);
            
            % Define projection matrices
            projectionMatrices = [];
            for iCamera = 1:cameraIndices.numberVectors
                % Obtain projection matrix
                projectionMatrixData = self.projectionMatrix.projectionMatrix ...
                                     + cameraIndices.u(iCamera) .* self.incrementMatrix_x ...
                                     + cameraIndices.v(iCamera) .* self.incrementMatrix_y;

                % Add it to the collection
                projectionMatrices = cat( 1, projectionMatrices ...
                                           , camera.models.pinhole.Pinhole( projectionMatrixData ));
            end
        end
        
        function rays_ijkl = project(self,worldPoints,cameraIndices)
            %
            % Obtain projection of world points in image sensor.
            %
            % INPUTS: 
            %   1. worldPoints - points defined in the world coordinate
            %   system. The points should not be given in homogeneous 
            %   coordinates.
            %   2. cameraIndices - camera indices to obtain the
            %   corresponding projections.
            % 
            narginchk(3,3);
            
            % Transform world points and camera indices to vector class 
            worldPoints   = utils.enums.Classes.POINT().convert(worldPoints);
            cameraIndices = utils.enums.Classes.PIXEL().convert(cameraIndices);
            
            % Define observation matrix and projection matrix
            observationMatrix_x  = arrayfun(@(x) x .* eye(3,3),cameraIndices.u,'UniformOutput',false);
            observationMatrix_y  = arrayfun(@(x) x .* eye(3,3),cameraIndices.v,'UniformOutput',false);
            observationMatrix_1  = repmat(eye(3,3),cameraIndices.numberVectors,1);
            observationMatrix_x  = [observationMatrix_x{:}]';
            observationMatrix_y  = [observationMatrix_y{:}]';
            observationMatrix    = [observationMatrix_1,observationMatrix_x,observationMatrix_y];
            projection           = [ self.projectionMatrix.projectionMatrix ...
                                   ; self.incrementMatrix_x ...
                                   ; self.incrementMatrix_y ];
            
            % Project all points in each camera given
            worldPoints = worldPoints.setHomogeneousCoordinates();
            points_uv   = observationMatrix * projection * worldPoints.data;

            % Order pixels (u,v)
            points_u = points_uv(1:3:end,:);
            points_v = points_uv(2:3:end,:);
            points_1 = points_uv(3:3:end,:);
            
            % Segregate pixels by points
            rays_ijkl = [];
            for iPoint = 1:worldPoints.numberVectors
                % Obtain pixels (u,v) removing homogeneous coordinates
                pixels_uv = image.Pixel( [ points_u(:,iPoint)' ...
                                         ; points_v(:,iPoint)' ...
                                         ; points_1(:,iPoint)' ], true );
                pixels_uv = pixels_uv.removeHomogeneousCoordinates();

                % Obtain coordinates (i,j,k,l)
                pointRays_ijkl = lightfield.ImageRay( [ cameraIndices.u; cameraIndices.v ...
                                                      ; pixels_uv.u; pixels_uv.v ], false );
                
                rays_ijkl = cat(1,rays_ijkl,pointRays_ijkl);
            end
        end
        
        function worldPoint = reconstruct(self,rays_ijkl,cameraType)
            %
            % Reconstruct world point from pixels in image sensor.
            %
            % INPUTS: 
            %   1. rays_ijkl - image ray coordinates for a given point.
            %   2. cameraType - define plenoptic camera type. Since the
            %   implementation of the homography is made based on the
            %   rectangular structure of the viewpoints. The default is
            %   viewpoint type.
            % 
            narginchk(2,3);
            
            if nargin <= 2
                cameraType = camera.models.plenoptic.enums.CameraTypes.VIEWPOINT();
            end
            
            % Transform world points and camera indices to vector class 
            rays_ijkl = utils.enums.Classes.IMAGE_RAY().convert(rays_ijkl);
            
            % If one is considering the microlens camera type switch
            % coordinates of the image rays.
            if cameraType == camera.models.plenoptic.enums.CameraTypes.MICROLENS
                rays_ijkl.data = rays_ijkl.data(cameraType.encoding,:);
            end
            
            % Initialize matrices
            cameraMatrix      = zeros(2 * rays_ijkl.numberVectors,3 * 2 * rays_ijkl.numberVectors);
            observationMatrix = zeros(2 * 3 * rays_ijkl.numberVectors,3);
            observationVector = zeros(2 * 3 * rays_ijkl.numberVectors,1);
            
            iRow    = 0;
            iColumn = 0;
            for iCamera = 1:rays_ijkl.numberVectors
                % Define camera matrix
                cameraMatrix(iRow + 1,iColumn + 1:iColumn + 3) = [1,rays_ijkl.i(iCamera),rays_ijkl.j(iCamera)];
                cameraMatrix(iRow + 2,iColumn + 4:iColumn + 6) = [1,rays_ijkl.i(iCamera),rays_ijkl.j(iCamera)];
                
                % Obtain projection matrices
                projection = cat(3,self.projectionMatrix.projectionMatrix,self.incrementMatrix_x);
                projection = cat(3,projection,self.incrementMatrix_y);
                
                % Define observation matrix and observation vector
                % x-coordinate
                observationMatrix(iRow + 1:iRow + 3,1) = projection(3,1,:) .* rays_ijkl.k(iCamera) - projection(1,1,:);
                observationMatrix(iRow + 1:iRow + 3,2) = projection(3,2,:) .* rays_ijkl.k(iCamera) - projection(1,2,:);
                observationMatrix(iRow + 1:iRow + 3,3) = projection(3,3,:) .* rays_ijkl.k(iCamera) - projection(1,3,:);
                observationVector(iRow + 1:iRow + 3,1) = projection(3,4,:) .* rays_ijkl.k(iCamera) - projection(1,4,:);
                
                % y-coordinate
                observationMatrix(iRow + 4:iRow + 6,1) = projection(3,1,:) .* rays_ijkl.l(iCamera) - projection(2,1,:);
                observationMatrix(iRow + 4:iRow + 6,2) = projection(3,2,:) .* rays_ijkl.l(iCamera) - projection(2,2,:);
                observationMatrix(iRow + 4:iRow + 6,3) = projection(3,3,:) .* rays_ijkl.l(iCamera) - projection(2,3,:);
                observationVector(iRow + 4:iRow + 6,1) = projection(3,4,:) .* rays_ijkl.l(iCamera) - projection(2,4,:);
                
                iRow    = iRow    + 6;
                iColumn = iColumn + 6;
            end            
            
            % Obtain estimate for world point
            worldPoint = mldivide(cameraMatrix * observationMatrix,-cameraMatrix * observationVector);
            worldPoint = math.Point(worldPoint,false);
        end
    end
    
    methods (Static)
        function self = PinholeFromProjectionMatrix(projectionMatrix)
            % 
            % Obtain plenoptic pinhole from plenoptic projection matrix.
            %
            % INPUTS: 
            %   1. projectionMatrix - plenoptic projection matrix instance.
            %
            narginchk(1,1);
            
            % Obtain coordinates of camera to get incremental entries for
            % pinhole.
            cameraIndices = image.Pixel([1,0;0,1],false);
            
            % Obtain projection matrices
            projectionMatrices = projectionMatrix.obtainProjectionMatrices(cameraIndices);
            
            % Create pinhole instance
            self = abstract.TemplatePlenopticPinhole();
            self.projectionMatrix  = projectionMatrix.projectionMatrix;
            self.incrementMatrix_x = projectionMatrices(1).projectionMatrix - projectionMatrix.projectionMatrix;
            self.incrementMatrix_y = projectionMatrices(2).projectionMatrix - projectionMatrix.projectionMatrix;
        end
        
        function [self,rmse,sse,residuals] = PinholeFromPointCorrespondences( worldPoints, imageRays ...
                                                                            , observationMatrixMethod ...
                                                                            , cameraType )
            %
            % Obtain parameters to describe a camera array using 
            % direct linear transform from point correspondences.
            %
            % INPUTS:
            %   1. worldPoints - points defined in the world coordinate
            %   system. Each point should be defined in a different column.
            %   The points should not be given in homogeneous coordinates.
            %   2. imageRays - projected image rays in the image plane
            %   corresponding to the world points given as input. The image
            %   rays should be given as a cell.
            %   3. observationMatrixMethod - method to construct
            %   observation matrix. Default is the cross product method.
            %   4. cameraType - define plenoptic camera type. The default 
            %   is viewpoint type.
            %
            narginchk(2,4);
            
            if nargin <= 2
                observationMatrixMethod = camera.models.pinhole.enums.ObservationMatrixMethods.CROSS_PRODUCT();
            end
            
            if nargin <= 3
                cameraType = camera.models.plenoptic.enums.CameraTypes.VIEWPOINT;
            end
            
            if nargout >= 2
				[self,rmse,sse,residuals] = abstract.TemplatePlenopticPinhole.PinholeFromPointCorrespondences( ...
								worldPoints, imageRays ...
							  , observationMatrixMethod ...
							  , cameraType );
            else
				self = abstract.TemplatePlenopticPinhole.PinholeFromPointCorrespondences( ...
								worldPoints, imageRays ...
							  , observationMatrixMethod ...
							  , cameraType );
            end
        end
        
        function [self,rmse,sse,residuals] = PinholeFromPointCorrespondences( worldPoints, imageRays ...
                                                                            , observationMatrixMethod ...
                                                                             cameraType )
            %
            % Obtain parameters to describe a camera array using 
            % direct linear transform from point correspondences.
            %
            % INPUTS:
            %   1. worldPoints - points defined in the world coordinate
            %   system. Each point should be defined in a different column.
            %   The points should not be given in homogeneous coordinates.
            %   2. imageRays - projected image rays in the image plane
            %   corresponding to the world points given as input. The image
            %   rays should be given as a cell.
            %   3. observationMatrixMethod - method to construct
            %   observation matrix. Default is the cross product method.
            %   4. cameraType - define plenoptic camera type. Since the
            %   implementation of the viewpoint is made based on the
            %   rectangular structure of the viewpoints. The default is
            %   viewpoint type.
            %
            narginchk(2,4);
            
            if nargin <= 2
                observationMatrixMethod = camera.models.pinhole.enums.ObservationMatrixMethods.CROSS_PRODUCT();
            end
            
            if nargin <= 3
                cameraType = camera.models.plenoptic.enums.CameraTypes.VIEWPOINT;
            end
            
            % Number of unknown parameters
            PARAMETERS_TO_ESTIMATE = 20;

            % Transform world points and image rays to vector class 
            worldPoints = utils.enums.Classes.POINT().convert(worldPoints);
            imageRays   = utils.enums.Classes.IMAGE_RAY().convert(imageRays);

            % If one is considering the microlens camera type switch
            % coordinates of the image rays.
            if cameraType == camera.models.plenoptic.enums.CameraTypes.MICROLENS
                imageRays.data = imageRays.data(cameraType.encoding,:);
            end
            
            % Normalize lightfield coordinates and world points
			% Obtain (i,j) coordinates and treat them has pixels
			pixels_ij           = image.Pixel( [ [imageRays.i] ...
											   ; [imageRays.j] ],false );

			% Obtain (k,l) coordinates and treat them has pixels
			microlenses_kl = image.Pixel( [ [imageRays.k] ...
										  ; [imageRays.l] ],false );
			
			% Obtain normalization matrices
			raysNormalizationMatrix_ij = pixels_ij.obtainNormalizationMatrix();
			raysNormalizationMatrix_kl = microlenses_kl.obtainNormalizationMatrix();
			pointsNormalizationMatrix  = worldPoints.obtainNormalizationMatrix();

			% Obtain normalized points
			worldPoints = worldPoints.normalize;

			% Obtain ray normalization matrix
			raysNormalizationMatrix = eye(5,5);
			raysNormalizationMatrix(1:2,1:2) = raysNormalizationMatrix_ij(1:2,1:2);
			raysNormalizationMatrix(1:2,5)   = raysNormalizationMatrix_ij(1:2,3);
			raysNormalizationMatrix(3:4,3:4) = raysNormalizationMatrix_kl(1:2,1:2);
			raysNormalizationMatrix(3:4,5)   = raysNormalizationMatrix_kl(1:2,3);
            
            % Set world points and image rays to homogeneous coordinates
            worldPoints = worldPoints.setHomogeneousCoordinates();
            imageRays   = imageRays.setHomogeneousCoordinates();
                
            % Normalize image rays
            imageRays.data = raysNormalizationMatrix * imageRays.data;

            %
            % Remember that in this implementation (i,j) are the camera
            % indices and (k,l) are the pixels.
            %
            % Construct direct linear transformation calibration matrix
            % considering (k,l) as the pixel coordinates on the
            % viewpoint images, (i,j) as the pixel coordinates on the
            % microlens images, and M as the world point coordinates:
            %
            %   - Linear Solution
            %       A = [ M^T 0   -k*M^T i*M^T 0 
            %           ; 0   M^T -l*M^T 0     j*M^T ]
            if observationMatrixMethod == camera.models.pinhole.enums.ObservationMatrixMethods.LINEAR_SYSTEM
                % Initialize direct linear transformation matrix
                % Each tuple image ray and point gives two equations
                dltMatrix = zeros(2 * imageRays.numberVectors,PARAMETERS_TO_ESTIMATE);

                % Obtain point coordinates for each projection
                worldPointsData = worldPoints.data';

                % Point data
                dltMatrix(1:imageRays.numberVectors,01:04)       = worldPointsData;
                dltMatrix(imageRays.numberVectors + 1:end,05:08) = worldPointsData;

                % Projection data
                dltMatrix(1:imageRays.numberVectors,09:12)       = -imageRays.k' .* worldPointsData;
                dltMatrix(imageRays.numberVectors + 1:end,09:12) = -imageRays.l' .* worldPointsData;
                dltMatrix(1:imageRays.numberVectors,13:16)       =  imageRays.i' .* worldPointsData;
                dltMatrix(imageRays.numberVectors + 1:end,17:20) =  imageRays.j' .* worldPointsData;

            %   - Cross-Product
            %       A = [ M^T     0     -k*M^T  i*M^T     0 
            %           ; 0      -M^T    l*M^T  0        -j*M^T 
            %           ; -l*M^T  k*M^T  0     -l*i*M^T   k*j*M^T ]
            elseif observationMatrixMethod == camera.models.pinhole.enums.ObservationMatrixMethods.CROSS_PRODUCT
                % Initialize direct linear transformation matrix
                % Each tuple image ray and point gives three equations
                dltMatrix = zeros(3 * imageRays.numberVectors,PARAMETERS_TO_ESTIMATE);

                % Obtain point coordinates for each projection
                worldPointsData = worldPoints.data';

                % First row
                % Point data
                dltMatrix(1:imageRays.numberVectors,01:04) =  worldPointsData;
                % Projection data
                dltMatrix(1:imageRays.numberVectors,09:12) = -imageRays.k' .* worldPointsData;
                dltMatrix(1:imageRays.numberVectors,13:16) =  imageRays.i' .* worldPointsData;

                % Second row
                % Point data
                dltMatrix( imageRays.numberVectors + 1 ...
                         : 2 * imageRays.numberVectors,05:08) = -worldPointsData;
                % Projection data
                dltMatrix( imageRays.numberVectors + 1 ...
                         : 2 * imageRays.numberVectors,09:12) =  imageRays.l' .* worldPointsData;
                dltMatrix( imageRays.numberVectors + 1 ...
                         : 2 * imageRays.numberVectors,17:20) = -imageRays.j' .* worldPointsData;

                % Third row
                % Projection data
                dltMatrix( 2 * imageRays.numberVectors + 1:end,01:04) = -imageRays.l' .* worldPointsData;
                dltMatrix( 2 * imageRays.numberVectors + 1:end,05:08) =  imageRays.k' .* worldPointsData;
                dltMatrix( 2 * imageRays.numberVectors + 1:end,13:16) = -imageRays.l' .* imageRays.i' .* worldPointsData;
                dltMatrix( 2 * imageRays.numberVectors + 1:end,17:20) =  imageRays.k' .* imageRays.j' .* worldPointsData;
            end

            % Since the number of observations is huge, perform 
            % QR-decomposition to obtain an orthogonal matrix (Q) and an 
            % upper triangular matrix (R)
            [~,upperTriangular] = qr(dltMatrix,0);

            % Detect if matrix is ill-conditioned
            if rank(upperTriangular) < min(size(dltMatrix)) - 1
                error = MException( 'TemplatePlenopticPinhole:illConditionedRectangular' ...
                                  , [ 'The observation matrix for the homography being estimated ' ...
                                      'is ill-conditioned. The solution is not correct.' ] );
                error.throw();
            end
            
            % The solution corresponds to the null space of the upper
            % triangular matrix
            [~,~,singularVectors] = svd(upperTriangular);
            nullSpaceSolution     = singularVectors(:, end);
            
            % Obtain residuals and errors for least squares problem
            if nargout > 1
                residuals = abs(dltMatrix * nullSpaceSolution);
                sse  = sum(residuals.^2);               % Sum of squared errors
                rmse = sqrt(sse / length(residuals));   % Root mean squared error
            end
            
            % Obtain viewpoint projection matrix. The projection matrix has
            % size (3 x 4)
            projectionMatrixData = reshape(nullSpaceSolution(01:12) , 4, 3)';

            % Obtain matrix with same size of projection matrix to obtain
            % projection matrices from other viewpoints.
            incrementMatrixData = [reshape(nullSpaceSolution(13:end), 4, 2)';zeros(1,4)];
            
            % This problem has two solutions (x and -x). Choose the
            % solution that has the homogeneous coordinate positive.
            homogeneousCoordinate = projectionMatrixData(3,:) * worldPoints.data;
            if any(homogeneousCoordinate < 0)
                projectionMatrixData = -projectionMatrixData;
                incrementMatrixData  = -incrementMatrixData;
            end
            
            % Obtain unnormalized matrices
            projectionMatrixData = mldivide(raysNormalizationMatrix_kl,projectionMatrixData * pointsNormalizationMatrix);
            incrementMatrixData  = mldivide(raysNormalizationMatrix_kl,incrementMatrixData  * pointsNormalizationMatrix);

            % Since projection matrix is defined up to a scale factor
            % normalize projection matrix with ||[p31 p32 p33]|| = 1.
            projectionMatrixScaleFactor = norm(projectionMatrixData(3,1:3));
            projectionMatrixData        = projectionMatrixData ./ projectionMatrixScaleFactor ;
            incrementMatrixData         = incrementMatrixData  ./ projectionMatrixScaleFactor;
            
            % Now let us correct the matrix due to the normalization of the
            % pixels (i,j) of the microlens
            projectionMatrixData = projectionMatrixData ...
                                 + diag(raysNormalizationMatrix_ij(:,3)) * incrementMatrixData;
            incrementMatrixData  = diag( [ raysNormalizationMatrix_ij(1,1) ...
                                         , raysNormalizationMatrix_ij(2,2) ...
                                         , 1 ] ) * incrementMatrixData;

            % Let us isolate the contribution of each coordinate (i,j)
            incrementMatrixData_i      = zeros(3,4);
            incrementMatrixData_i(1,:) = incrementMatrixData(1,:);
            incrementMatrixData_j      = zeros(3,4);
            incrementMatrixData_j(2,:) = incrementMatrixData(2,:);
            
            % Create template pinhole instance
            self = abstract.TemplatePlenopticPinhole();
            self.projectionMatrix  = camera.models.pinhole.Pinhole(projectionMatrixData);
            self.incrementMatrix_x = incrementMatrixData_i;
            self.incrementMatrix_y = incrementMatrixData_j;
        end
    end
end
