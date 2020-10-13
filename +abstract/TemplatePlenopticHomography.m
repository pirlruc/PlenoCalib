classdef TemplatePlenopticHomography
    %TEMPLATEPLENOPTICHOMOGRAPHY
    %   Template plenoptic homography matrix instance.
    
    properties
        homographyMatrix  = eye(3,3)         % Homography matrix
        incrementMatrix_x = zeros(3,3)       % Increment matrix data for x-coordinate
        incrementMatrix_y = zeros(3,3)       % Increment matrix data for y-coordinate
    end
    
    methods
        function self = TemplatePlenopticHomography(varargin)
            %
            % Create template plenoptic homography instance.
            %
            % INPUTS:
            %   1. homographyMatrix  - homography matrix data
            %   2. incrementMatrix_x - increment matrix data for
            %   x-coordinate.
            %   3. incrementMatrix_y - increment matrix data for
            %   y-coordinate.
            %
            narginchk(0,3);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.homographyMatrix = varargin{1};
                end
                
                if nargin >= 2
                    self.incrementMatrix_x = varargin{2};
                end
                
                if nargin >= 3
                    self.incrementMatrix_y = varargin{3};
                end
            end
        end
    end
    
    methods
        function homographyMatrices = obtainHomographyMatrices(self,cameraIndices)
            %
            % Obtain homography matrices for cameras.
            %
            % INPUTS:
            %   1. cameraIndices - camera indices to obtain homography
            %      matrices. Each camera should be provided in different
            %      columns.
            %
            narginchk(1,2);
            
            if nargin <= 1
                [camera_x,camera_y] = meshgrid(1:4:10,1:4:10);
                cameraIndices       = [camera_x(:),camera_y(:)]';
            end
            cameraIndices = utils.enums.Classes.PIXEL().convert(cameraIndices);
            
            % Define homography matrices
            homographyMatrices = [];
            for iCamera = 1:cameraIndices.numberVectors
                % Obtain homography matrix
                homographyMatrixData = self.homographyMatrix ...
                                     + cameraIndices.u(iCamera) .* self.incrementMatrix_x ...
                                     + cameraIndices.v(iCamera) .* self.incrementMatrix_y;

                % Add it to the collection
                homographyMatrices = cat( 1, homographyMatrices ...
                                           , camera.models.pinhole.Homography( homographyMatrixData ));
            end
        end
        
        function homography = obtainHomographyMatrixForTargetCamera( self, sourceCameraIndex ...
                                                                   , targetCameraIndex, targetHomographyArray )
            %
            % Obtain homography matrix between cameras source and target
            % cameras.
            %
            % INPUTS:
            %   1. sourceCameraIndex - coordinates of the source camera.
            %   Default is camera (5,5).
            %   2. targetCameraIndex - coordinates of the target camera.
            %   Default is (+1,+1) camera relatively to the source camera.
            %   3. targetHomographyArray - target camera array homography
            %   object. This assumes that the current class is the source
            %   homography camera array. If not provided, the homography is
            %   computed considering a unique pose.
            %
            narginchk(1,4);
            
            if nargin <= 1
                sourceCameraIndex = [5;5];
            end
            sourceCameraIndex = utils.enums.Classes.PIXEL().convert(sourceCameraIndex);
            
            if nargin <= 2
                targetCameraIndex = sourceCameraIndex.data + 1;
            end
            targetCameraIndex = utils.enums.Classes.PIXEL().convert(targetCameraIndex);
            
            if nargin <= 3
                targetHomographyArray = self;
            end
            
            % Obtain homography matrix of each camera
            sourceHomographyMatrix = self.obtainHomographyMatrices(sourceCameraIndex);
            targetHomographyMatrix = targetHomographyArray.obtainHomographyMatrices(targetCameraIndex);
            
            % Obtain homography matrix
            homography = mldivide( sourceHomographyMatrix.homographyMatrix' ...
                                 , targetHomographyMatrix.homographyMatrix' )';

            % Normalize homography matrix
            homography = homography ./ homography(3,3);
            homography = camera.models.pinhole.Homography(homography);
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
            homographyProjection = [ self.homographyMatrix ...
                                   ; self.incrementMatrix_x ...
                                   ; self.incrementMatrix_y ];
            
            % Project all points in each camera given
            worldPixels = image.Pixel([worldPoints.x;worldPoints.y],false);
            worldPixels = worldPixels.setHomogeneousCoordinates();
            points_uv   = observationMatrix * homographyProjection * worldPixels.data;

            % Order pixels (u,v)
            points_u = points_uv(1:3:end,:);
            points_v = points_uv(2:3:end,:);
            points_1 = points_uv(3:3:end,:);
            
            % Segregate pixels by points
            rays_ijkl = [];
            for iPoint = 1:worldPixels.numberVectors
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
            % Reconstruct world point from pixels in image sensor. This
            % will give the (x,y)-coordinates of the world point.
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
            observationMatrix = zeros(2 * 3 * rays_ijkl.numberVectors,2);
            observationVector = zeros(2 * 3 * rays_ijkl.numberVectors,1);
            
            iRow    = 0;
            iColumn = 0;
            for iCamera = 1:rays_ijkl.numberVectors
                % Define camera matrix
                cameraMatrix(iRow + 1,iColumn + 1:iColumn + 3) = [1,rays_ijkl.i(iCamera),rays_ijkl.j(iCamera)];
                cameraMatrix(iRow + 2,iColumn + 4:iColumn + 6) = [1,rays_ijkl.i(iCamera),rays_ijkl.j(iCamera)];
                
                % Obtain homography matrices
                homography = cat(3,self.homographyMatrix,self.incrementMatrix_x);
                homography = cat(3,homography,self.incrementMatrix_y);
                
                % Define observation matrix and observation vector
                % x-coordinate
                observationMatrix(iRow + 1:iRow + 3,1) = homography(3,1,:) .* rays_ijkl.k(iCamera) - homography(1,1,:);
                observationMatrix(iRow + 1:iRow + 3,2) = homography(3,2,:) .* rays_ijkl.k(iCamera) - homography(1,2,:);
                observationVector(iRow + 1:iRow + 3,1) = homography(3,3,:) .* rays_ijkl.k(iCamera) - homography(1,3,:);
                
                % y-coordinate
                observationMatrix(iRow + 4:iRow + 6,1) = homography(3,1,:) .* rays_ijkl.l(iCamera) - homography(2,1,:);
                observationMatrix(iRow + 4:iRow + 6,2) = homography(3,2,:) .* rays_ijkl.l(iCamera) - homography(2,2,:);
                observationVector(iRow + 4:iRow + 6,1) = homography(3,3,:) .* rays_ijkl.l(iCamera) - homography(2,3,:);
                
                iRow    = iRow    + 6;
                iColumn = iColumn + 6;
            end            
            
            % Obtain estimate for world point
            worldPoint = mldivide(cameraMatrix * observationMatrix,-cameraMatrix * observationVector);
            worldPoint = math.Point([worldPoint;0],false);
        end
    end

    methods (Static)
        function [self,rmse,sse,residuals] = HomographyFromPointCorrespondences( worldPoints, imageRays ...
                                                                               , observationMatrixMethod ...
                                                                               , cameraType )
            %
            % Obtain parameters to describe a camera array homography using
            % direct linear transform from point correspondences.
            %
            % There should be a point for each image ray even if the point
            % is the same.
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
            
            % Number of unknown parameters
            PARAMETERS_TO_ESTIMATE = 15;

            % Transform world points and image rays to vector class 
            worldPoints = utils.enums.Classes.POINT().convert(worldPoints);
            imageRays   = utils.enums.Classes.IMAGE_RAY().convert(imageRays);

            % If one is considering the microlens camera type switch
            % coordinates of the image rays.
            if cameraType == camera.models.plenoptic.enums.CameraTypes.MICROLENS
                imageRays.data = imageRays.data(cameraType.encoding,:);
            end
            
            % Since for the world points we are going to consider only
            % the (x,y) coordinates, these can be considered as pixels
            worldPixels = image.Pixel([worldPoints.x;worldPoints.y],false);

            % Normalize lightfield coordinates and world points
			% Obtain (i,j) coordinates and treat them has pixels
			pixels_ij      = image.Pixel( [ [imageRays.i] ...
										  ; [imageRays.j] ],false );

			% Obtain (k,l) coordinates and treat them has pixels
			microlenses_kl = image.Pixel( [ [imageRays.k] ...
										  ; [imageRays.l] ],false );
			
			% Obtain normalization matrices
			raysNormalizationMatrix_ij = pixels_ij.obtainNormalizationMatrix();
			raysNormalizationMatrix_kl = microlenses_kl.obtainNormalizationMatrix();
			pointsNormalizationMatrix  = worldPixels.obtainNormalizationMatrix();

			% Obtain normalized points
			worldPixels = worldPixels.normalize;

			% Obtain ray normalization matrix
			raysNormalizationMatrix = eye(5,5);
			raysNormalizationMatrix(1:2,1:2) = raysNormalizationMatrix_ij(1:2,1:2);
			raysNormalizationMatrix(1:2,5)   = raysNormalizationMatrix_ij(1:2,3);
			raysNormalizationMatrix(3:4,3:4) = raysNormalizationMatrix_kl(1:2,1:2);
			raysNormalizationMatrix(3:4,5)   = raysNormalizationMatrix_kl(1:2,3);
            
            % Set world points and image rays to homogeneous coordinates
            worldPixels = worldPixels.setHomogeneousCoordinates();
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
                worldPointsData = worldPixels.data';

                % Point data
                dltMatrix(1:imageRays.numberVectors,01:03)       = worldPointsData;
                dltMatrix(imageRays.numberVectors + 1:end,04:06) = worldPointsData;

                % Projection data
                dltMatrix(1:imageRays.numberVectors,07:09)       = -imageRays.k' .* worldPointsData;
                dltMatrix(imageRays.numberVectors + 1:end,07:09) = -imageRays.l' .* worldPointsData;
                dltMatrix(1:imageRays.numberVectors,10:12)       =  imageRays.i' .* worldPointsData;
                dltMatrix(imageRays.numberVectors + 1:end,13:15) =  imageRays.j' .* worldPointsData;

            %   - Cross-Product
            %       A = [ M^T     0     -k*M^T  i*M^T    0 
            %           ; 0      -M^T    l*M^T  0       -j*M^T 
            %           ; -l*M^T  k*M^T  0     -l*i*M^T  k*j*M^T ]
            elseif observationMatrixMethod == camera.models.pinhole.enums.ObservationMatrixMethods.CROSS_PRODUCT
                % Initialize direct linear transformation matrix
                % Each tuple image ray and point gives three equations
                dltMatrix = zeros(3 * imageRays.numberVectors,PARAMETERS_TO_ESTIMATE);

                % Obtain point coordinates for each projection
                worldPointsData = worldPixels.data';

                % First row
                % Point data
                dltMatrix(1:imageRays.numberVectors,01:03) =  worldPointsData;
                % Projection data
                dltMatrix(1:imageRays.numberVectors,07:09) = -imageRays.k' .* worldPointsData;
                dltMatrix(1:imageRays.numberVectors,10:12) =  imageRays.i' .* worldPointsData;

                % Second row
                % Point data
                dltMatrix( imageRays.numberVectors + 1 ...
                         : 2 * imageRays.numberVectors,04:06) = -worldPointsData;
                % Projection data
                dltMatrix( imageRays.numberVectors + 1 ...
                         : 2 * imageRays.numberVectors,07:09) =  imageRays.l' .* worldPointsData;
                dltMatrix( imageRays.numberVectors + 1 ...
                         : 2 * imageRays.numberVectors,13:15) = -imageRays.j' .* worldPointsData;

                % Third row
                % Projection data
                dltMatrix( 2 * imageRays.numberVectors + 1:end,01:03) = -imageRays.l' .* worldPointsData;
                dltMatrix( 2 * imageRays.numberVectors + 1:end,04:06) =  imageRays.k' .* worldPointsData;
                dltMatrix( 2 * imageRays.numberVectors + 1:end,10:12) = -imageRays.l' .* imageRays.i' .* worldPointsData;
                dltMatrix( 2 * imageRays.numberVectors + 1:end,13:15) =  imageRays.k' .* imageRays.j' .* worldPointsData;
            end

            % Since the number of observations is huge, perform
            % QR-decomposition to obtain an orthogonal matrix (Q) and an 
            % upper triangular matrix (R)
            [~,upperTriangular] = qr(dltMatrix,0);
            
            % Detect if matrix is ill-conditioned
            if rank(upperTriangular) < min(size(dltMatrix)) - 1
                error = MException( 'ViewpointHomography:illConditioned' ...
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
            
            % Obtain homography matrix. The homography matrix has size 
            % (3 x 3)
            homographyMatrixData = reshape(nullSpaceSolution(01:09) , 3, 3)';

            % Obtain matrix with same size of homography matrix to obtain
            % homography matrices from other viewpoints.
			incrementMatrixData = [reshape(nullSpaceSolution(10:end), 3, 2)';zeros(1,3)];

            % Obtain unnormalized matrices
            homographyMatrixData = mldivide(raysNormalizationMatrix_kl,homographyMatrixData * pointsNormalizationMatrix);
            incrementMatrixData  = mldivide(raysNormalizationMatrix_kl,incrementMatrixData  * pointsNormalizationMatrix);

			% Ensure the homography has the same signal
			homographyMatrixScaleFactor = sign(homographyMatrixData(3,3));
            homographyMatrixData = homographyMatrixData ./ homographyMatrixScaleFactor;
            incrementMatrixData  = incrementMatrixData  ./ homographyMatrixScaleFactor;

            % Now let us correct the matrix due to the normalization of the
            % pixels (i,j) of the microlens
            homographyMatrixData = homographyMatrixData ...
                                 + diag(raysNormalizationMatrix_ij(:,3)) * incrementMatrixData;
            incrementMatrixData  = diag( [ raysNormalizationMatrix_ij(1,1) ...
                                         , raysNormalizationMatrix_ij(2,2) ...
                                         , 1 ] ) * incrementMatrixData;
            
            % Let us isolate the contribution of each coordinate (i,j)
            incrementMatrixData_i      = zeros(3,3);
            incrementMatrixData_i(1,:) = incrementMatrixData(1,:);
            incrementMatrixData_j      = zeros(3,3);
            incrementMatrixData_j(2,:) = incrementMatrixData(2,:);

            % Create template pinhole instance
            self = abstract.TemplatePlenopticHomography();
            self.homographyMatrix  = homographyMatrixData;
            self.incrementMatrix_x = incrementMatrixData_i;
            self.incrementMatrix_y = incrementMatrixData_j;
        end
    end
end
