classdef TemplatePlenopticProjectionMatrix
    %TEMPLATEPLENOPTICPROJECTIONMATRIX
    %   Template plenoptic projection matrix camera instance
    
    properties (Constant)
        POSE_CODE_SCALE      = 100  % Scale for creating pose code
        LENS_CODE_SCALE      = 10   % Scale for creating lens code
        DIMENSION_CODE_SCALE = 1    % Scale for creating dimension code
    end
    
    properties
        intrinsicMatrix = camera.models.pinhole.IntrinsicMatrix().intrinsicMatrix
                % Intrinsic matrix that allows to project points in the camera coordinate system
        extrinsicMatrix = camera.models.ExtrinsicMatrix().extrinsicMatrix
                % Extrinsic matrix that allows to convert points from the world to the camera 
                % coordinate system
        incrementIntrinsicMatrix_x  = zeros(3,3)   % Increment intrinsic matrix data for x-coordinate
        incrementIntrinsicMatrix_y  = zeros(3,3)   % Increment intrinsic matrix data for y-coordinate
        incrementExtrinsicMatrix_x  = zeros(3,4)   % Increment extrinsic matrix data for x-coordinate
        incrementExtrinsicMatrix_y  = zeros(3,4)   % Increment extrinsic matrix data for y-coordinate
    end
    
    properties (Dependent)
        projectionMatrix    % Projection matrix that allows to project a point in the world 
                            % coordinate system to the image plane
        baseline            % Baseline between consecutive (0,0) and (1,1) camera
    end
    
    methods
        function self = TemplatePlenopticProjectionMatrix(varargin)
            %
            % Create microlens projection matrix instance.
            %
            % INPUTS:
            %   1. intrinsicMatrix - intrinsic matrix data.
            %   2. extrinsicMatrix - extrinsic matrix data.
            %   3. incrementIntrinsicMatrix_x - increment intrinsic matrix 
            %   data for x-coordinate.
            %   4. incrementIntrinsicMatrix_y - increment intrinsic matrix 
            %   data for y-coordinate.
            %   5. incrementExtrinsicMatrix_x - increment extrinsic matrix 
            %   data for x-coordinate.
            %   6. incrementExtrinsicMatrix_y - increment extrinsic matrix 
            %   data for y-coordinate.
            %
            narginchk(0,6);
            
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
        
        function projection = get.projectionMatrix(self)
            projection = self.intrinsicMatrix * self.extrinsicMatrix;
            % Since projection matrix is defined up to a scale factor
            % normalize projection matrix with ||[p31 p32 p33]|| = 1.
            projection = projection ./ norm(projection(3,1:3));
        end
        
        function baseline = get.baseline(self)
            % The extrinsic matrix difference corresponds to the baseline.
            baseline = math.Point( self.incrementExtrinsicMatrix_x(:,end) ...
                                 + self.incrementExtrinsicMatrix_y(:,end), false );
        end
    end
    
    methods
        function [intrinsicParameters,extrinsicParameters] = encode(self)
            %
            % Encode intrinsic and extrinsic parameters of camera array.
            %
            
            % Obtain intrinsic parameters
            intrinsicParameters  = [ self.intrinsicMatrix(:) ...
                                   ; self.incrementIntrinsicMatrix_x(:) ...
                                   ; self.incrementIntrinsicMatrix_y(:) ...
                                   ; self.incrementExtrinsicMatrix_x(:) ...
                                   ; self.incrementExtrinsicMatrix_y(:) ]';

            % Obtain extrinsic parameters
            extrinsics = camera.models.ExtrinsicMatrix.ExtrinsicMatrixFromMatrix(self.extrinsicMatrix);
            extrinsicParameters = extrinsics.encode();
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

            % Convert projection matrix to pinhole
            pinhole = abstract.TemplatePlenopticPinhole.PinholeFromProjectionMatrix(self);

            % Obtain projections
            rays_ijkl = pinhole.project(worldPoints,cameraIndices);
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
            
            % Convert projection matrix to pinhole
            pinhole = abstract.TemplatePlenopticPinhole.PinholeFromProjectionMatrix(self);

            % Reconstruct point
            worldPoint = pinhole.reconstruct(rays_ijkl,cameraType);
        end
        
        function projectionMatrices = obtainProjectionMatrices(self,cameraIndices)
            %
            % Obtain projection matrices for cameras.
            %
            % INPUTS:
            %   1. cameraIndices - camera indices to obtain projection
            %      matrices. Each microlens should be provided in different
            %      columns.
            %
            narginchk(1,2);
            
            if nargin <= 1
                [camera_x,camera_y] = meshgrid(1:4:10,1:4:10);
                cameraIndices       = [camera_x(:),camera_y(:)]';
            end
            cameraIndices = utils.enums.Classes.PIXEL().convert(cameraIndices);
            
            % Define projection matrices
            projectionMatrices = [];
            for iCamera = 1:cameraIndices.numberVectors
                % Obtain intrinsic and extrinsic matrix
                intrinsicMatrixData  = self.intrinsicMatrix ...
                                     + cameraIndices.u(iCamera) .* self.incrementIntrinsicMatrix_x ...
                                     + cameraIndices.v(iCamera) .* self.incrementIntrinsicMatrix_y;
                extrinsicMatrixData  = self.extrinsicMatrix ...
                                     + cameraIndices.u(iCamera) * self.incrementExtrinsicMatrix_x ...
                                     + cameraIndices.v(iCamera) * self.incrementExtrinsicMatrix_y;

                % Obtain projection matrix
                cameraProjectionMatrix = camera.models.pinhole.ProjectionMatrix();
                cameraProjectionMatrix.intrinsicMatrix = intrinsicMatrixData;
                cameraProjectionMatrix.extrinsicMatrix = extrinsicMatrixData;
                                 
                % Add it to the collection
                projectionMatrices = cat( 1, projectionMatrices ...
                                           , cameraProjectionMatrix );
            end
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
            
            INTRINSIC_MATRIX_PARAMETERS = 9;
            EXTRINSIC_MATRIX_PARAMETERS = 12;
            
            % Create projection matrix instance
            self = abstract.TemplatePlenopticProjectionMatrix();
            
            % Decode intrinsic parameters
            iParameter = 0;
            for iIntrinsic = 1:3
                intrinsicParameters = parameters(iParameter + 1:iParameter + INTRINSIC_MATRIX_PARAMETERS);
                if iIntrinsic == 1
                    self.intrinsicMatrix = reshape(intrinsicParameters,3,3);
                elseif iIntrinsic == 2
                    self.incrementIntrinsicMatrix_x = reshape(intrinsicParameters,3,3);
                else
                    self.incrementIntrinsicMatrix_y = reshape(intrinsicParameters,3,3);
                end
                iParameter = iParameter + INTRINSIC_MATRIX_PARAMETERS;
            end
            
            % Decode extrinsic parameters
            for iExtrinsic = 1:2
                extrinsicParameters = parameters(iParameter + 1:iParameter + EXTRINSIC_MATRIX_PARAMETERS);
                if iExtrinsic == 1
                    self.incrementExtrinsicMatrix_x = reshape(extrinsicParameters,3,4);
                else
                    self.incrementExtrinsicMatrix_y = reshape(extrinsicParameters,3,4);
                end
                iParameter = iParameter + EXTRINSIC_MATRIX_PARAMETERS;
            end
        end
        
        function [intrinsicVariables,numberParameters,constraintVariables,dependencyVariables] = ...
                                            intrisicParametersToOptimize( cameraType ...
                                                                        , pinholes )
            %
            % Obtain intrinsic parameters to be optimized on camera array.
            %
            % INPUTS:
            %   1. cameraType - camera type for calibration.
            %   2. pinholes   - pinholes to be ensured during calibration.
            %
            narginchk(2,2);
            
            INTRINSIC_MATRIX = [1 0 1;0 1 1;0 0 0];
            INCREMENT_INTRINSIC_MATRIX_X = [0 0 1;0 0 0;0 0 0];
            INCREMENT_INTRINSIC_MATRIX_Y = [0 0 0;0 0 1;0 0 0];
            INCREMENT_EXTRINSIC_MATRIX_X = [0 0 0 1;0 0 0 0;0 0 0 0];
            INCREMENT_EXTRINSIC_MATRIX_Y = [0 0 0 0;0 0 0 1;0 0 0 0];
            PINHOLE_CONSTRAINT           = [0 0 0 0;0 0 0 0;0 0 0 0];
            DEPENDENCIES_INTRINSIC       = [1 0 0;0 1 0;0 0 0];
            
			incrementExtrinsicMatrix_x = INCREMENT_EXTRINSIC_MATRIX_X;
			incrementExtrinsicMatrix_y = INCREMENT_EXTRINSIC_MATRIX_Y;
            
            % Obtain parameters to optimize with and without constraint
            default    = [ INTRINSIC_MATRIX(:) ...
                         ; INCREMENT_INTRINSIC_MATRIX_X(:) ...
                         ; INCREMENT_INTRINSIC_MATRIX_Y(:) ...
                         ; INCREMENT_EXTRINSIC_MATRIX_X(:) ...
                         ; INCREMENT_EXTRINSIC_MATRIX_Y(:) ]';
            parameters = [ INTRINSIC_MATRIX(:) ...
                         ; INCREMENT_INTRINSIC_MATRIX_X(:) ...
                         ; INCREMENT_INTRINSIC_MATRIX_Y(:) ...
                         ; incrementExtrinsicMatrix_x(:) ...
                         ; incrementExtrinsicMatrix_y(:) ]';
            intrinsicVariables  = find(parameters);
            constraintVariables = find(default - parameters);
            numberParameters    = length(default);
            
            % Obtain dependencies for constraints
            if ~isempty(constraintVariables)
                dependencyVariables = [ DEPENDENCIES_INTRINSIC(:) ...
                                      ; INCREMENT_INTRINSIC_MATRIX_X(:) ...
                                      ; INCREMENT_INTRINSIC_MATRIX_Y(:) ...
                                      ; incrementExtrinsicMatrix_x(:) ...
                                      ; incrementExtrinsicMatrix_y(:) ]';
                dependencyVariables = find(dependencyVariables);
            else
                dependencyVariables = [];
            end
        end
        
        function [dependency_x,dependency_y] = intrisicParametersDependencies()
            %
            % Obtain intrinsic parameters dependencies.
            %
            INTRINSIC_MATRIX_X = [1 1 1;0 0 0;0 0 0];
            INTRINSIC_MATRIX_Y = [0 0 0;1 1 1;0 0 0];
            INCREMENT_INTRINSIC_MATRIX_X = [1 1 1;0 0 0;0 0 0];
            INCREMENT_INTRINSIC_MATRIX_Y = [0 0 0;1 1 1;0 0 0];
            INCREMENT_EXTRINSIC_MATRIX_X = [1 1 1 1;0 0 0 0;0 0 0 0];
            INCREMENT_EXTRINSIC_MATRIX_Y = [0 0 0 0;1 1 1 1;0 0 0 0];
            
            % Obtain dependencies
            dependency_x = [ INTRINSIC_MATRIX_X(:) ...
                           ; INCREMENT_INTRINSIC_MATRIX_X(:) ...
                           ; INCREMENT_INTRINSIC_MATRIX_X(:) ...
                           ; INCREMENT_EXTRINSIC_MATRIX_X(:) ...
                           ; INCREMENT_EXTRINSIC_MATRIX_X(:) ]';
            dependency_y = [ INTRINSIC_MATRIX_Y(:) ...
                           ; INCREMENT_INTRINSIC_MATRIX_Y(:) ...
                           ; INCREMENT_INTRINSIC_MATRIX_Y(:) ...
                           ; INCREMENT_EXTRINSIC_MATRIX_Y(:) ...
                           ; INCREMENT_EXTRINSIC_MATRIX_Y(:) ]';
            dependency_x = find(dependency_x);
            dependency_y = find(dependency_y);
        end
        
        function [intrinsicMatrix,rmse,sse,residuals] = IntrinsicsFromHomographies( homographies ...
                                                                                  , cameras_xy )
            %
            % Estimate intrinsic matrix that characterizes the camera array
            % from homographies considering rectangular sampling optionally
            % removing the principal point shift.
            %
            % This estimation can be applied to either microlens or
            % viewpoint cameras.
            %
            % INPUTS:
            %   1. homographies - collection of camera array homographies.
            %   2. cameras_xy - Indicate cameras to perform intrinsics
            %   estimation.
            %
            narginchk(1,2);
            
            if nargin <= 1
                cameras_xy  = image.Pixel([1,0,2,0;0,1,0,2],false);
            end
            cameras_xy = utils.enums.Classes.PIXEL().convert(cameras_xy);
            
            DEGREES_FREEDOM = 11;
            
            % Obtain observation matrix based on homographies for x = y = 0
            % for x = 1 and y = 0, and for x = 0 and y = 1. This allows to 
            % obtain a symmetric positive definite matrix B that 
            % corresponds to 
            %        B = D_0 + x * E_x + y * E_y + x^2 * F_x + y^2 * F_y = 
            %          = (K_0 + C_xy)^-T * (K_0 + C_xy)^T
            % The matrix B has 11 degrees of freedom if we do not consider 
            % skew.
            iRow = 0;
            dltMatrix   = zeros(2 * cameras_xy.numberVectors * length(homographies),DEGREES_FREEDOM);
            for iHomography = 1:length(homographies)
                % Camera x = y = 0
                % Obtain homography matrix and corresponding entries
                homography = homographies(iHomography).homographyMatrix;
                h11        = homography(1,1);
                h12        = homography(1,2);
                h21        = homography(2,1);
                h22        = homography(2,2);
                h31        = homography(3,1);
                h32        = homography(3,2);
                
                % Define observation matrix entries
                dltMatrix(iRow + 1,1:5) = [ h11 * h12 ...
                                          , h11 * h32 + h12 * h31 ...
                                          , h21 * h22 ...
                                          , h21 * h32 + h22 * h31 ...
                                          , h31 * h32 ];
                dltMatrix(iRow + 2,1:5) = [ h11^2 - h12^2 ...
                                          , 2 * ( h11 * h31 - h12 * h32 ) ...
                                          , h21^2 - h22^2 ...
                                          , 2 * ( h21 * h31 - h22 * h32 ) ...
                                          , h31^2 - h32^2 ];
                iRow = iRow + 2;

                % Obtain the remainder homography matrices
                homographyMatrices = homographies(iHomography).obtainHomographyMatrices(cameras_xy);

                for iCamera = 1:cameras_xy.numberVectors
                    % Obtain homography matrix and corresponding entries
                    camera_xy  = cameras_xy.obtainVectors(iCamera);
                    homography = homographyMatrices(iCamera).homographyMatrix;
                    h11        = homography(1,1);
                    h12        = homography(1,2);
                    h21        = homography(2,1);
                    h22        = homography(2,2);
                    h31        = homography(3,1);
                    h32        = homography(3,2);

                    % Define observation matrix entries
                    dltMatrix(iRow + 1,1:11) = [ h11 * h12 ...
                                               , h11 * h32 + h12 * h31 ...
                                               , h21 * h22 ...
                                               , h21 * h32 + h22 * h31 ...
                                               , h31 * h32 ...
                                               , camera_xy.u   * (h11 * h32 + h12 * h31) ...
                                               , camera_xy.u   * (h31 * h32) ...
                                               , camera_xy.v   * (h21 * h32 + h22 * h31) ...
                                               , camera_xy.v   * (h31 * h32) ...
                                               , camera_xy.u^2 * (h31 * h32) ...
                                               , camera_xy.v^2 * (h31 * h32) ];
                    dltMatrix(iRow + 2,1:11) = [ h11^2 - h12^2 ...
                                               , 2 * ( h11 * h31 - h12 * h32 ) ...
                                               , h21^2 - h22^2 ...
                                               , 2 * ( h21 * h31 - h22 * h32 ) ...
                                               , h31^2 - h32^2 ...
                                               , camera_xy.u * 2 * (h11 * h31 - h12 * h32) ...
                                               , camera_xy.u     * (h31^2 - h32^2) ...
                                               , camera_xy.v * 2 * (h21 * h31 - h22 * h32) ...
                                               , camera_xy.v     * (h31^2 - h32^2) ...
                                               , camera_xy.u^2   * (h31^2 - h32^2) ...
                                               , camera_xy.v^2   * (h31^2 - h32^2) ];
                    iRow = iRow + 2;
                end
            end
            
            % The matrix B corresponds to the null space of the observation
            % matrix
            solver = optimization.Linear();
            solver.observationMatrix      = dltMatrix;
            [nullSpaceSolution,errorFlag] = solver.solve();
            
            % Obtain residuals and errors for least squares problem
            if nargout > 1
                residuals = abs(dltMatrix * nullSpaceSolution);
                sse  = sum(residuals.^2);               % Sum of squared errors
                rmse = sqrt(sse / length(residuals));   % Root mean squared error
            end                

            % See if solution is correct, otherwise consider the solution
            % with 10 degrees of freedom assuming that the elements b11 and
            % b22 are equal
            decodeSolution = 11;
            if errorFlag < 0 || nullSpaceSolution(1) * nullSpaceSolution(5) <= 0
                % Obtain observation matrix
                dltMatrix_DOF10      = dltMatrix(:,[1,2,4:end]);
                dltMatrix_DOF10(:,1) = dltMatrix_DOF10(:,1) + dltMatrix(:,3);
                dltMatrix_DOF10      = dltMatrix_DOF10 ./ max(abs(dltMatrix_DOF10(:)));
                
                % Obtain solution with 10 degrees of freedom
                solver.observationMatrix = dltMatrix_DOF10;
                [nullSpaceSolution_DOF10,errorFlag] = solver.solve();
                
                % See if solution is correct, otherwise consider the
                % solution with 9 degrees of freedom assuming that the
                % elements b13 and b23 are equal
                decodeSolution = 9;
                if errorFlag < 0 || nullSpaceSolution_DOF10(1) * nullSpaceSolution_DOF10(4) <= 0
                    % Obtain observation matrix
                    dltMatrix_DOF9      = dltMatrix_DOF10(:,[1,2,4:end]);
                    dltMatrix_DOF9(:,2) = dltMatrix_DOF9(:,2) + dltMatrix_DOF10(:,3);

                    % Obtain solution with 9 degrees of freedom
                    solver.observationMatrix = dltMatrix_DOF9;
                    [nullSpaceSolution_DOF9,errorFlag] = solver.solve();

                    % See if solution is correct.
                    % 9 degrees of freedom solution
                    if errorFlag == 0 && nullSpaceSolution_DOF9(1) * nullSpaceSolution_DOF9(3) > 0
                        decodeSolution = 9;
                    % 10 degrees of freedom
                    elseif nullSpaceSolution_DOF10(1) * nullSpaceSolution_DOF10(4) > 0
                        decodeSolution = 10;
                    % 9 degrees of freedom
                    elseif nullSpaceSolution_DOF9(1) * nullSpaceSolution_DOF9(3) > 0
                        decodeSolution = 9;
                    % 11 degrees of freedom
                    else 
                        decodeSolution = 11;
                    end
                end
            end
            
            % 10 degrees of freedom
            if decodeSolution == 8
                nullSpaceSolution = [ nullSpaceSolution_DOF10(1:2) ...
                                    ; nullSpaceSolution_DOF10(1) ...
                                    ; nullSpaceSolution_DOF10(3:end) ];
            % 9 degrees of freedom
            elseif decodeSolution == 7
                nullSpaceSolution = [ nullSpaceSolution_DOF9(1:2) ...
                                    ; nullSpaceSolution_DOF9(1:2) ...
                                    ; nullSpaceSolution_DOF9(3:end) ];
            end
            
            % Obtain entries from matrix D_0 and E_ij and assemble matrices
            d11   = nullSpaceSolution(1); 
            d13   = nullSpaceSolution(2);
            d22   = nullSpaceSolution(3);
            d23   = nullSpaceSolution(4);
            d33   = nullSpaceSolution(5);
			e13_x = nullSpaceSolution(6);
			e33_x = nullSpaceSolution(7);
			e23_y = nullSpaceSolution(8);
			e33_y = nullSpaceSolution(9);
			f33_x = nullSpaceSolution(10);
			f33_y = nullSpaceSolution(11);
            intrinsicLikeMatrix_xy_0 = [d11 0 d13; 0 d22 d23; d13 d23 d33];
            incrementLikeMatrix_x    = [0 0 e13_x; 0 0 0; e13_x 0 e33_x];
            incrementLikeMatrix_x2   = [0 0 0; 0 0 0; 0 0 f33_x];
            incrementLikeMatrix_y    = [0 0 0; 0 0 e23_y; 0 e23_y e33_y];
            incrementLikeMatrix_y2   = [0 0 0; 0 0 0; 0 0 f33_y];
            
            % Obtain intrinsic like matrix for camera index different than
            % zero
            camera_xy = image.Pixel([max(cameras_xy.u);max(cameras_xy.v)],false);
            intrinsicLikeMatrix_xy_1 = intrinsicLikeMatrix_xy_0 ...
                                     + camera_xy.u   * incrementLikeMatrix_x ...
                                     + camera_xy.u^2 * incrementLikeMatrix_x2 ...
                                     + camera_xy.v   * incrementLikeMatrix_y ...
                                     + camera_xy.v^2 * incrementLikeMatrix_y2;
            
            % Obtain intrinsic matrix using Cholesky decomposition. The
            % Cholesky decomposition gives the inverse of the intrinsic
            % matrix since it returns the matrix C from C^T * C = A.
            % This problem has two solutions (x and -x). Choose the
            % solution that gives positive diagonal elements. D_0 matrix is
            % positive definite.
            try
                intrinsicMatrixData_xy_0 = inv(chol(intrinsicLikeMatrix_xy_0,'upper'));
            catch error
                if strcmp(error.identifier,'MATLAB:posdef') > 0
                    intrinsicMatrixData_xy_0 = inv(chol(-intrinsicLikeMatrix_xy_0,'upper'));
                else
                    error.rethrow;
                end
            end
            
            % Intrinsic matrix data is defined up to a scale factor, so
            % normalize intrinsic matrix with k33 = 1.
            intrinsicMatrixData_xy_0 = intrinsicMatrixData_xy_0 ./ intrinsicMatrixData_xy_0(3,3);
            
            % Obtain intrinsic matrix for camera with coordinates u = v = 2
            obtainedIntrinsicMatrix = true;
            try
                intrinsicMatrixData_xy_1 = inv(chol(intrinsicLikeMatrix_xy_1,'upper'));
            catch error
                if strcmp(error.identifier,'MATLAB:posdef') > 0
                    try
                        intrinsicMatrixData_xy_1 = inv(chol(-intrinsicLikeMatrix_xy_1,'upper'));
                    catch error
                        if strcmp(error.identifier,'MATLAB:posdef') > 0
                            % Obtain principal point shift using ratios of entries
                            obtainedIntrinsicMatrix  = false;
                            principalPointShift_x    = -e13_x / d11;
                            principalPointShift_y    = -e23_y / d22;
                            intrinsicMatrixData_xy_1 = intrinsicMatrixData_xy_0;
                        else
                            error.rethrow;
                        end
                    end
                else
                    error.rethrow;
                end
            end
            
            % Intrinsic matrix data is defined up to a scale factor, so
            % normalize intrinsic matrix with k33 = 1.
            intrinsicMatrixData_xy_1 = intrinsicMatrixData_xy_1 ./ intrinsicMatrixData_xy_1(3,3);

            % Obtain mean value for scaling factor since this should be the
            % same among different viewpoint images.
            intrinsicMatrixData_xy_0(1,1) = (intrinsicMatrixData_xy_0(1,1) + intrinsicMatrixData_xy_1(1,1)) / 2;
            intrinsicMatrixData_xy_0(2,2) = (intrinsicMatrixData_xy_0(2,2) + intrinsicMatrixData_xy_1(2,2)) / 2;
            
            % Obtain increment matrix data
            incrementMatrixData = intrinsicMatrixData_xy_1 - intrinsicMatrixData_xy_0;
            
            % Let us isolate the contribution of each coordinate (x,y).
            % These contributions correspond to the principal point shift
            incrementMatrixData_x      = zeros(3,3);
            incrementMatrixData_x(1,3) = incrementMatrixData(1,3) ./ camera_xy.u;
            incrementMatrixData_y      = zeros(3,3);
            incrementMatrixData_y(2,3) = incrementMatrixData(2,3) ./ camera_xy.v;
            if obtainedIntrinsicMatrix == false
                incrementMatrixData_x(1,3) = principalPointShift_x;
                incrementMatrixData_y(2,3) = principalPointShift_y;
            end                
            
            % Define intrinsic matrix for array
            intrinsicMatrix = abstract.TemplatePlenopticProjectionMatrix();
            intrinsicMatrix.intrinsicMatrix = intrinsicMatrixData_xy_0;
            intrinsicMatrix.incrementIntrinsicMatrix_x = incrementMatrixData_x;
            intrinsicMatrix.incrementIntrinsicMatrix_y = incrementMatrixData_y;
        end
        
        function [projectionMatrices,rmse,sse,residuals] = ExtrinsicMatrixFromHomographies( intrinsics ...
                                                                                          , homographies )
            %
            % Estimate extrinsic matrix that characterizes the camera array
            % from homographies considering rectangular sampling.
            %
            % INPUTS:
            %   1. intrinsics   - intrinsic matrix obtained from
            %   homographies.
            %   2. homographies - collection of camera array homographies.
            %
            narginchk(2,2);
            
            VARIABLE_EXTRINSIC_DEGREES_FREEDOM   = 9;
            ADDITIONAL_EXTRINSIC_DEGREES_FREEDOM = 2;

            % Define cameras to be used
            camera_xy_1  = image.Pixel([1;1],false);

            % Obtain intrinsic matrix for camera x = y = 0, and for x = y =
            % 1.
            intrinsicMatrixData_xy_0 = intrinsics.intrinsicMatrix;
            intrinsicMatrixData_xy_1 = intrinsics.obtainProjectionMatrices(camera_xy_1).intrinsicMatrix;
            
            % Obtain extrinsic parameters from intrinsic matrix and 
            % homography matrix. Remember that the extrinsic parameters
            % between the several cameras differ only on the translation
            % components, namely in (x,y) (baseline). 
            % Vectorization puts 1st column, 2nd column and the 3rd
            % column in an unique vector.

            % Estimate extrinsic parameters for each homography and 
            % obtain projection matrix
            residuals          = [];
            projectionMatrices = [];
            for iHomography = 1:length(homographies)
                homographyMatrixData = [ homographies(iHomography).homographyMatrix(:) ...
                                       ; homographies(iHomography).obtainHomographyMatrices(camera_xy_1).homographyMatrix(:) ];
                intrinsicMatrixData  = zeros(length(homographyMatrixData), VARIABLE_EXTRINSIC_DEGREES_FREEDOM ...
                                                                         + ADDITIONAL_EXTRINSIC_DEGREES_FREEDOM );

                % Camera x = y = 0
                intrinsicMatrixData(01:03,01:03) = intrinsicMatrixData_xy_0;
                intrinsicMatrixData(04:06,04:06) = intrinsicMatrixData_xy_0;
                intrinsicMatrixData(07:09,07:09) = intrinsicMatrixData_xy_0;

                % Camera x = y = 1
                intrinsicMatrixData(10:12,01:03) = intrinsicMatrixData_xy_1;
                intrinsicMatrixData(13:15,04:06) = intrinsicMatrixData_xy_1;
                intrinsicMatrixData(16:18,07:09) = intrinsicMatrixData_xy_1;
                intrinsicMatrixData(16:18,10:11) = intrinsicMatrixData_xy_1(:,1:2);

                % Obtain extrinsic matrix data
                extrinsicMatrixData = mldivide( intrinsicMatrixData, homographyMatrixData );

                if nargout > 1
                    residuals = cat(1,residuals,abs(intrinsicMatrixData * extrinsicMatrixData) - homographyMatrixData);
                end
                
                % Isolate each extrinsic parameter and obtain scale factor
                % for correcting extrinsic parameters estimate. The scale
                % factor comes from the intrinsic matrix and homography
                % matrix being defined up to a scale factor
                rotationMatrix_r1      = extrinsicMatrixData(01:03);
                rotationMatrix_r2      = extrinsicMatrixData(04:06);
                translationVector_xy_0 = extrinsicMatrixData(07:09);
                baseline_xy_1          = extrinsicMatrixData(10:11);
                translationVector_xy_1 = translationVector_xy_0 + [baseline_xy_1;0];

                % Obtain rotation matrix using 2 vectors and scale
                % translation
                [rotationMatrix,scaleFactor] = math.RotationMatrix.RotationMatrixFrom2Vectors( rotationMatrix_r1 ...
                                                                                             , rotationMatrix_r2 );
                translationVector_xy_0 = math.Point(translationVector_xy_0 ./ scaleFactor,false);
                translationVector_xy_1 = math.Point(translationVector_xy_1 ./ scaleFactor,false);

                % Correct extrinsic parameters to have the z-coordinate
                % positive.
                if translationVector_xy_0.z < 0
                    rotationMatrix.data    = [ -rotationMatrix.data(:,1:2) ...
                                             ,  rotationMatrix.data(:,3) ];
                    translationVector_xy_0 = -translationVector_xy_0.data;
                    translationVector_xy_1 = -translationVector_xy_1.data;
                end

                % Camera x = y = 0
                % Create extrinsic matrix
                extrinsic_xy_0 = camera.models.ExtrinsicMatrix();
                extrinsic_xy_0.translationVector = translationVector_xy_0;
                extrinsic_xy_0.rotationMatrix    = rotationMatrix;

                % Camera x = y = 1
                % Create extrinsic matrix
                extrinsic_xy_1 = camera.models.ExtrinsicMatrix();
                extrinsic_xy_1.translationVector = translationVector_xy_1;
                extrinsic_xy_1.rotationMatrix    = rotationMatrix;

                % Obtain increment extrinsic matrix
                incrementExtrinsic = extrinsic_xy_1.extrinsicMatrix - extrinsic_xy_0.extrinsicMatrix;

                % Isolate the contribution of each coordinate
                incrementExtrinsic_x = zeros(3,4);
                incrementExtrinsic_y = zeros(3,4);
                incrementExtrinsic_x(1,:) = incrementExtrinsic(1,:);
                incrementExtrinsic_y(2,:) = incrementExtrinsic(2,:);

                % Create camera projection matrix array
                cameraArray = intrinsics;
                cameraArray.extrinsicMatrix = extrinsic_xy_0.extrinsicMatrix;
                cameraArray.incrementExtrinsicMatrix_x = incrementExtrinsic_x;
                cameraArray.incrementExtrinsicMatrix_y = incrementExtrinsic_y;

                % Update projection matrices
                projectionMatrices = cat(2,projectionMatrices,cameraArray);
            end
            
            if nargout > 1
                sse  = sum(residuals.^2);               % Sum of squared errors
                rmse = sqrt(sse / length(residuals));   % Root mean squared error
            end
        end
    end
end
