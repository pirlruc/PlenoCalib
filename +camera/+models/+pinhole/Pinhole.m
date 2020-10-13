classdef Pinhole
    %PINHOLE
    %   Pinhole camera model simulator for projection and calibration.
    
    properties
        projectionMatrix = camera.models.pinhole.ProjectionMatrix().projectionMatrix
                % Projection matrix that allows to project a point in the world coordinate 
                % system to the image plane
    end
    
    properties (Dependent)
        % Decomposition of projection matrix based on:
        %       R. Hartley, A. Zisserman, "Multiple View Geometry",
        %       Cambridge University Press 2000, page 150.
        extrinsicMatrix         % Extrinsic matrix from projection matrix
        intrinsicMatrix         % Intrinsic matrix from projection matrix
        principalAxis           % Principal axis
        projectionCenter        % Projection center associated with projection matrix
    end
    
    methods
        function self = Pinhole(varargin)
            %
            % Create pinhole camera model instance.
            %
            % INPUTS:
            %   1. projectionMatrix - projection matrix.
            %
            narginchk(0,1);

            if ~isempty(varargin)
                if nargin >= 1
                    self.projectionMatrix = varargin{1};
                end
            end
        end
        
        function self = set.projectionMatrix(self,newProjectionMatrixData)
            self.projectionMatrix = newProjectionMatrixData;
        end
        
        function matrix = get.extrinsicMatrix(self)
            % Decompose projection matrix 
            [~,rotationMatrix,translationVector] = camera.models.pinhole.utils.proj_decomp(self.projectionMatrix);
            
            % Create extrinsic matrix instance
            matrix = camera.models.ExtrinsicMatrix();
            matrix.rotationMatrix    = rotationMatrix;
            matrix.translationVector = translationVector;
        end
        
        function matrix = get.intrinsicMatrix(self)
            % Decompose projection matrix. The zeros are added to obtain a
            % matrix with the same structure of the object
            % camera.models.pinhole.IntrisicMatrix.
            matrix = camera.models.pinhole.utils.proj_decomp(self.projectionMatrix);
        end
        
        function axis = get.principalAxis(self)
            % Obtain principal axis from projection matrix. According to
            % the formulas presented in Hartley and Zisserman.
            
            % Obtain matrix that results from multiplying the intrinsics
            % matrix with the rotation matrix: M = K * R.
            matrix = self.intrinsicMatrix * self.extrinsicMatrix.rotationMatrix.data;
            
            % Obtain principal axis: det(M) * m_3. m_3 corresponds to the
            % 3rd row of matrix M. The axis is returned as a column vector.
            axis = det(matrix) * matrix(end,:)';
            
            % Normalize principal axis
            axis = axis ./ norm(axis);
        end
        
        function point = get.projectionCenter(self)
            % Obtain projection center. The projection center corresponds
            % to the null space of the projection matrix.
            [~,~,nullSpaceSolution] = svd(self.projectionMatrix);
            nullSpaceSolution       = nullSpaceSolution(:,end);
            
            % Solution is given in homogeneous coordinates
            point = math.Point( nullSpaceSolution(1:end-1) ./ nullSpaceSolution(end), false );
        end
    end
    
    methods
        function pixels = project(self,worldPoints,pixelsToNearestInteger)
            %
            % Obtain projection of points in the world coordinate system.
            %
            % INPUTS:
            %   1. worldPoints - points in the world coordinate system. 
            %      Each point should be defined in a different column. The 
            %      points should not be given in homogeneous coordinates.
            %   2. pixelsToNearestInteger - flag that indicates if pixels
            %      must be rounded to the nearest integer. The default is
            %      true.
            %
            narginchk(1,3);
            
            % If no points are provided, generate random points.
            if nargin <= 1
                worldPoints = rand(3,50);
            end
            
            % If round flag is not provided, assume that the pixels should
            % be returned rounded to the nearest integer.
            if nargin <= 2
                pixelsToNearestInteger = true;
            end
            
            % If points are provided in rows instead of columns, transpose
            % the points
            worldPoints = utils.enums.Classes.POINT().convert(worldPoints);

            % Define world points in homogeneous coordinates
            worldPoints = worldPoints.setHomogeneousCoordinates();
            
            % Obtain projected points
            pixels = image.Pixel(self.projectionMatrix * worldPoints.data,true);
            pixels = pixels.removeHomogeneousCoordinates();
            
            % Round pixels to nearest integer
            if pixelsToNearestInteger == true
                pixels.data = round(pixels.data);
            end
        end
    end
    
    methods (Static)
        function [self,rmse,sse,residuals] = PinholeFromPointCorrespondences( worldPoints,pixels ...
                                                                            , normalize ...
                                                                            , observationMatrixMethod )
            %
            % Create pinhole camera model instance. The projection matrix
            % is obtained from calibration using direct linear
            % transformation from point correspondences.
            %
            % INPUTS:
            %   1. worldPoints - points defined in the world coordinate
            %   system. Each point should be defined in a different column.
            %   The points should not be given in homogeneous coordinates.
            %   2. pixels    - projected points in the image plane
            %   corresponding to the world points given as input.
            %   3. normalize - normalize projection and world coordinates. 
            %   Default is true.
            %   4. observationMatrixMethod - method to construct
            %   observation matrix. Default is the cross product method.
            %
            narginchk(2,4);
            
            if nargin <= 2
                normalize = true;
            end
            
            if nargin <= 3
                observationMatrixMethod = camera.models.pinhole.enums.ObservationMatrixMethods.CROSS_PRODUCT();
            end
            
            % Transform world points and image pixels to vector class
            worldPoints = utils.enums.Classes.POINT().convert(worldPoints);
            pixels      = utils.enums.Classes.PIXEL().convert(pixels);
            
            % Normalize projections
            if normalize == true
                % Obtain normalization matrices
                pixelsNormalizationMatrix  = pixels.obtainNormalizationMatrix();
                pointsNormalizationMatrix  = worldPoints.obtainNormalizationMatrix();

                % Obtain normalized pixels and points
                pixels      = pixels.normalize;
                worldPoints = worldPoints.normalize;
            else
                pixelsNormalizationMatrix = eye(3,3);
                pointsNormalizationMatrix = eye(4,4);
            end
            
            % Transform points to homogeneous coordinates
            worldPoints = worldPoints.setHomogeneousCoordinates();
            pixels      = pixels.setHomogeneousCoordinates();
           
            % Construct direct linear transformation calibration matrix
            % considering (u,v) as the pixel coordinates and M as the point
            % in world coordinates:
            %
            %   - Linear Solution
            %       A = [ M^T 0 -u*M^T ; 0 M^T -v*M^T ]
            if observationMatrixMethod == camera.models.pinhole.enums.ObservationMatrixMethods.LINEAR_SYSTEM
                dltMatrix = [ worldPoints.data' ...
                            , zeros(worldPoints.numberVectors,worldPoints.numberComponents) ...
                            , -repmat(pixels.u',1,worldPoints.numberComponents) .* worldPoints.data' ...
                            ; zeros(worldPoints.numberVectors,worldPoints.numberComponents) ...
                            , worldPoints.data' ...
                            , -repmat(pixels.v',1,worldPoints.numberComponents) .* worldPoints.data' ];
            
            %   - Cross-Product
            %       A = [ M^T 0 -u*M^T ; 0 -M^T v*M^T ; -v*M^T u*M^T 0 ]
            elseif observationMatrixMethod == camera.models.pinhole.enums.ObservationMatrixMethods.CROSS_PRODUCT
                dltMatrix = [ worldPoints.data' ...
                            , zeros(worldPoints.numberVectors,worldPoints.numberComponents) ...
                            , -repmat(pixels.u',1,worldPoints.numberComponents) .* worldPoints.data' ...
                            ; zeros(worldPoints.numberVectors,worldPoints.numberComponents) ...
                            , -worldPoints.data' ...
                            , repmat(pixels.v',1,worldPoints.numberComponents) .* worldPoints.data' ...
                            ; -repmat(pixels.v',1,worldPoints.numberComponents) .* worldPoints.data' ...
                            , repmat(pixels.u',1,worldPoints.numberComponents) .* worldPoints.data' ...
                            , zeros(worldPoints.numberVectors,worldPoints.numberComponents) ];
            end
            
            % The projection matrix corresponds to the null space of this
            % matrix
            [~,~,eigenvectors] = svd(dltMatrix);
            nullSpaceSolution  = eigenvectors(:, end);
            
            % Obtain residuals and errors for least squares problem
            if nargout > 1
                residuals = abs(dltMatrix * nullSpaceSolution);
                sse  = sum(residuals.^2);               % Sum of squared errors
                rmse = sqrt(sse / length(residuals));   % Root mean squared error
            end                
            
            % Obtain projection matrix from null space solution
            projectionMatrixData = reshape( nullSpaceSolution ...
                                          , worldPoints.numberComponents ...
                                          , pixels.numberComponents)';
                                 
            % This problem has two solutions (x and -x). Choose the
            % solution that has the homogeneous coordinate positive.
            homogeneousCoordinate = projectionMatrixData(3,:) * worldPoints.data;
            if any(homogeneousCoordinate < 0)
                projectionMatrixData = -projectionMatrixData;
            end
            
            % Obtain unnormalized projection matrix
            projectionMatrixData = mldivide(pixelsNormalizationMatrix,projectionMatrixData * pointsNormalizationMatrix);
            
            % Since projection matrix is defined up to a scale factor
            % normalize projection matrix with ||[p31 p32 p33]|| = 1.
            projectionMatrixData = projectionMatrixData ./ norm(projectionMatrixData(3,1:3));
            
            % Create pinhole camera model instance with the projection
            % matrix obtained from calibration
            self = camera.models.pinhole.Pinhole(projectionMatrixData);
        end
        
        function [self,rmse,sse,residuals] = PinholeFromPointLineCorrespondences( worldPoints, lines )
            %
            % Create pinhole camera model instance. The projection matrix
            % is obtained from calibration using direct linear
            % transformation from point to line correspondences.
            %
            % INPUTS:
            %   1. worldPoints - points defined in the world coordinate
            %   system. Since a line has several points, the points must be
            %   given as a cell. Each cell contains the points of a line.
            %   The points should not be given in homogeneous coordinates.
            %   2. lines - line parameters obtained from image plane. The
            %   line parameters are given in a vector (a,b,c).
            %
            narginchk(2,2);
            
            % Transform to cell
            if ~iscell(worldPoints)
                worldPoints = {worldPoints};
            end
            
            % Transform world and line parameters to vector class
            worldPoints = cellfun(@(x) utils.enums.Classes.POINT().convert(x),worldPoints,'UniformOutput',false);
            lines       = utils.enums.Classes.VECTOR().convert(lines);
            
            % Transform points to homogeneous coordinates
            worldPoints = cellfun(@(x) x.setHomogeneousCoordinates(), worldPoints,'UniformOutput',false);

            % Initialize dlt matrix
            numberPoints = [worldPoints{:}];
            numberPoints = [numberPoints.numberVectors];
            dltMatrix    = zeros(sum(numberPoints),12);
            
            % Construct direct linear transformation calibration matrix
            % considering (a,b,c) as the image line parameters and M as the
            % point in world coordinates:
            %
            %   - Linear Solution
            %       A = [ a*M^T b*M^T c*M^T ]
            iRow = 0;
            for iLine = 1:lines.numberVectors
                % Obtain line and point correspondences
                line       = lines.obtainVectors(iLine);
                linePoints = worldPoints{iLine};
                
                dltMatrix(iRow + 1:iRow + linePoints.numberVectors,:) = ...
                            [ line.data(1) .* linePoints.data' ...
                            , line.data(2) .* linePoints.data' ...
                            , line.data(3) .* linePoints.data' ];
                        
                % Update row index
                iRow = iRow + linePoints.numberVectors;
            end
            
            % The projection matrix corresponds to the null space of this
            % matrix
            [~,~,eigenvectors] = svd(dltMatrix);
            nullSpaceSolution  = eigenvectors(:, end);
            
            % Obtain residuals and errors for least squares problem
            if nargout > 1
                residuals = abs(dltMatrix * nullSpaceSolution);
                sse  = sum(residuals.^2);               % Sum of squared errors
                rmse = sqrt(sse / length(residuals));   % Root mean squared error
            end                
            
            % Obtain projection matrix from null space solution, Projection
            % matrix corresponds to a 3 x 4 matrix
            projectionMatrixData = reshape( nullSpaceSolution, 4 ,3)';
                                 
            % This problem has two solutions (x and -x). Choose the
            % solution that has the homogeneous coordinate positive.
            homogeneousCoordinate = projectionMatrixData(3,:) * linePoints.data;
            if any(homogeneousCoordinate < 0)
                projectionMatrixData = -projectionMatrixData;
            end
            
            % Since projection matrix is defined up to a scale factor
            % normalize projection matrix with ||[p31 p32 p33]|| = 1.
            projectionMatrixData = projectionMatrixData ./ norm(projectionMatrixData(3,1:3));
            
            % Create pinhole camera model instance with the projection
            % matrix obtained from calibration
            self = camera.models.pinhole.Pinhole(projectionMatrixData);
        end
    end
end
