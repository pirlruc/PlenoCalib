classdef ViewpointPinhole < abstract.TemplatePlenopticPinhole
    %VIEWPOINTPINHOLE
    %   Viewpoint pinhole camera instance
    
    properties (Dependent)
        incrementMatrix_i       % Increment matrix data for i-coordinate
        incrementMatrix_j       % Increment matrix data for j-coordinate
    end
    
    methods
        function self = ViewpointPinhole(varargin)
            %
            % Create viewpoint pinhole instance.
            %
            % INPUTS:
            %   1. projectionMatrix - projection matrix data
            %   2. incrementMatrix_i - increment matrix data for
            %   i-coordinate.
            %   3. incrementMatrix_j - increment matrix data for
            %   j-coordinate.
            %
            narginchk(0,3);
            
            self = self@abstract.TemplatePlenopticPinhole;

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
        
        function increment_i = get.incrementMatrix_i(self)
            increment_i = self.incrementMatrix_x;
        end
        
        function increment_j = get.incrementMatrix_j(self)
            increment_j = self.incrementMatrix_y;
        end
    end
    
    methods
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
            projectionMatrices = self.obtainProjectionMatrices@abstract.TemplatePlenopticPinhole(pixels_ij);
        end
        
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
            
            rays_ijkl = self.project@abstract.TemplatePlenopticPinhole( worldPoints ...
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
            rays_ijkl   = self.reconstruct@abstract.TemplatePlenopticPinhole( rays_ijkl ...
                                                                            , CAMERA_TYPE );
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
            
            % Obtain plenoptic pinhole structure
            temp = abstract.TemplatePlenopticPinhole.PinholeFromProjectionMatrix(projectionMatrix);
            
            % Create structure for viewpoint pinhole
            self = camera.models.plenoptic.standard.ViewpointPinhole();
            self.projectionMatrix  = temp.projectionMatrix;
            self.incrementMatrix_x = temp.incrementMatrix_x;
            self.incrementMatrix_y = temp.incrementMatrix_y;
        end
        
        function [self,rmse,sse,residuals] = PinholeFromPointCorrespondences( worldPoints, imageRays ...
                                                                            , observationMatrixMethod )
            %
            % Obtain parameters to describe a viewpoint camera array using 
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
            %
            narginchk(2,3);
            
            if nargin <= 2
                observationMatrixMethod = camera.models.pinhole.enums.ObservationMatrixMethods.CROSS_PRODUCT();
            end
            
            CAMERA_TYPE = camera.models.plenoptic.enums.CameraTypes.VIEWPOINT();
            
            % Obtain information for viewpoint pinhole
            if nargout >= 2
                [temp,rmse,sse,residuals] = PinholeFromPointCorrespondences@abstract.TemplatePlenopticPinhole( ...
                                                    worldPoints, imageRays ...
                                                  , observationMatrixMethod ...
                                                  , CAMERA_TYPE );
            else
                temp = PinholeFromPointCorrespondences@abstract.TemplatePlenopticPinhole( ...
                                                    worldPoints, imageRays ...
                                                  , observationMatrixMethod ...
                                                  , CAMERA_TYPE );
            end
            
            % Create viewpoint pinhole class
            self = camera.models.plenoptic.standard.ViewpointPinhole();
            self.projectionMatrix  = temp.projectionMatrix;
            self.incrementMatrix_x = temp.incrementMatrix_x;
            self.incrementMatrix_y = temp.incrementMatrix_y;
        end
        
        function [self,rmse,sse,residuals] = ViewpointPinholeFromPointLineCorrespondences( worldPoints ...
                                                                                         , viewpointsLines )
            %
            % Obtain parameters to describe a viewpoint camera array using 
            % direct linear transform from point and line correspondences.
            %
            % INPUTS:
            %   1. worldPoints - points defined in the world coordinate
            %   system. Since in a lightfield a line has several points and
            %   a point has several lines, the points must be given as a 
            %   cell. Each cell contains the points of a line.
            %   The points should not be given in homogeneous coordinates.
            %   2. viewpointsLines - line parameters identified on a given
            %   viewpoint image. This is a vector containing (i,j,a,b,c).
            %
            narginchk(2,2);
            
            % Number of unknown parameters
            PARAMETERS_TO_ESTIMATE = 20;

            % Transform to cell
            if ~iscell(worldPoints)
                worldPoints = num2cell(worldPoints);
            end
            
            % Transform world points and line parameters to vector class 
            worldPoints     = cellfun(@(x) utils.enums.Classes.POINT().convert(x),worldPoints,'UniformOutput',false);
            viewpointsLines = utils.enums.Classes.VECTOR().convert(viewpointsLines);

            % Transform points to homogeneous coordinates
            worldPoints = cellfun(@(x) x.setHomogeneousCoordinates(), worldPoints,'UniformOutput',false);

            % Initialize dlt matrix
            numberPoints = [worldPoints{:}];
            numberPoints = [numberPoints.numberVectors];
            dltMatrix    = zeros(sum(numberPoints),PARAMETERS_TO_ESTIMATE);

            % Construct direct linear transformation calibration matrix
            % considering (i,j,a,b,c) as the image line parameters and M 
            % as the point in world coordinates:
            %
            %   - Linear Solution
            %       A = [ a*M^T b*M^T c*M^T a*i*M^T b*j*M^T ]
            iRow = 0;
            for iLine = 1:viewpointsLines.numberVectors
                % Obtain line and point correspondences
                line       = viewpointsLines.obtainVectors(iLine);
                linePoints = worldPoints{iLine};
                
                dltMatrix(iRow + 1:iRow + linePoints.numberVectors,:) = ...
                            [ line.data(3) .* linePoints.data' ...
                            , line.data(4) .* linePoints.data' ...
                            , line.data(5) .* linePoints.data' ...
                            , line.data(3) .* line.data(1) .* linePoints.data' ...
                            , line.data(4) .* line.data(2) .* linePoints.data' ];
                        
                % Update row index
                iRow = iRow + linePoints.numberVectors;
            end
            % Since the number of observations is huge, perform 
            % QR-decomposition to obtain an orthogonal matrix (Q) and an 
            % upper triangular matrix (R)
            [~,upperTriangular] = qr(dltMatrix,0);
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
            homogeneousCoordinate = projectionMatrixData(3,:) * linePoints.data;
            if any(homogeneousCoordinate < 0)
                projectionMatrixData = -projectionMatrixData;
                incrementMatrixData  = -incrementMatrixData;
            end
            
            % Since projection matrix is defined up to a scale factor
            % normalize projection matrix with ||[p31 p32 p33]|| = 1.
            projectionMatrixScaleFactor = norm(projectionMatrixData(3,1:3));
            projectionMatrixData        = projectionMatrixData ./ projectionMatrixScaleFactor ;
            incrementMatrixData         = incrementMatrixData  ./ projectionMatrixScaleFactor;
            
            % Let us isolate the contribution of each coordinate (i,j)
            incrementMatrixData_i      = zeros(3,4);
            incrementMatrixData_i(1,:) = incrementMatrixData(1,:);
            incrementMatrixData_j      = zeros(3,4);
            incrementMatrixData_j(2,:) = incrementMatrixData(2,:);

            % Create viewpoint pinhole instance
            self = camera.models.plenoptic.standard.ViewpointPinhole();
            self.projectionMatrix  = camera.models.pinhole.Pinhole(projectionMatrixData);
            self.incrementMatrix_x = incrementMatrixData_i;
            self.incrementMatrix_y = incrementMatrixData_j;
        end
    end
end
