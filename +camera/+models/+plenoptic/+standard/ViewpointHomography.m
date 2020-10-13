classdef ViewpointHomography < abstract.TemplatePlenopticHomography
    %VIEWPOINTHOMOGRAPHY
    %   Viewpoint homography matrix instance
    
    properties (Dependent)
        incrementMatrix_i       % Increment matrix data for i-coordinate
        incrementMatrix_j       % Increment matrix data for j-coordinate
    end
    
    methods
        function self = ViewpointHomography(varargin)
            %
            % Create viewpoint homography instance.
            %
            % INPUTS:
            %   1. homographyMatrix  - homography matrix data
            %   2. incrementMatrix_i - increment matrix data for
            %   i-coordinate.
            %   3. incrementMatrix_j - increment matrix data for
            %   j-coordinate.
            %
            narginchk(0,3);
            
            self = self@abstract.TemplatePlenopticHomography;

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
        
        function increment_i = get.incrementMatrix_i(self)
            increment_i = self.incrementMatrix_x;
        end
        
        function increment_j = get.incrementMatrix_j(self)
            increment_j = self.incrementMatrix_y;
        end
    end
    
    methods
        function homographyMatrices = obtainHomographyMatrices(self,pixels_ij)
            %
            % Obtain homography matrices for viewpoint cameras.
            %
            % INPUTS:
            %   1. pixels_ij - pixel indices to obtain homography matrices.
            %      Each pixel should be provided in different columns.
            %
            narginchk(1,2);
            
            % If pixels are not provided, obtain a regular sampling of
            % pixels
            if nargin <= 1
                [pixels_i,pixels_j] = meshgrid(1:4:10,1:4:10);
                pixels_ij = [pixels_i(:),pixels_j(:)]';
            end
            
            % Obtain homography matrices
            homographyMatrices = self.obtainHomographyMatrices@abstract.TemplatePlenopticHomography(pixels_ij);
        end
        
        function homography = obtainHomographyMatrixForTargetCamera( self, sourcePixel_ij ...
                                                                   , targetPixel_ij, targetHomographyArray )
            %
            % Obtain homography matrix between viewpoint cameras source and
            % target viewpoint cameras.
            %
            % INPUTS:
            %   1. sourcePixel_ij - coordinates of the source viewpoint
            %   camera. Default is viewpoint camera i = j = 5.
            %   2. targetPixel_ij - coordinates of the target viewpoint
            %   camera. Default is (+1,+1) viewpoint camera relatively to
            %   the source viewpoint camera.
            %   3. targetHomographyArray - target camera array homography
            %   object. This assumes that the current class is the source
            %   homography camera array. If not provided, the homography is
            %   computed considering a unique pose.
            %
            narginchk(1,4);
            
            % If pixels are not provided, consider the viewpoint camera i =
            % j = 5
            if nargin <= 1
                sourcePixel_ij = [5;5];
            end
            sourcePixel_ij = utils.enums.Classes.PIXEL().convert(sourcePixel_ij);
            
            % If pixels are not provided, consider the next viewpoint
            % camera relatively to the source viewpoint camera.
            if nargin <= 2
                targetPixel_ij = sourcePixel_ij.data + 1;
            end
            targetPixel_ij = utils.enums.Classes.PIXEL().convert(targetPixel_ij);
            
            if nargin <= 3
                targetHomographyArray = self;
            end
            
            % Obtain homography matrices
            homography = self.obtainHomographyMatrixForTargetCamera@abstract.TemplatePlenopticHomography( ...
                            sourcePixel_ij, targetPixel_ij, targetHomographyArray );
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
            
            rays_ijkl = self.project@abstract.TemplatePlenopticHomography( worldPoints ...
                                                                         , pixels_ij );
        end
        
        function worldPoint = reconstruct(self,rays_ijkl)
            %
            % Reconstruct world point from pixels in image sensor. This
            % will give the (x,y)-coordinates of the world point.
            %
            % INPUTS: 
            %   1. rays_ijkl - image ray coordinates for a given point.
            % 
            narginchk(2,2);
            
            CAMERA_TYPE = camera.models.plenoptic.enums.CameraTypes.VIEWPOINT(); 
            worldPoint  = self.reconstruct@abstract.TemplatePlenopticHomography( rays_ijkl ...
                                                                               , CAMERA_TYPE );
        end
    end

    methods (Static)
        function self = ViewpointHomographyFromStandardPlenopticCamera(plenopticCamera)
            %
            % Obtain parameters to describe a viewpoint homography from 
            % a standard plenoptic camera. In this situation we are 
            % considering that (i,j) is fixed and (k,l) can have different 
            % values.
            %
            narginchk(1,1);
            
            % Convert to plenoptic camera instance
            plenopticCamera = utils.enums.Classes.STANDARD_PLENOPTIC_CAMERA().convert(plenopticCamera);
            
            % Obtain viewpoint camera array from standard plenoptic camera
            viewpointCameraArray = camera.models.plenoptic.standard.ViewpointProjectionMatrix.ViewpointProjectionMatrixFromStandardPlenopticCamera(plenopticCamera);
            
            % Obtain projection matrix for viewpoint i = j = 0 and for i =
            % j = 1.
            pixels_ij          = image.Pixel([0,1;0,1],false);
            projectionMatrices = viewpointCameraArray.obtainProjectionMatrices(pixels_ij);
            
            % Homography matrix for i = j = 0
            pinhole    = projectionMatrices(1).projectionMatrix;
            homography = camera.models.pinhole.Homography.HomographyFromPinhole(pinhole);
            
            % Homography matrix for i = j = 1
            pinhole   = projectionMatrices(2).projectionMatrix;
            increment = camera.models.pinhole.Homography.HomographyFromPinhole(pinhole);

            % Obtain increment matrix
            incrementMatrix = increment.homographyMatrix - homography.homographyMatrix;
            
            % Let us isolate the contribution of each coordinate (i,j)
            incrementMatrix_i      = zeros(3,3);
            incrementMatrix_i(1,:) = incrementMatrix(1,:);
            incrementMatrix_j      = zeros(3,3);
            incrementMatrix_j(2,:) = incrementMatrix(2,:);
            
            % Create viewpoint homography instance
            self = camera.models.plenoptic.standard.ViewpointHomography();
            self.homographyMatrix  = homography.homographyMatrix;
            self.incrementMatrix_x = incrementMatrix_i;
            self.incrementMatrix_y = incrementMatrix_j;
        end
        
        function [self,rmse,sse,residuals] = HomographyFromPointCorrespondences( worldPoints, imageRays ...
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
            % Obtain information for viewpoint homography
            if nargout >= 2
                [temp,rmse,sse,residuals] = HomographyFromPointCorrespondences@abstract.TemplatePlenopticHomography( ...
                                                    worldPoints, imageRays ...
                                                  , observationMatrixMethod ...
                                                  , CAMERA_TYPE );
            else
                temp = HomographyFromPointCorrespondences@abstract.TemplatePlenopticHomography( ...
                                                    worldPoints, imageRays ...
                                                  , observationMatrixMethod ...
                                                  , CAMERA_TYPE );
            end
            
            % Create viewpoint homography class
            self = camera.models.plenoptic.standard.ViewpointHomography();
            self.homographyMatrix  = temp.homographyMatrix;
            self.incrementMatrix_x = temp.incrementMatrix_x;
            self.incrementMatrix_y = temp.incrementMatrix_y;
        end
        
        function [self,rmse,sse,residuals] = ViewpointHomographyFromPointLineCorrespondences( worldPoints ...
                                                                                            , viewpointsLines  )
            %
            % Obtain parameters to describe a viewpoint camera array 
            % homography using direct linear transform from point and line 
            % correspondences.
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
            PARAMETERS_TO_ESTIMATE = 15;

            % Transform to cell
            if ~iscell(worldPoints)
                worldPoints = num2cell(worldPoints);
            end
            
            % Transform world points and image rays to vector class 
            worldPoints     = cellfun(@(x) utils.enums.Classes.POINT().convert(x),worldPoints,'UniformOutput',false);
            viewpointsLines = utils.enums.Classes.VECTOR().convert(viewpointsLines);

            % Since for the world points we are going to consider only
            % the (x,y) coordinates, these can be considered as pixels
            worldPixels = cellfun(@(p) image.Pixel([p.x;p.y],false),worldPoints,'UniformOutput',false);

            % Transform points to homogeneous coordinates
            worldPixels = cellfun(@(x) x.setHomogeneousCoordinates(), worldPixels,'UniformOutput',false);

            % Initialize dlt matrix
            numberPoints = [worldPixels{:}];
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
                linePoints = worldPixels{iLine};
                
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
            
            % Obtain viewpoint homography matrix. The homography matrix has
            % size (3 x 3)
            homographyMatrixData = reshape(nullSpaceSolution(01:09) , 3, 3)';

            % Obtain matrix with same size of homography matrix to obtain
            % homography matrices from other viewpoints.
            incrementMatrixData = [reshape(nullSpaceSolution(10:end), 3, 2)';zeros(1,3)];
            
            % Ensure the homography has the same signal
            homographyMatrixScaleFactor = sign(homographyMatrixData(3,3));
            homographyMatrixData        = homographyMatrixData ./ homographyMatrixScaleFactor;
            incrementMatrixData         = incrementMatrixData  ./ homographyMatrixScaleFactor;

            % Let us isolate the contribution of each coordinate (i,j)
            incrementMatrixData_i      = zeros(3,3);
            incrementMatrixData_i(1,:) = incrementMatrixData(1,:);
            incrementMatrixData_j      = zeros(3,3);
            incrementMatrixData_j(2,:) = incrementMatrixData(2,:);
            
            % Create viewpoint homography instance
            self = camera.models.plenoptic.standard.ViewpointHomography();
            self.homographyMatrix  = homographyMatrixData;
            self.incrementMatrix_x = incrementMatrixData_i;
            self.incrementMatrix_y = incrementMatrixData_j;
        end
    end
end

