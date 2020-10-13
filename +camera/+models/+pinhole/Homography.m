classdef Homography
    %HOMOGRAPHY
    %   Homography estimation utilities.
    
    properties
        homographyMatrix = eye(3,3)         % Homography matrix 
    end
    
    methods
        function self = Homography(varargin)
            %
            % Create homography instance.
            %
            % INPUTS:
            %   1. homographyMatrix - homography matrix.
            %
            narginchk(0,1);

            if ~isempty(varargin)
                if nargin >= 1
                    self.homographyMatrix = varargin{1};
                end
            end
        end
    end
    
    methods (Static)
        function self = HomographyFromPinhole(pinhole)
            %
            % Obtain homography from calibrated pinhole camera.
            %
            % INPUTS:
            %   1. pinhole - pinhole camera model object.
            %
            narginchk(1,1);
            
            % Convert to pinhole object
            pinhole = utils.enums.Classes.PINHOLE().convert(pinhole);
            
            % Obtain homography matrix
            homographyMatrix = pinhole.intrinsicMatrix ...
                             * [ pinhole.extrinsicMatrix.rotationMatrix.data(:,1:2) ...
                               , pinhole.extrinsicMatrix.translationVector.data ];

            % Since homography matrix is defined up to a scale factor,
            % normalize homography matrix with entry h33 = 1.
            homographyMatrix = homographyMatrix ./ homographyMatrix(3,3);
            
            % Create homography matrix with the matrix obtained from 
            % calibration
            self = camera.models.pinhole.Homography(homographyMatrix);
        end
        
        function [self,rmse,sse,residuals] = HomographyFromPointCorrespondences( worldPoints,pixels ...
                                                                               , normalize ...
                                                                               , observationMatrixMethod ...
                                                                               , normalizeHomographies )
            %
            % Create homography matrix. The homography matrix is obtained 
            % from calibration using direct linear transformation from 
            % point correspondences.
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
            %   5. normalizeHomographies   - flag to indicate if
            %   homographies should be normalized to h33 = 1. Default is
            %   false.
            %
            narginchk(2,5);
            
            if nargin <= 2
                normalize = true;
            end
            
            if nargin <= 3
                observationMatrixMethod = camera.models.pinhole.enums.ObservationMatrixMethods.CROSS_PRODUCT();
            end
            
            if nargin <= 4
                normalizeHomographies = false;
            end
            
            % Transform world points and image pixels to vector class
            worldPoints = utils.enums.Classes.POINT().convert(worldPoints);
            pixels      = utils.enums.Classes.PIXEL().convert(pixels);
            
            % Since for the world points we are going to consider only
            % the (x,y) coordinates, these can be considered as pixels
            worldPixels = image.Pixel([worldPoints.x;worldPoints.y],false);

            % Normalize projections
            if normalize == true
                % Obtain normalization matrices
                pixelsNormalizationMatrix = pixels.obtainNormalizationMatrix();
                pointsNormalizationMatrix = worldPixels.obtainNormalizationMatrix();
                
                % Obtain normalized pixels and points
                pixels      = pixels.normalize;
                worldPixels = worldPixels.normalize;
            else
                pixelsNormalizationMatrix = eye(3,3);
                pointsNormalizationMatrix = eye(3,3);
            end
            
            % Transform points to homogeneous coordinates
            worldPixels     = worldPixels.setHomogeneousCoordinates();
            pixels          = pixels.setHomogeneousCoordinates();
            
            % Construct direct linear transformation homography matrix
            % considering (u,v) as the pixel coordinates and M as the point
            % in world coordinates using only (x,y) coordinates:
            %
            %   - Linear Solution
            %       A = [ M_xy^T 0 -u*M_xy^T ; 0 M_xy^T -v*M_xy^T ]
            if observationMatrixMethod == camera.models.pinhole.enums.ObservationMatrixMethods.LINEAR_SYSTEM
                dltMatrix = [ worldPixels.data' ...
                            , zeros(worldPixels.numberVectors,worldPixels.numberComponents) ...
                            , -repmat(pixels.u',1,worldPixels.numberComponents) .* worldPixels.data' ...
                            ; zeros(worldPixels.numberVectors,worldPixels.numberComponents) ...
                            , worldPixels.data' ...
                            , -repmat(pixels.v',1,worldPixels.numberComponents) .* worldPixels.data' ];
            
            %   - Cross-Product
            %       A = [ M_xy^T 0 -u*M_xy^T ; 0 -M_xy^T v*M_xy^T ; -v*M_xy^T u*M_xy^T 0 ]
            elseif observationMatrixMethod == camera.models.pinhole.enums.ObservationMatrixMethods.CROSS_PRODUCT
                dltMatrix = [ worldPixels.data' ...
                            , zeros(worldPixels.numberVectors,worldPixels.numberComponents) ...
                            , -repmat(pixels.u',1,worldPixels.numberComponents) .* worldPixels.data' ...
                            ; zeros(worldPixels.numberVectors,worldPixels.numberComponents) ...
                            , -worldPixels.data' ...
                            , repmat(pixels.v',1,worldPixels.numberComponents) .* worldPixels.data' ...
                            ; -repmat(pixels.v',1,worldPixels.numberComponents) .* worldPixels.data' ...
                            , repmat(pixels.u',1,worldPixels.numberComponents) .* worldPixels.data' ...
                            , zeros(worldPixels.numberVectors,worldPixels.numberComponents) ];
            end
            
            % The homography matrix corresponds to the null space of this
            % matrix
            [~,~,eigenvectors] = svd(dltMatrix);
            nullSpaceSolution  = eigenvectors(:, end);
            
            % Obtain residuals and errors for least squares problem
            if nargout > 1
                residuals = abs(dltMatrix * nullSpaceSolution);
                sse  = sum(residuals.^2);               % Sum of squared errors
                rmse = sqrt(sse / length(residuals));   % Root mean squared error
            end                
            
            % Obtain homography matrix from null space solution. The matrix
            % is 3 x 3.
            homographyMatrix = reshape( nullSpaceSolution, 3, 3)';
                                 
            % Obtain unnormalized projection matrix
            homographyMatrix = mldivide(pixelsNormalizationMatrix,homographyMatrix * pointsNormalizationMatrix);
            
            % Since homography matrix is defined up to a scale factor,
            % normalize homography matrix with entry h33 = 1.
            if normalizeHomographies == true
                homographyMatrixScaleFactor = homographyMatrix(3,3);
            else
                % Even if homography is not normalized, ensure the
                % homography has the same signal
                homographyMatrixScaleFactor = sign(homographyMatrixData(3,3));
            end
            homographyMatrix = homographyMatrix ./ homographyMatrixScaleFactor;
            
            % Create homography matrix with the matrix obtained from 
            % calibration
            self = camera.models.pinhole.Homography(homographyMatrix);
        end
        
        function [self,rmse,sse,residuals] = HomographyFromPointLineCorrespondences( worldPoints ...
                                                                                   , lines ...
                                                                                   , normalizeHomographies )
            %
            % Create homography matrix. The homography matrix is obtained 
            % from calibration using direct linear transformation from 
            % point to line correspondences.
            %
            % INPUTS:
            %   1. worldPoints - points defined in the world coordinate
            %   system. Each point should be defined in a different column.
            %   The points should not be given in homogeneous coordinates.
            %   2. lines - line parameters obtained from image plane.
            %   3. normalizeHomographies - flag to indicate if homographies
            %   should be normalized to h33 = 1. Default is false.
            %
            narginchk(2,3);
            
            if nargin <= 2
                normalizeHomographies = false;
            end
            
            % Transform to cell
            if ~iscell(worldPoints)
                worldPoints = {worldPoints};
            end
            
            % Transform world and line parameters to vector class
            worldPoints = cellfun(@(x) utils.enums.Classes.POINT().convert(x),worldPoints,'UniformOutput',false);
            lines       = utils.enums.Classes.VECTOR().convert(lines);
            
            % Since for the world points we are going to consider only
            % the (x,y) coordinates, these can be considered as pixels
            worldPixels = cellfun(@(p) image.Pixel([p.x;p.y],false),worldPoints,'UniformOutput',false);

            % Transform points to homogeneous coordinates
            worldPixels = cellfun(@(x) x.setHomogeneousCoordinates(),worldPixels,'UniformOutput',false);

            % Initialize dlt matrix
            numberPoints = [worldPixels{:}];
            numberPoints = [numberPoints.numberVectors];
            dltMatrix    = zeros(sum(numberPoints),9);
            
            % Construct direct linear transformation calibration matrix
            % considering (a,b,c) as the image line parameters and M as the
            % point in world coordinates using only (x,y) coordinates:
            %
            %   - Linear Solution
            %       A = [ a*M^T b*M^T c*M^T ]
            iRow = 0;
            for iLine = 1:lines.numberVectors
                % Obtain line and point correspondences
                line       = lines.obtainVectors(iLine);
                linePoints = worldPixels{iLine};
                
                dltMatrix(iRow + 1:iRow + linePoints.numberVectors,:) = ...
                            [ line.data(1) .* linePoints.data' ...
                            , line.data(2) .* linePoints.data' ...
                            , line.data(3) .* linePoints.data' ];
                        
                % Update row index
                iRow = iRow + linePoints.numberVectors;
            end
            
            % The homography matrix corresponds to the null space of this
            % matrix
            [~,~,eigenvectors] = svd(dltMatrix);
            nullSpaceSolution  = eigenvectors(:, end);
            
            % Obtain residuals and errors for least squares problem
            if nargout > 1
                residuals = abs(dltMatrix * nullSpaceSolution);
                sse  = sum(residuals.^2);               % Sum of squared errors
                rmse = sqrt(sse / length(residuals));   % Root mean squared error
            end                
            
            % Obtain homography matrix from null space solution. The matrix
            % is 3 x 3.
            homographyMatrix = reshape( nullSpaceSolution, 3, 3)';
                                 
            % Since homography matrix is defined up to a scale factor,
            % normalize homography matrix with entry h33 = 1.
            if normalizeHomographies == true
                homographyMatrixScaleFactor = homographyMatrix(3,3);
            else
                % Even if homography is not normalized, ensure the
                % homography has the same signal
                homographyMatrixScaleFactor = sign(homographyMatrixData(3,3));
            end
            homographyMatrix = homographyMatrix ./ homographyMatrixScaleFactor;
            
            % Create homography matrix with the matrix obtained from 
            % calibration
            self = camera.models.pinhole.Homography(homographyMatrix);
        end
    end
end
