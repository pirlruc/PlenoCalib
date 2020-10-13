classdef ProjectionMatrix
    %PROJECTIONMATRIX
    %   Projection matrix for pinhole camera model.
    
    properties
        intrinsicMatrix = camera.models.pinhole.IntrinsicMatrix().intrinsicMatrix
                % Intrinsic matrix that allows to project points in the camera coordinate system
        extrinsicMatrix = camera.models.ExtrinsicMatrix().extrinsicMatrix 
                % Extrinsic matrix that allows to convert points from the world to the camera 
                % coordinate system
    end
    
    properties (Dependent)
        principalAxis       % Principal axis
        projectionCenter    % Projection center associated with projection matrix
        projectionMatrix    % Projection matrix that allows to project a point in the world 
                            % coordinate system to the image plane
    end
    
    methods
        function self = ProjectionMatrix(varargin)
            %
            % Create projection matrix instance.
            %
            % INPUTS:
            %   1. intrinsicMatrix - intrinsic matrix data
            %   2. extrinsicMatrix - extrinsic matrix data
            %
            narginchk(0,2);
            
            if ~isempty(varargin)
                if nargin >= 2
                    self.extrinsicMatrix = varargin{2};
                end
                
                if nargin >= 1
                    self.intrinsicMatrix = varargin{1};
                end
            end
        end 
        
        function self = set.intrinsicMatrix(self,newIntrinsicMatrixData)
            self.intrinsicMatrix = newIntrinsicMatrixData;
        end
        
        function self = set.extrinsicMatrix(self,newExtrinsicMatrixData)
            self.extrinsicMatrix = newExtrinsicMatrixData;
        end
        
        function projection = get.projectionMatrix(self)
            projection = self.intrinsicMatrix * self.extrinsicMatrix;
            % Since projection matrix is defined up to a scale factor
            % normalize projection matrix with ||[p31 p32 p33]|| = 1.
            projection = projection ./ norm(projection(3,1:3));
        end
        
        function axis = get.principalAxis(self)
            % Obtain principal axis from projection matrix. According to
            % the formulas presented in Hartley and Zisserman.
            
            % Obtain rotation matrix from extrinsic matrix not considering
            % homogeneous coordinates
            extrinsics     = camera.models.ExtrinsicMatrix.ExtrinsicMatrixFromMatrix(self.extrinsicMatrix);
            rotationMatrix = extrinsics.rotationMatrix.data;
            
            % Obtain matrix that results from multiplying the intrinsics
            % matrix with the rotation matrix: M = K * R.
            matrix = self.intrinsicMatrix * rotationMatrix;
            
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
    
    methods (Static)
        function [matrices,rmse,sse,residuals] = ProjectionMatrixFromHomographies( homographies ...
                                                                                 , withSkew )
            %
            % Estimate intrinsic and extrinsic parameters from
            % homographies. This method obtains a collection of projection
            % matrices corresponding to the different homographies.
            %
            % INPUTS:
            %   1. homographies - collection of homographies.
            %   2. withSkew     - flag indicating if skew should be
            %   considered while computing the projection matrices. Default
            %   is false.
            %
            narginchk(1,2);
            
            if nargin <= 1
                withSkew = false;
            end
            
            if withSkew == true
                degreesFreedom = 6;
            else
                degreesFreedom = 5;
            end
            
            % Obtain observation matrix based on homographies. This allows
            % to obtain a symmetric positive definite matrix B that
            % corresponds to (K * K^T)^-1 = K^-T * K^T. The matrix B has 6
            % degrees of freedom considering skew and 5 degrees of freedom
            % not considering skew.
            iRow = 0;
            dltMatrix = zeros(2 * length(homographies),degreesFreedom);
            for iHomography = 1:length(homographies)
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
                
                if withSkew == true
                    dltMatrix(iRow + 1,6) = h11 * h22 + h12 * h21;
                    dltMatrix(iRow + 2,6) = 2 * ( h11 * h21 - h12 * h22 );
                end
                
                iRow = iRow + 2;
            end
            
            % The matrix B corresponds to the null space of the observation
            % matrix
            [~,~,eigenvectors] = svd(dltMatrix);
            nullSpaceSolution  = eigenvectors(:, end);
            
            % Obtain residuals and errors for least squares problem
            if nargout > 1
                residuals = abs(dltMatrix * nullSpaceSolution);
                sse  = sum(residuals.^2);               % Sum of squared errors
                rmse = sqrt(sse / length(residuals));   % Root mean squared error
            end                
            
            % Obtain entries from matrix B and assemble matrix
            b11 = nullSpaceSolution(1); 
            b13 = nullSpaceSolution(2);
            b22 = nullSpaceSolution(3);
            b23 = nullSpaceSolution(4);
            b33 = nullSpaceSolution(5);
            
            if withSkew == true
                b12 = nullSpaceSolution(6);
            else
                b12 = 0;
            end
            
            intrinsicLikeMatrix = [b11 b12 b13; b12 b22 b23; b13 b23 b33];
            
            % Obtain intrinsic matrix using Cholesky decomposition. The
            % Cholesky decomposition gives the inverse of the intrinsic
            % matrix since it returns the matrix C from C^T * C = A.
            % This problem has two solutions (x and -x). Choose the
            % solution that gives positive diagonal elements. B matrix is
            % positive definite.
            try
                intrinsicMatrixData = inv(chol(intrinsicLikeMatrix,'upper'));
            catch error
                if strcmp(error.identifier,'MATLAB:posdef') > 0
                    intrinsicMatrixData = inv(chol(-intrinsicLikeMatrix,'upper'));
                else
                    error.rethrow;
                end
            end
            
            % Intrinsic matrix data is defined up to a scale factor, so
            % normalize intrinsic matrix with k33 = 1.
            intrinsicMatrixData = intrinsicMatrixData ./ intrinsicMatrixData(3,3);
            
            % Estimate extrinsic parameters for each homography and obtain
            % projection matrix
            matrices = [];
            for iHomography = 1:length(homographies)
                % Obtain extrinsic parameters from intrinsic matrix and
                % homography matrix. This returns the 1st and 2nd column of
                % the rotation matrix and the translation vector
                extrinsicParametersData = mldivide( intrinsicMatrixData ...
                                                  , homographies(iHomography).homographyMatrix );
                                              
                % Isolate each extrinsic parameter and obtain scale factor
                % for correcting extrinsic parameters estimate. The scale
                % factor comes from the intrinsic matrix and homography
                % matrix being defined up to a scale factor
                rotationMatrix_r1 = extrinsicParametersData(:,1);
                rotationMatrix_r2 = extrinsicParametersData(:,2);
                translationVector = extrinsicParametersData(:,3);
                
                % Obtain rotation matrix using 2 vectors and scale
                % translation
                [rotationMatrix,scaleFactor] = math.RotationMatrix.RotationMatrixFrom2Vectors( rotationMatrix_r1 ...
                                                                                             , rotationMatrix_r2 );
                translationVector = math.Point(translationVector ./ scaleFactor,false);

                % Obtain the 2nd and 3rd column of the rotation matrix
                rotationMatrix_r3 = cross(rotationMatrix_r1,rotationMatrix_r2);
                scaleFactor_r3    = sqrt(sum(rotationMatrix_r3.^2));
                rotationMatrix_r3 = rotationMatrix_r3./ scaleFactor_r3;
                rotationMatrix_r2 = cross(rotationMatrix_r3,rotationMatrix_r1);
                
                if translationVector.z < 0
                    rotationMatrix.data = [ -rotationMatrix.data(:,1:2) ...
                                          ,  rotationMatrix.data(:,3) ];
                    translationVector   = -translationVector.data;
                end
                
                % Create extrinsic matrix
                extrinsic = camera.models.ExtrinsicMatrix();
                extrinsic.translationVector = translationVector;
                extrinsic.rotationMatrix    = rotationMatrix;
                                                
                % Obtain projection matrix
                projectionMatrix = camera.models.pinhole.ProjectionMatrix();
                projectionMatrix.intrinsicMatrix = intrinsicMatrixData;
                projectionMatrix.extrinsicMatrix = extrinsic.extrinsicMatrix;
                
                % Update projection matrices
                matrices = cat(2,matrices,projectionMatrix);
            end
        end
    end
end
