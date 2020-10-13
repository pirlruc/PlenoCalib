classdef RotationMatrix
    %ROTATIONMATRIX
    %   Utility to generate rotation matrices.
    
    properties
        data = []           % Rotation matrix data
    end
    
    properties (Dependent)
        vector              % Rotation vector associated with rotation matrix
        matrixLogarithm     % Rotation matrix logarithm. This format allows to sum
                            % and add rotation matrices. To get back a rotation 
                            % matrix use expm function
        matrixExponential   % Rotation matrix exponential
        eulerAngles         % Euler angles representation assuming ZYX order
        quaternion          % Quaternion representation
    end
    
    methods
        function self = RotationMatrix(varargin)
            %
            % Rotation matrix instance.
            %
            % INPUTS:
            %   1. rotationMatrixData - rotation matrix.
            %
            narginchk(0,1);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
        
        function self = set.data(self,newRotationMatrixData)
            self.data = newRotationMatrixData;
        end
        
        function vector = get.vector(self)
            %
            % Obtain rotation vector associated with rotation matrix using
            % the rodrigues formula.
            %
            
            % If data is empty, return empty vector
            if isempty(self.data)
                vector = [];
            else
                vector = rodrigues(self.data);
            end
        end
        
        function matrix = get.matrixLogarithm(self)
            matrix = logm(self.data);
        end
        
        function matrix = get.matrixExponential(self)
            matrix = expm(self.data);
        end
        
        function angles = get.eulerAngles(self)
            angles = rotm2eul(self.data,'ZYX');
        end
        
        function quaternion = get.quaternion(self)
            % Ensure the rotation matrix is orthonormal
            quaternion = rotm2quat(rodrigues(self.vector));
        end        
    end
    
    methods (Static)
        function self = RotationMatrixFromAngles(angle_x, angle_y, angle_z)
            %
            % Obtain the rotation matrix as a product of the matrices
            % representing each individual rotation.
            %
            % INPUTS:
            %   1. angle_x - rotation relatively to the x-axis.
            %   2. angle_y - rotation relatively to the y-axis.
            %   3. angle_z - rotation relatively to the z-axis.
            %
            
            % Obtain individual rotation matrices
            rotation_x = [1, 0, 0; 0, cos(angle_x), sin(angle_x); 0, -sin(angle_x), cos(angle_x)];
            rotation_y = [cos(angle_y), 0, -sin(angle_y); 0, 1, 0; sin(angle_y), 0, cos(angle_y)];
			rotation_z = [cos(angle_z), sin(angle_z), 0; -sin(angle_z), cos(angle_z), 0; 0, 0, 1];
            
            % Compute rotation matrix data and obtain matrix
            self = math.RotationMatrix(rotation_z * rotation_y * rotation_x);
        end
        
        function self = RotationMatrixFromRotationVector(rotationVector)
            %
            % Obtain rotation matrix from rotation vector using the
            % rodrigues formula.
            %
            self = math.RotationMatrix(rodrigues(rotationVector));
        end
        
        function [self,scaleFactor] = RotationMatrixFrom2Vectors(vector_r1,vector_r2)
            %
            % Obtain rotation matrix from 2 vectors assuming that the
            % vectors are not normalized.
            % 
            % INPUTS: 
            %   1. vector_r1 - rotation vector r1.
            %   2. vector_r2 - rotation vector r2.
            %
            narginchk(2,2);
            
            % Obtain scale factor for vector r1 and correct vector
            scaleFactor_r1    = sqrt(sum(vector_r1.^2));
            scaleFactor_r2    = sqrt(sum(vector_r2.^2));
            scaleFactor       = 0.5 * (scaleFactor_r1 + scaleFactor_r2);
            rotationMatrix_r1 = vector_r1 ./ scaleFactor;
            rotationMatrix_r2 = vector_r2 ./ scaleFactor;
            
            % Obtain rotation vector r3 and correct scaling
            rotationMatrix_r3 = cross(rotationMatrix_r1,rotationMatrix_r2);
            
            % Create rotation matrix
            self = math.RotationMatrix();
            self.data = [rotationMatrix_r1, rotationMatrix_r2, rotationMatrix_r3];
            
            % Obtain consistent rotation matrix
            self = math.RotationMatrix.RotationMatrixFromRotationVector(self.vector);
            
            % Display warning if rotation matrix is complex
            if any(~isreal(self.data))
                warning( [ 'RotationMatrix:RotationMatrixFrom2Vectors - ' ...
                         , 'Rotation matrix is complex, only real values are stored.' ]);
                self.data = real(self.data);
            end
        end
        
        function self = RotationMatrixFromEulerAngles(eulerAngles,rotationSequence)
            %
            % Obtain rotation matrix from euler angles.
            %
            % INPUTS:
            %   1. eulerAngles - euler angles.
            %   2. rotationSequence - rotation sequence for euler angles.
            %   Default sequence is 'ZYX'.
            %
            narginchk(1,2);
            
            % Default rotation sequence
            if nargin <= 1
                rotationSequence = 'ZYX';
            end
            
            % Obtain rotation matrix
            self = math.RotationMatrix( eul2rotm(eulerAngles,rotationSequence) );
        end
        
        function self = RotationMatrixFromQuaternion(quaternion)
            %
            % Obtain rotation matrix from quaternion. The quaternion
            % representation corresponds to the euler parameters.
            %
            % INPUTS:
            %   1. quaternion - quaternion data.
            %
            narginchk(1,1);
            
            % Obtain rotation matrix
            self = math.RotationMatrix( quat2rotm(quaternion) );
        end
        
        function self = RotationMatricesFromSO3LieGroup(number,mean,covariance)
            %
            % Obtain random rotation matrices by randomly combining SO3 
            % Lie group generators. The matrix is generated assuming a
            % normal distribution with a given mean and variance.
            %
            % INPUTS:
            %   1. number - number of rotation matrices to generate.
            %   2. mean   - mean considered to generate random rotation
            %   matrices.
            %   3. covariance - covariance considered to generate random
            %   rotation matrices.
            %
            narginchk(0,3);
            
            if nargin <= 0
                number = 1;
            end
            
            if nargin <= 1
                mean = zeros(3,1);
            end
            
            if nargin <= 1
                covariance = eye(3,3);
            end
            
            % Obtain multivariate normal distribution. Since this is a
            % rotation matrix, we need 3 variables one for each axis.
            normal = matlabshared.tracking.internal.NormalDistribution(3);
            normal.Mean       = mean;
            normal.Covariance = covariance;
            
            % Generate random weights for each generator
            randomWeight = math.Point(normal.sample(number),false);

            % Obtain rotation matrix by summing the generators of the SO3
            % Lie group.
            self = math.RotationMatrix( expm( [               0, -randomWeight.z,  randomWeight.y ...
                                              ;  randomWeight.z,               0, -randomWeight.x ...
                                              ; -randomWeight.y,  randomWeight.x,               0 ] ) );
        end
        
        function self = RotationMatricesFromAverage(rotationMatrices)
            %
            % Obtain average of rotation matrices.
            %
            % INPUTS: 
            %   1. rotationMatrices - rotation matrices
            %
            narginchk(1,1);
            
            % Obtain number of rotation matrices
            numberRotationMatrices = length(rotationMatrices);
            
            % Obtain quaternion list 
            quaternionList = reshape([rotationMatrices.quaternion],4,numberRotationMatrices)';

            % Obtain average quaternion and rotation matrix
            quaternion = avg_quaternion_markley(quaternionList)';
            self       = math.RotationMatrix.RotationMatrixFromQuaternion(quaternion);
        end
    end
end
