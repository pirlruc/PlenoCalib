classdef ExtrinsicMatrix
    %EXTRINSICMATRIX
    %   Extrinsic matrix that allows to switch the world to the camera
    %   coordinate system.
    
    properties
        rotationMatrix    = math.RotationMatrix(eye(3,3))
                                % Rotation matrix that represent the relative orientations 
                                % between the two coordinate systems
        translationVector = math.Point(zeros(3,1),false)
                                % Translation vector that define the relative position of
                                % the origin of each coordinate system
    end
    
    properties (Dependent)
        extrinsicMatrix         % Extrinsic matrix that allows to convert points in the
                                % world coordinate system to the camera coordinate system
    end
    
    methods
        function self = ExtrinsicMatrix(varargin)
            %
            % Create extrinsic matrix instance.
            %
            % INPUTS:
            %   1. rotationMatrix    - rotation matrix to transform world
            %      to camera coordinate system.
            %   2. translationVector - translation vector to transform 
            %      world to camera coordinate system.
            %
            narginchk(0,2);
            
            if ~isempty(varargin)
                if nargin >= 2
                    self.translationVector = varargin{2};
                end
            
                if nargin >= 1
                    self.rotationMatrix = varargin{1};
                end
            end
        end
        
        function self = set.rotationMatrix(self,newRotationMatrixData)
            self.rotationMatrix = utils.enums.Classes.ROTATION_MATRIX().convert(newRotationMatrixData);
        end
        
        function self = set.translationVector(self,newTranslationVectorData)
            % Define translation vector as a point
            self.translationVector = utils.enums.Classes.POINT().convert(newTranslationVectorData);
        end
        
        function extrinsic = get.extrinsicMatrix(self)
            % The extrinsic matrix is provided in homogeneous coordinates
            extrinsic = [ self.rotationMatrix.data, self.translationVector.data ];
        end
    end
    
    methods
        function parameters = encode(self)
            %
            % Encode extrinsic parameters into 6 parameters.
            %
            
            % Obtain translation parameters 
            translationParameters = self.translationVector.data';
            
            % Obtain rotation vector for rotation matrix
            rotationParameters = self.rotationMatrix.vector';
            
            parameters = [translationParameters,rotationParameters];
        end
        
        function cameraPoints = obtainCameraPoints(self,worldPoints)
            %
            % Convert points in the world coordinate system to the camera
            % coordinate system.
            %
            % INPUTS:
            %   1. worldPoints - points in the world coordinate system. 
            %      Each point should be defined in a different column. The 
            %      points should not be given in homogeneous coordinates.
            %
            narginchk(1,2);
            
            % If points are not provided, generate random points
            if nargin <= 1
                worldPoints = rand(3,50);
            end
            
            % If points are provided in rows instead of columns, transpose
            % the points
            worldPoints = utils.enums.Classes.POINT().convert(worldPoints);
            
            % Obtain the point in homogeneous coordinates
            worldPoints = worldPoints.setHomogeneousCoordinates();
            
            cameraPoints = math.Point(self.extrinsicMatrix * worldPoints.data,false);
        end
        
        function worldPoints = obtainWorldPoints(self,cameraPoints)
            %
            % Convert points in the camera coordinate system to the world
            % coordinate system.
            %
            % INPUTS:
            %   1. cameraPoints - points in the camera coordinate system. 
            %      Each point should be defined in a different column. The 
            %      points should not be given in homogeneous coordinates.
            %
            narginchk(1,2);
            
            % If points are not provided, generate random points
            if nargin <= 1
                cameraPoints = rand(3,50);
            end
            
            % If points are provided in rows instead of columns, transpose
            % the points
            cameraPoints = utils.enums.Classes.POINT().convert(cameraPoints);
            
            % Obtain the point in homogeneous coordinates
            cameraPoints = cameraPoints.setHomogeneousCoordinates();
            
            % Obtain extrinsic matrix for converting camera to world
            % coordinates
            cameraToWorldMatrix = [ self.rotationMatrix.data' ...
                                  , -self.rotationMatrix.data' * self.translationVector.data ];
            
            worldPoints = math.Point(cameraToWorldMatrix * cameraPoints.data,false);
        end
    end
    
    methods (Static)
        function self = ExtrinsicMatrixFromMatrix(extrinsicMatrix)
            %
            % Obtain components of extrinsic matrix data and create
            % extrinsic matrix object.
            %
            % INPUTS:
            %   1. extrinsicMatrix - extrinsic matrix data.
            % 
            narginchk(1,1);
            
            % Initialize extrinsic matrix data
            self = camera.models.ExtrinsicMatrix();
            self.rotationMatrix    = extrinsicMatrix(:,1:3);
            self.translationVector = extrinsicMatrix(:,end);
        end
        
        function self = ExtrinsicMatrixFromEncodedParameters(parameters)
            %
            % Obtain extrinsic matrix from rotation and translation
            % parameters.
            %
            % INPUTS: 
            %   1. parameters - extrinsic parameters in vectorized form (6
            %   parameters - 3 for rotation + 3 for translation).
            %
            narginchk(1,1);
            
            % Obtain rotation and translation vectors
            translationParameters = parameters(1:3);
            rotationParameters    = parameters(4:6);
            
            % Create extrinsic matrix object
            self = camera.models.ExtrinsicMatrix();
            self.rotationMatrix    = math.RotationMatrix.RotationMatrixFromRotationVector(rotationParameters);
            self.translationVector = translationParameters;
        end
        
        function [self,procrustesError,meanDistance] = ExtrinsicMatrixUsingProcrustes( worldPoints ...
                                                                                     , cameraPoints ...
                                                                                     , showResults )
            %
            % Obtain extrinsic matrix to convert between world and camera
            % coordinate system.
            %
            % INPUTS:
            %   1. worldPoints  - points in the world coordinate system.
            %   2. cameraPoints - points in the camera coordinate system.
            %   3. showResults  - plot transformation results. Default is
            %   false.
            %
            narginchk(0,3);
            
            if nargin <= 0
                worldPoints  = randn(3,100);
            end
            
            if nargin <= 1
                cameraPoints = randn(3,100);
            end
            
            if nargin <= 2
                showResults = false;
            end
        
            % If points are provided in rows instead of columns, transpose
            % the points
            worldPoints  = utils.enums.Classes.POINT().convert(worldPoints);
            cameraPoints = utils.enums.Classes.POINT().convert(cameraPoints);
            
            % Obtain transformation from camera to world coordinates
            % T - gives the rotation
            % c - gives the translation
            [procrustesError,worldPointsAfterProcrustes,extrinsics] = procrustes( worldPoints.data' ...
                                                                                , cameraPoints.data' ...
                                                                                , 'scaling',false );

            % Initialize extrinsic matrix and obtain transformation from
            % world to camera coordinates:
            %       m_c = R m_w - R c
            %
            self = camera.models.ExtrinsicMatrix();
            self.rotationMatrix    = extrinsics.T;
            self.translationVector = -extrinsics.T * extrinsics.c(1,:)';
            
            % Obtain mean distance
            if nargout > 2
                distances    = math.Point(worldPoints.data' - worldPointsAfterProcrustes,false);
                meanDistance = sum(distances.norm.data,2) ./ distances.numberVectors;
            end

            if showResults == true
                % Transform to point instance
                worldPointsAfterProcrustes = math.Point(worldPointsAfterProcrustes,false);
            
                % Show results of procrustes estimation
                figure();
                hold on
                scatter3(worldPoints.x,worldPoints.z,worldPoints.y,'b.')
                scatter3(worldPointsAfterProcrustes.x,worldPointsAfterProcrustes.z,worldPointsAfterProcrustes.y,'ro')
                grid on
                xlabel('x'); ylabel('z'); zlabel('y');
                view([90,0]);
                hold off
            end
        end
    end
end

