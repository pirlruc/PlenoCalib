classdef IntrinsicMatrix
    %INTRINSICMATRIX
    %   Intrinsic matrix for standard plenoptic camera.
    %
    %   Remember that the indices (i,j,k,l) are one-based indices.
    
    properties
        numberPixels     = 10                       % Number of pixels per microlens image
        pixelOffset      = 6                        % Pixel offset
        spatialFrequency = 716790                   % Spatial frequency
        spatialOffset    = 1645.3                   % Spatial frequency offset
        directionalFrequency = 71950                % Directional frequency
        directionalOffset    = 164.7                % Directional frequency offset
        microlensFocalLength = 0.025e-3             % Microlenses focal length
        distanceMicrolensToMainLens = 6.6506e-3     % Distance from the microlens plane to the main lens plane
        mainLensFocalLength  = 6.45e-3              % Main lens focal length
        worldPlaneDistance   = 0.3                  % Distance between a plane in the world and the main lens position 
        distancePlanes_st_uv = 1                    % Distance between planes (s,t) and (u,v)
        sensorRotation       = 0                    % Rotation between sensor and microlens plane
        sensorTranslation    = image.Pixel([0;0])   % Translation between sensor and microlens plane
        samplingMatrix       = eye(5,5)             % Sampling matrix to consider hexagonal or rectangular sampling
    end
    
    properties (Dependent)
        relativeToAbsoluteMatrix    % Convert the local reference frame of the microimages to a global 
                                    % reference frame
        misalignmentMatrix          % Realign microlenses and sensor planes
        samplesToMetricMatrix       % Convert the pixel and microlenses samples to metric units
        metricToDirectionMatrix     % Convert the two-plane parameterization from 2 points for a point
                                    % and a direction
        imageToMicrolensMatrix      % Propagate the point in the ray from the image to the microlens plane
        microlensToMainLensMatrix   % Propagate the point in the ray from the microlens to the main lens plane
        mainLensRefractionMatrix    % Modifies the direction of the ray due to refraction on the main lens
        mainLensToWorldMatrix       % Propagate the point in the ray from the main lens to the world plane
        worldTo2PointsMatrix        % Propagate the direction to a point in the (u,v) plane
        
        imagePlaneMatrix            % Back-project image rays to the image plane in metric units
        microlensPlaneMatrix        % Back-project image rays to the microlens plane in metric units
        mainLensPlaneMatrix         % Back-project image rays to the main lens plane in metric units
        intrinsicMatrixDirection             
            % Back-project image rays to the world plane in metric units.
            % Applying this intrinsic matrix generates a two-plane 
            % parameterization with one point and one direction.
        intrinsicMatrixPoints
            % Back-project image rays to the world plane in metric units.
            % Applying this intrinsic matrix generates a two-plane 
            % parameterization with two points.
    end
    
    methods
        function self = IntrinsicMatrix(varargin)
            %
            % Create intrinsic matrix instance.
            %
            % INPUTS:
            %   01. numberPixels     - number of pixels per microlens
            %   image.
            %   02. pixelOffset      - pixel offset.
            %   03. spatialFrequency - spatial frequency.
            %   04. spatialOffset    - spatial frequency offset.
            %   05. directionalFrequency - directional frequency.
            %   06. directionalOffset    - directional frequency offset.
            %   07. microlensFocalLength - microlenses focal length.
            %   08. distanceMicrolensToMainLens - distance from the 
            %   microlens plane to the main lens plane.
            %   09. mainLensFocalLength  - main lens focal length.
            %   10. worldPlaneDistance   - distance between a plane in the 
            %   world and the main lens position.
            %   11. distancePlanes_st_uv - distance between planes (s,t) 
            %   and (u,v)
            %   12. sensorRotation    - angle of misalignment between
            %   microlens plane and sensor.
            %   13. sensorTranslation - translation between microlens plane
            %   and sensor.
            %   14. hexagonalSampling - flag to indicate if hexagonal
            %   sampling should be considered. Default is false.
            %
            narginchk(0,14);
            
            if nargin <= 13
                hexagonalSampling = false;
            else
                hexagonalSampling = varargin{14};
            end

            if hexagonalSampling == true
                self.samplingMatrix = [ 1 0       0         0 0 ...
                                      ; 0 1       0         0 0 ...
                                      ; 0 0 sqrt(3) sqrt(3)/2 0 ...
                                      ; 0 0       0       3/2 0 ...
                                      ; 0 0       0         0 1 ];
            else
                self.samplingMatrix = eye(5,5);
            end
            
            if ~isempty(varargin)
                if nargin >= 13
                    self.sensorTranslation = varargin{13};
                end
                
                if nargin >= 12
                    self.sensorRotation = varargin{12};
                end
                
                if nargin >= 11
                    self.distancePlanes_st_uv = varargin{11};
                end
                
                if nargin >= 10
                    self.worldPlaneDistance = varargin{10};
                end
                
                if nargin >= 9
                    self.mainLensFocalLength = varargin{9};
                end
                
                if nargin >= 8
                    self.distanceMicrolensToMainLens = varargin{8};
                end
                
                if nargin >= 7
                    self.microlensFocalLength = varargin{7};
                end
                
                if nargin >= 6
                    self.directionalOffset = varargin{6};
                end
                
                if nargin >= 5
                    self.directionalFrequency = varargin{5};
                end
                
                if nargin >= 4
                    self.spatialOffset = varargin{4};
                end
                
                if nargin >= 3
                    self.spatialFrequency = varargin{3};
                end
                
                if nargin >= 2
                    self.pixelOffset = varargin{2};
                end
                
                if nargin >= 1
                    self.numberPixels = varargin{1};
                end
            end
        end
        
        function self = set.numberPixels(self,newNumberPixels)
            self.numberPixels = newNumberPixels;
        end
        
        function self = set.pixelOffset(self,newPixelOffset)
            self.pixelOffset = newPixelOffset;
        end
        
        function self = set.spatialFrequency(self,newSpatialFrequency)
            self.spatialFrequency = newSpatialFrequency;
        end
        
        function self = set.spatialOffset(self,newSpatialOffset)
            self.spatialOffset = newSpatialOffset;
        end
        
        function self = set.directionalFrequency(self,newDirectionalFrequency)
            self.directionalFrequency = newDirectionalFrequency;
        end
        
        function self = set.directionalOffset(self,newDirectionalOffset)
            self.directionalOffset = newDirectionalOffset;
        end
        
        function self = set.microlensFocalLength(self,newMicrolensFocalLength)
            self.microlensFocalLength = newMicrolensFocalLength;
        end
        
        function self = set.distanceMicrolensToMainLens(self,newDistanceMicrolensToMainLens)
            self.distanceMicrolensToMainLens = newDistanceMicrolensToMainLens;
        end
        
        function self = set.mainLensFocalLength(self,newMainLensFocalLength)
            self.mainLensFocalLength = newMainLensFocalLength;
        end
        
        function self = set.worldPlaneDistance(self,newWorldPlaneDistance)
            self.worldPlaneDistance = newWorldPlaneDistance;
        end
        
        function self = set.distancePlanes_st_uv(self,newDistancePlanes_st_uv)
            self.distancePlanes_st_uv = newDistancePlanes_st_uv;
        end
        
        function self = set.sensorTranslation(self,newTranslation)
            self.sensorTranslation = utils.enums.Classes.PIXEL().convert(newTranslation);
        end
        
        function relativeToAbsolute = get.relativeToAbsoluteMatrix(self)
            relativeToAbsolute = [ 1, 0, self.numberPixels,                 0, -self.pixelOffset ...
                                 ; 0, 1,                 0, self.numberPixels, -self.pixelOffset ...
                                 ; 0, 0,                 1,                 0,                 0 ...
                                 ; 0, 0,                 0,                 1,                 0 ...
                                 ; 0, 0,                 0,                 0,                 1 ];
        end
        
        function misalignment = get.misalignmentMatrix(self)
            misalignment = [  cos(self.sensorRotation), sin(self.sensorRotation), 0, 0, -self.sensorTranslation.u ...
                           ; -sin(self.sensorRotation), cos(self.sensorRotation), 0, 0, -self.sensorTranslation.v ...
                           ;                         0,                        0, 1, 0,                       0 ...
                           ;                         0,                        0, 0, 1,                       0 ...
                           ;                         0,                        0, 0, 0,                       1 ];
        end
        
        function samplesToMetric = get.samplesToMetricMatrix(self)
            samplesToMetric    = [ 1/self.spatialFrequency,                       0,                           0,                           0,         -self.spatialOffset/self.spatialFrequency ...
                                 ;                       0, 1/self.spatialFrequency,                           0,                           0,         -self.spatialOffset/self.spatialFrequency ...
                                 ;                       0,                       0, 1/self.directionalFrequency,                           0, -self.directionalOffset/self.directionalFrequency ...
                                 ;                       0,                       0,                           0, 1/self.directionalFrequency, -self.directionalOffset/self.directionalFrequency ...
                                 ;                       0,                       0,                           0,                           0,                                                 1 ];
        end

        function metricToDirection = get.metricToDirectionMatrix(self)
            metricToDirection = [                            1,                            0,                           0,                           0, 0 ...
                                ;                            0,                            1,                           0,                           0, 0 ...
                                ; -1/self.microlensFocalLength,                            0, 1/self.microlensFocalLength,                           0, 0 ...
                                ;                            0, -1/self.microlensFocalLength,                           0, 1/self.microlensFocalLength, 0 ...
                                ;                            0,                            0,                           0,                           0, 1 ];
        end
        
        function imageToMicrolens = get.imageToMicrolensMatrix(self)
            imageToMicrolens = [ 1, 0, self.microlensFocalLength,                         0, 0 ...
                               ; 0, 1,                         0, self.microlensFocalLength, 0 ...
                               ; 0, 0,                         1,                         0, 0 ...
                               ; 0, 0,                         0,                         1, 0 ...
                               ; 0, 0,                         0,                         0, 1 ];
        end
            
        function microlensToMainLens = get.microlensToMainLensMatrix(self)
            microlensToMainLens = [ 1, 0, self.distanceMicrolensToMainLens,                                0, 0 ...
                                  ; 0, 1,                                0, self.distanceMicrolensToMainLens, 0 ...
                                  ; 0, 0,                                1,                                0, 0 ...
                                  ; 0, 0,                                0,                                1, 0 ...
                                  ; 0, 0,                                0,                                0, 1 ];
        end

        function mainLensRefraction = get.mainLensRefractionMatrix(self)
            mainLensRefraction = [                           1,                           0, 0, 0, 0 ...
                                 ;                           0,                           1, 0, 0, 0 ...
                                 ; -1/self.mainLensFocalLength,                           0, 1, 0, 0 ...
                                 ;                           0, -1/self.mainLensFocalLength, 0, 1, 0 ...
                                 ;                           0,                           0, 0, 0, 1 ];
        end
        
        function mainLensToWorld = get.mainLensToWorldMatrix(self)
            mainLensToWorld = [ 1, 0, self.worldPlaneDistance,                       0, 0 ...
                              ; 0, 1,                       0, self.worldPlaneDistance, 0 ...
                              ; 0, 0,                       1,                       0, 0 ...
                              ; 0, 0,                       0,                       1, 0 ...
                              ; 0, 0,                       0,                       0, 1 ];
        end
        
        function worldTo2Points = get.worldTo2PointsMatrix(self)
            worldTo2Points = [ 1, 0,                         0,                         0, 0 ...
                             ; 0, 1,                         0,                         0, 0 ...
                             ; 1, 0, self.distancePlanes_st_uv,                         0, 0 ...
                             ; 0, 1,                         0, self.distancePlanes_st_uv, 0 ...
                             ; 0, 0,                         0,                         0, 1 ];
        end
        
        function imagePlane = get.imagePlaneMatrix(self)
            imagePlane = self.samplesToMetricMatrix * self.misalignmentMatrix * self.relativeToAbsoluteMatrix * self.samplingMatrix;
        end

        function microlensPlane = get.microlensPlaneMatrix(self)
            microlensPlane = self.imageToMicrolensMatrix * self.metricToDirectionMatrix * self.imagePlaneMatrix;
        end

        function mainLensPlane = get.mainLensPlaneMatrix(self)
            mainLensPlane = self.microlensToMainLensMatrix * self.microlensPlaneMatrix;
        end

        function intrinsic = get.intrinsicMatrixDirection(self)
            intrinsic = self.mainLensToWorldMatrix * self.mainLensRefractionMatrix * self.mainLensPlaneMatrix;
        end
        
        function intrinsic = get.intrinsicMatrixPoints(self)
            intrinsic = self.worldTo2PointsMatrix * self.mainLensRefractionMatrix * self.mainLensPlaneMatrix;
        end
    end    
end
