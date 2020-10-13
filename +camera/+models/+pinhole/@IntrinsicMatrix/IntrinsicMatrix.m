classdef IntrinsicMatrix
    %INTRINSICMATRIX
    %   Intrinsic matrix for pinhole camera model.
    
    properties
        pixelSize      = image.Pixel([1;1])  % Pixel size to convert pixel to metric units
        principalPoint = image.Pixel([0;0])  % Principal point position to define the 
                                             % intersection of the optical axis with the 
                                             % image plane. This point is defined in pixels
        skewFactor     = 0                   % Skew factor
        imagePlaneDistance = 1      % Distance between the image plane and the pinhole position. Normally,
                                    % corresponds to the focal length
    end
    
    properties (Dependent)
        pixelToMetricMatrix     % Matrix that allows to convert pixels to metric units
        perspectiveMatrix       % Perspective projection matrix
        intrinsicMatrix         % Intrinsic matrix that allows to project points in the camera
                                % coordinate system to the image plane
    end
    
    methods
        function self = IntrinsicMatrix(varargin)
            %
            % Create intrinsic matrix instance.
            %
            % INPUTS:
            %   1. pixelSize - pixel size to convert pixel to metric units.
            %   2. principalPoint     - principal point position.
            %   3. imagePlaneDistance - distance between the image plane 
            %      and the pinhole position.
            %   4. skewFactor - skew factor value
            %
            narginchk(0,4);
            
            if ~isempty(varargin)
                if nargin >= 4
                    self.skewFactor = varargin{4};
                end
                
                if nargin >= 3
                    self.imagePlaneDistance = varargin{3};
                end
                
                if nargin >= 2
                    self.principalPoint = varargin{2};
                end
                
                if nargin >= 1
                    self.pixelSize = varargin{1};
                end
            end
        end
        
        function self = set.pixelSize(self,newPixelSize)
            self.pixelSize = utils.enums.Classes.PIXEL().convert(newPixelSize);
        end
        
        function self = set.principalPoint(self,newPrincipalPoint)
            self.principalPoint = utils.enums.Classes.PIXEL().convert(newPrincipalPoint);
        end
        
        function self = set.imagePlaneDistance(self,newImagePlaneDistance)
            self.imagePlaneDistance = newImagePlaneDistance;
        end
        
        function pixelToMetric = get.pixelToMetricMatrix(self)
            pixelToMetric = [ self.pixelSize.u,  self.skewFactor, self.principalPoint.u ...
                            ;                0, self.pixelSize.v, self.principalPoint.v ...
                            ;                0,                0,                     1 ];
        end
        
        function perspective = get.perspectiveMatrix(self)
            perspective = [ self.imagePlaneDistance,                       0, 0 ...
                          ;                       0, self.imagePlaneDistance, 0 ...
                          ;                       0,                       0, 1 ];
        end
        
        function intrinsic = get.intrinsicMatrix(self)
            intrinsic = self.pixelToMetricMatrix * self.perspectiveMatrix;
        end
    end
    
    methods
        function pixels = project(self,cameraPoints,pixelsToNearestInteger)
            %
            % Obtain projection of points in the camera coordinate system.
            %
            % INPUTS:
            %   1. cameraPoints - points in the camera coordinate system.
            %      Each point should be defined in a different column. The 
            %      points should not be given in homogeneous coordinates.
            %   2. pixelsToNearestInteger - flag that indicates if pixels
            %      must be rounded to the nearest integer. The default is
            %      true.
            %
            narginchk(1,3);
            
            % If no points are provided, generate random points
            if nargin <= 1
                cameraPoints = rand(3,50);
            end
            
            % If round flag is not provided, assume that the pixels should
            % be returned rounded to the nearest integer.
            if nargin <= 2
                pixelsToNearestInteger = true;
            end
            
            % If points are provided in rows instead of columns, transpose
            % the points
            cameraPoints = utils.enums.Classes.POINT().convert(cameraPoints);
            
            % Obtain projected points
            pixels = image.Pixel(self.intrinsicMatrix * cameraPoints.data,true);
            pixels = pixels.removeHomogeneousCoordinates();
            
            % Round pixels to nearest integer
            if pixelsToNearestInteger == true
                pixels.data = round(pixels.data);
            end
        end
    end
end

