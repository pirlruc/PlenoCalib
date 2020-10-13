classdef Point < math.Vector
    %POINT
    %   Utility to represent point data.
    
    properties (Constant)
        NUMBER_COMPONENTS = 3       % Number of components in point
    end
    
    properties (Dependent)
        x       % x-component for point
        y       % y-component for point
        z       % z-component for point
    end
    
    methods
        function self = Point(varargin)
            %
            % Point instance.
            %
            % INPUTS:
            %   1. pointData - point entries. To represent several
            %      points, use a matrix with each point being represented
            %      as a new column.
            %   2. homogeneous - homogeneous flag to indicate if point is
            %      in homogeneous coordinates or not.
            %
            narginchk(0,2);
            
            % Create super class instance
            self = self@math.Vector();
            
            if ~isempty(varargin)
                if nargin >= 2
                    self.homogeneous = varargin{2};
                end
                
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
        
        function data = get.x(self)
            % If data is empty, return empty list
            if isempty(self.data)
                data = [];
            else
                data = self.data(1,:);
            end
        end
        
        function data = get.y(self)
            % If data is empty, return empty list
            if isempty(self.data)
                data = [];
            else
                data = self.data(2,:);
            end
        end
        
        function data = get.z(self)
            % If data is empty, return empty list
            if isempty(self.data)
                data = [];
            else
                data = self.data(3,:);
            end
        end
    end
    
    methods (Access = protected)
        function data = correctInputData(self,newData)
            %
            % Correct input data to have points defined in each column.
            %
            % INPUT:
            %   1. newData - input data to point.
            %
            
            % Different points should be provided in different columns
            newData = math.Vector(newData);
            if     self.homogeneous == false ...
                && newData.numberComponents ~= self.NUMBER_COMPONENTS
                data = newData.data';
            elseif self.homogeneous == true ...
                && newData.numberComponents ~= (self.NUMBER_COMPONENTS + 1)
                data = newData.data';
            else
                data = newData.data;
            end
        end
    end
    
    methods
        function points = obtainVectors(self,pointPositions)
            %
            % Obtain points from the list of points in current object.
            %
            % INPUTS:
            %   1. pointPositions - positions of the points in the list of
            %      points.
            %
            
            % Obtain vector
            vectors = obtainVectors@math.Vector(self,pointPositions);
            
            % Convert to point structure
            points      = math.Point();
            points.homogeneous = vectors.homogeneous;
            points.data        = vectors.data;
        end
        
        function self = normalize(self)
            %
            % Normalize points in order to have the points centered at the
            % origin with an average distance of sqrt(3).
            %
            
            % Set homogeneous coordinates
            self = self.setHomogeneousCoordinates();

            % Normalize points
            self.data = self.obtainNormalizationMatrix * self.data;
            
            % Remove homogeneous coordinates
            self = self.removeHomogeneousCoordinates();
        end
        
        function normalizationMatrix = obtainNormalizationMatrix(self)
            %
            % Obtain normalization matrix to have the points centered at 
            % the origin with an average distance of sqrt(3).
            %
            % Obtain mean for points in the object space
            shift_x = mean(self.x);
            shift_y = mean(self.y);
            shift_z = mean(self.z);

            % Scale points to have root mean squared distance equal to
            % sqrt(3).
            scalePoints = sqrt(3) / sqrt(mean( (self.x - shift_x).^2 ...
                                             + (self.y - shift_y).^2 ...
                                             + (self.z - shift_z).^2 ));

            % Obtain normalization matrix
            normalizationMatrix  = [ scalePoints 0 0 -scalePoints * shift_x ...
                                   ; 0 scalePoints 0 -scalePoints * shift_y ...
                                   ; 0 0 scalePoints -scalePoints * shift_z ...
                                   ; 0 0 0 1 ];
        end
    end
end

