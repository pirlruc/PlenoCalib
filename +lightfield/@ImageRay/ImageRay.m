classdef ImageRay < math.Vector
    %IMAGERAY
    %   Utility to represent image rays.
    
    properties (Constant)
        NUMBER_COMPONENTS = 4       % Number of components in image rays
    end
    
    properties (Dependent)
        i       % i-component for image ray (microlens pixel i)
        j       % j-component for image ray (microlens pixel j)
        k       % k-component for image ray (microlens k)
        l       % l-component for image ray (microlens j)
    end
    
    methods
        function self = ImageRay(varargin)
            %
            % Image ray instance.
            %
            % INPUTS:
            %   1. imageRayData - image ray entries. To represent several
            %      image rays, use a matrix with each image ray being 
            %      represented as a new column.
            %   2. homogeneous - homogeneous flag to indicate if image ray
            %      is in homogeneous coordinates or not.
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
        
        function data = get.i(self)
            % If data is empty, return empty list
            if isempty(self.data)
                data = [];
            else
                data = self.data(1,:);
            end
        end
        
        function data = get.j(self)
            % If data is empty, return empty list
            if isempty(self.data)
                data = [];
            else
                data = self.data(2,:);
            end
        end
        
        function data = get.k(self)
            % If data is empty, return empty list
            if isempty(self.data)
                data = [];
            else
                data = self.data(3,:);
            end
        end
        
        function data = get.l(self)
            % If data is empty, return empty list
            if isempty(self.data)
                data = [];
            else
                data = self.data(4,:);
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
        function imageRays = obtainVectors(self,imageRayPositions)
            %
            % Obtain image rays from the list of image rays in current 
            % object.
            %
            % INPUTS:
            %   1. imageRayPositions - positions of the image rays in the 
            %      list of image rays.
            %
            
            % Obtain vector
            vectors = obtainVectors@math.Vector(self,imageRayPositions);
            
            % Convert to image ray structure
            imageRays      = lightfield.ImageRay();
            imageRays.homogeneous = vectors.homogeneous;
            imageRays.data        = vectors.data;
        end
    end
end

