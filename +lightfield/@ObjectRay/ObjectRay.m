classdef ObjectRay < math.Vector
    %OBJECTRAY
    %   Utility to represent object rays.
    
    properties (Constant)
        NUMBER_COMPONENTS = 4       % Number of components in object rays
    end
    
    properties (Dependent)
        s       % s-component for object ray (position s)
        t       % t-component for object ray (position t)
        u       % u-component for object ray (direction u)
        v       % v-component for object ray (direction v)
    end
    
    methods
        function self = ObjectRay(varargin)
            %
            % Object ray instance.
            %
            % INPUTS:
            %   1. objectRayData - object ray entries. To represent several
            %      object rays, use a matrix with each object ray being 
            %      represented as a new column.
            %   2. homogeneous - homogeneous flag to indicate if object ray
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
        
        function data = get.s(self)
            % If data is empty, return empty list
            if isempty(self.data)
                data = [];
            else
                data = self.data(1,:);
            end
        end
        
        function data = get.t(self)
            % If data is empty, return empty list
            if isempty(self.data)
                data = [];
            else
                data = self.data(2,:);
            end
        end
        
        function data = get.u(self)
            % If data is empty, return empty list
            if isempty(self.data)
                data = [];
            else
                data = self.data(3,:);
            end
        end
        
        function data = get.v(self)
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
        function objectRays = obtainVectors(self,objectRayPositions)
            %
            % Obtain object rays from the list of object rays in current 
            % object.
            %
            % INPUTS:
            %   1. objectRayPositions - positions of the object rays in the 
            %      list of object rays.
            %
            
            % Obtain vector
            vectors = obtainVectors@math.Vector(self,objectRayPositions);
            
            % Convert to object ray structure
            objectRays      = lightfield.ObjectRay();
            objectRays.homogeneous = vectors.homogeneous;
            objectRays.data        = vectors.data;
        end
    end
end

