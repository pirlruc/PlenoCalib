classdef Vector
    %VECTOR
    %   Utility to represent vector data.
    
    properties
        data = []               % Vector data
        homogeneous = false     % Flag indicating if vector is in homogeneous coordinates
    end
    
    properties (Dependent)
        numberComponents    % Number of components that compose the vector
        numberVectors       % Number of vectors. The class suports a matrix representation for multiple vectors
        norm                % L2-norm of each vector.
    end
    
    methods
        function self = Vector(varargin)
            %
            % Vector instance.
            %
            % INPUTS:
            %   1. vectorData  - vector entries. To represent several
            %      vectors, use a matrix with each vector being represented
            %      as a new column.
            %   2. homogeneous - homogeneous flag to indicate if vector is
            %      in homogeneous coordinates or not.
            %
            narginchk(0,2);
            
            if ~isempty(varargin)
                if nargin >= 2
                    self.homogeneous = varargin{2};
                end
                   
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
        
        function self = set.data(self,newVectorData)
            self.data = self.correctInputData(newVectorData);
        end
        
        function self = set.homogeneous(self,newHomogeneousFlag)
            self.homogeneous = newHomogeneousFlag;
        end
        
        function number = get.numberComponents(self)
            number = size(self.data,1);
        end
        
        function number = get.numberVectors(self)
            number = size(self.data,2);
        end
        
        function value = get.norm(self)
            if isempty(self.data)
                value = math.Vector([],false);
            else
                value = math.Vector( sqrt( sum(self.data.^2, 1) ), false );
            end
        end
    end
    
    methods (Access = protected)
        function data = correctInputData(~,newData)
            %
            % This method should be used to correct the input data if
            % necessary.
            %
            % INPUT:
            %   1. newData - input data to vector.
            %
            data = newData;
        end
    end
    
    methods
        function vectors = obtainVectors(self,vectorPositions)
            %
            % Obtain vectors from the list of vectors in current object.
            %
            % INPUTS:
            %   1. vectorPositions - positions of the vectors in the list 
            %      of vectors.
            %
            if isempty(self.data)
                vectors = math.Vector();
            else
                vectors = math.Vector(self.data(:,vectorPositions),self.homogeneous);
            end
        end
        
        function vectors = obtainComponents(self,componentPositions)
            %
            % Obtain components from the list of components in current 
            % object.
            %
            % INPUTS:
            %   1. componentPositions - positions of the components in the 
            %      list of components.
            %
            if isempty(self.data)
                vectors = math.Vector();
            else
                vectors = math.Vector(self.data(componentPositions,:),self.homogeneous);
            end
        end
        
        function self = setHomogeneousCoordinates(self)
            %
            % Represent the vector using homogenous coordinates. The vector
            % is obtained by adding a new row with entry equal to one.
            %
            
            % If number of vectors or components is zero, do not set 
            % homogenous coordinates
            if isempty(self.data)
                error = MException( 'vector:setHomogeneousCoordinates:noData' ...
                                  , 'No vector data provided...' );
                error.throw();
            end
            
            % If vector is already in homogeneous coordinates, do not set
            if self.homogeneous == true
                warning( 'vector:setHomogeneousCoordinates:alreadyHomogeneous' ...
                       , 'Vector is already in homogeneous coordinates...' );
            else
                self.homogeneous = true;
                self.data = [self.data;ones(1,self.numberVectors)];
            end
        end
        
        function self = removeHomogeneousCoordinates(self)
            %
            % Remove homogeneous coordinates representation. The new vector
            % is obtained by diving each of the components with the
            % corresponding entry of the last row.
            %

            % If number of vectors or components is zero, do not remove 
            % homogenous coordinates
            if isempty(self.data)
                error = MException( 'vector:removeHomogeneousCoordinates:noData' ...
                                  , 'No vector data provided...' );
                error.throw();
            end
            
            % If number of components is less than one, do not remove 
            % homogenous coordinates
            if self.numberComponents <= 1
                error = MException( 'vector:removeHomogeneousCoordinates:minimumComponents' ...
                                  , 'Vector data must have a minimum of 2 components...' );
                error.throw();
            end
            
            if self.homogeneous == false
                warning( 'vector:removeHomogeneousCoordinates:alreadyNotHomogeneous' ...
                       , 'Vector is not in homogeneous coordinates...' );
            else
                self.homogeneous = false;
                self.data =    self.data(1:self.numberComponents - 1,:) ...
                            ./ repmat(self.data(end,:),self.numberComponents - 1,1);
            end
        end
    end
end

