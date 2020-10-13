classdef (Abstract) TemplateSize
    %TEMPLATESIZE
    %   Size template.
    
    properties 
        data = []       % Array with information about the data dimensions
    end
    
    methods
        function self = TemplateSize(varargin)
            %
            % Template size instance.
            %
            % INPUTS:
            %   1. dataSize - array with the dimension of the data.
            %
            narginchk(0,1);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
        
        function self = set.data(self, newSizeData)
            self.data = newSizeData;
        end
    end
    
    methods (Access = protected)
        function number = numberItemsInDimension(self,dimension)
            %
            % Obtain number of items in a given dimension.
            %
            % INPUTS:
            %   1. dimension - dimension of the array to obtain the number
            %   of items.
            %
            
            % If data is empty, the number of items is zero for all
            % dimensions.
            if isempty(self.data)
                number = 0;
                
            % If the dimension is greater than the number of dimensions, 
            % the number of items is 1.
            elseif dimension > length(self.data)
                number = 1;
            
            % Otherwise, return the number of items in dimension.
            else
                number = self.data(dimension);
            end
        end
    end
end
