classdef StructureTensorSize < image.ImageSize
    %STRUCTURETENSORSIZE
    %   Utility to define the structure tensor size.
    
    properties (Dependent)
        numberDimensions    % Number of dimensions of structure tensor
    end
    
    methods
        function self = StructureTensorSize(varargin)
            %
            % Structure tensor size instance.
            %
            % INPUTS:
            %   1. structureTensorSize - array with information of 
            %   structure tensor size.
            %
            narginchk(0,1)
            
            self@image.ImageSize(varargin{:});
        end
        
        function number = get.numberDimensions(self)
            number = self.numberItemsInDimension(4);
        end
    end
end