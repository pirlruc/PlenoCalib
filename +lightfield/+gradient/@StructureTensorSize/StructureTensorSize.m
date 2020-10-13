classdef StructureTensorSize < lightfield.LightfieldSize
    %STRUCTURETENSORSIZE
    %   Utility to define the structure tensor size.
    
    properties (Dependent)
        numberEntries    % Number of independent entries of tensor
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
            
            self@lightfield.LightfieldSize(varargin{:});
        end
        
        function number = get.numberEntries(self)
            number = self.numberItemsInDimension(6);
        end
    end
end