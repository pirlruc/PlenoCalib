classdef FilteringAlgorithms
    %FILTERINGALGORITHMS
    %   Filtering algorithms available.
    
    properties (~Dependent, SetAccess = immutable, GetAccess = public)
        algorithm       % Filtering algorithm
        description     % Description of filtering algorithms
    end
    
    methods
        function self = FilteringAlgorithms(algorithm, description)
            %
            % Filtering algorithms enumeration instance.
            %
            self.algorithm   = algorithm;
            self.description = description;
        end
    end
    
    enumeration
        CORRELATION ('corr' , 'Multidimensional filtering using correlation.');
        CONVOLUTION ('conv' , 'Multidimensional filtering using convolution.');
    end
end