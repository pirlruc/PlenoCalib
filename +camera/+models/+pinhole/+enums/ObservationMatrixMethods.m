classdef ObservationMatrixMethods
    %OBSERVATIONMATRIXMETHODS
    %   Enumerates the observation matrix methods available.
    
    properties
        description         % Observation matrix method description
        method              % Method name for observation matrix
    end
    
    methods
        function self = ObservationMatrixMethods(description,method)
            %
            % Observation matrix method instance.
            %
            self.description = description;
            self.method      = method;
        end
    end
    
    enumeration
        LINEAR_SYSTEM( 'Observation matrix constructed using the solution of the projection equation.'  , 'linear' ) 
        CROSS_PRODUCT( 'Observation matrix constructed applying the cross-product to the projection equation.', 'cross' )
    end
end

