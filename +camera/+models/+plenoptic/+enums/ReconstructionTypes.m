classdef ReconstructionTypes
    %RECONSTRUCTIONTYPES
    %   Enumerates the reconstruction types available.
    
    properties
        description         % Reconstruction method description
        method              % Method name for reconstruction
    end
    
    methods
        function self = ReconstructionTypes(description,method)
            %
            % Reconstruction types instance.
            %
            self.description = description;
            self.method      = method;
        end
    end
    
    enumeration
        POINT_RECONSTRUCTION ( 'Point reconstruction'  , 'reconstruct' ) 
        POINT_RECONSTRUCTION_FITTING_LINES ( 'Point reconstruction correcting the observations after fitting lines to the observations' ...
                                           , 'reconstructByImposingLines' )
        LINE_PARAMETER_RECONSTRUCTION ( 'Reconstruction using line parameters from fitting lines to observations' ...
                                      , 'reconstructUsingLineParameters' )
    end
end

