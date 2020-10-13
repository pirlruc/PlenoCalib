classdef PinholeTypes
    %PINHOLETYPES
    %   Enumerate pinhole types that can be generated for the lytro camera.
    
    properties
        keyword         % Pinhole type keyword
        description     % Pinhole type description
        parameters      % Intrinsic parameters to estimate
    end
    
    methods
        function self = PinholeTypes(keyword, description, parameters)
            %
            % Pinhole types instance.
            %
            self.keyword     = keyword;
            self.description = description;
            self.parameters  = parameters;
        end
    end
    
    enumeration
        % Do not optimize h_tl
        VIEWPOINTS  ( 'viewpoints' , 'Ensure pinholes for viewpoint cameras by not optimizing h_tl.' ...
                    , [1 0 1 1 0;0 1 0 0 0;1 0 1 1 0;0 1 0 1 0;0 0 0 0 0] );
    end
end

