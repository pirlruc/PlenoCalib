classdef ProjectionTypes
    %PROJECTIONTYPES
    %   Enumerate projection types that can be considered for lytro camera.
    
    properties
        keyword         % Projection type keyword
        description     % Projection type description
    end
    
    methods
        function self = ProjectionTypes(keyword, description)
            %
            % Projection types instance.
            %
            self.keyword     = keyword;
            self.description = description;
        end
    end
    
    enumeration
        MICROLENSES ( 'microlenses', 'Projection considering that microlens coordinates (k,l) are integers.' );
        VIEWPOINTS  ( 'viewpoints' , 'Projection considering that viewpoint coordinates (i,j) are integers.' );
        MIXTURE     ( 'mixture'    , 'Projection considering either viewpoint or microlens coordinates as integers.' );
    end
end

