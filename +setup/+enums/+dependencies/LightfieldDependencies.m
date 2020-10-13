classdef LightfieldDependencies < setup.enums.dependencies.TemplateDependencies
    %LIGHTFIELDPENDENCIES
    %   Lightfield package dependencies.
    
    methods
        function self = LightfieldDependencies(url,filepath,packagePath,sourcePath)
            %
            % Lightfield package dependencies instance.
            %
            self@setup.enums.dependencies.TemplateDependencies(url,filepath,packagePath,sourcePath);
        end
    end
    
    enumeration
        LFTOOLBOX( '', '', 'LFToolbox', {} )
        BOK    ( 'https://drive.google.com/uc?export=download&id=0B2553ggh3QTcRkY5OWRIaURESjA' ...
               , 'LightField_GeoCalibration_ver2.zip' ...
               , 'LightField_GeoCalibration_ver2' ...
               , {} );
    end
end