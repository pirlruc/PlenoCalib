classdef ImageDependencies < setup.enums.dependencies.TemplateDependencies
    %IMAGEDEPENDENCIES
    %   Image package dependencies.
    
    methods
        function self = ImageDependencies(url,filepath,packagePath,sourcePath)
            %
            % Image package dependencies instance.
            %
            self@setup.enums.dependencies.TemplateDependencies(url,filepath,packagePath,sourcePath);
        end
    end
    
    enumeration
        BOUGUET( 'http://www.vision.caltech.edu/bouguetj/calib_doc/download/toolbox_calib.zip' ...
               , 'TOOLBOX_calib_2017_06_01.zip' ...
               , 'TOOLBOX_calib' ...
               , {'TOOLBOX_calib/TOOLBOX_calib'} );
        BOK    ( 'https://drive.google.com/uc?export=download&id=0B2553ggh3QTcRkY5OWRIaURESjA' ...
               , 'LightField_GeoCalibration_ver2.zip' ...
               , 'LightField_GeoCalibration_ver2' ...
               , {} );
    end
end