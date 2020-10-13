classdef MathDependencies < setup.enums.dependencies.TemplateDependencies
    %MATHDEPENDENCIES
    %   Math package dependencies.
    
    methods
        function self = MathDependencies(url,filepath,packagePath,sourcePath)
            %
            % Math package dependencies instance.
            %
            self@setup.enums.dependencies.TemplateDependencies(url,filepath,packagePath,sourcePath);
        end
    end
    
    enumeration
        BOUGUET( 'http://www.vision.caltech.edu/bouguetj/calib_doc/download/toolbox_calib.zip' ...
               , 'TOOLBOX_calib_2017_06_01.zip' ...
               , 'TOOLBOX_calib' ...
               , {'TOOLBOX_calib'} );
    end
end