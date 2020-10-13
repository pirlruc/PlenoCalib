classdef Classes
    %CLASSES
    %   Package to define classes available.
    
    properties (SetAccess = immutable)
        class           % Class name
        constructor     % Class constructor
        optional        % Optional class names
    end
   
    methods
        function self = Classes( class, constructor, optional )
            %
            % Classes instance.
            %
            self.class       = class;
            self.constructor = constructor;
            self.optional    = optional;
        end
    end
    
    methods 
        function newData = convert(self,data)
            % If data is of specified class, do not convert
            if isa(data,self.class)
                newData = data;
            else
                % Check if current class is a class from the list of 
                % optional classes
                optionalClass = false;
                for iClass = 1:length(self.optional)
                    className = self.optional{iClass};
                    if isa(data,className)
                        optionalClass = true;
                        break
                    end
                end
                
                % If it is not an optional class, convert data to specified
                % class
                if optionalClass == false
                    newData = eval(sprintf(self.constructor,'data'));
                
                % Otherwise, obtain common data field and convert to
                % specified class
                else
                    newData = eval(sprintf(self.constructor,'data.data'));
                end
            end
        end
    end
    
    enumeration
        STANDARD_PLENOPTIC_CAMERA( 'camera.models.plenoptic.standard.Camera' ...
                                 , 'camera.models.plenoptic.standard.Camera(%s)', {} )
        STRUCTURE_TENSOR_FIELD( 'image.gradient.StructureTensorField', 'image.gradient.StructureTensorField(%s)', {} )
        STRUCTURE_TENSOR( 'image.gradient.StructureTensor', 'image.gradient.StructureTensor(%s)', {} )
        LIGHTFIELD_SIZE ( 'lightfield.LightfieldSize', 'lightfield.LightfieldSize(%s)', {} )
        ROTATION_MATRIX ( 'math.RotationMatrix' , 'math.RotationMatrix(%s)'  , {} )
        LYTRO_CAMERA    ( 'camera.lytro.Camera', 'camera.lytro.Camera(%s)', {} ) 
        LIGHTFIELD( 'lightfield.Lightfield'     , 'lightfield.Lightfield(%s)', {} )
        STATISTICS( 'statistics.Statistics'     , 'statistics.Statistics(%s)', {} )
        FILTER    ( 'image.filters.Filter'      , 'image.filters.Filter(%s)' , {} )
        POINT     ( 'math.Point'  , 'math.Point(%s,false)' , {} )
        PIXEL     ( 'image.Pixel' , 'image.Pixel(%s,false)', {} ) 
        IMAGE     ( 'image.Image' , 'image.Image(%s)' ...
                  , {'abstract.TemplateImage', 'abstract.TemplateEpipolarPlaneImage'} )
        VIDEO     ( 'camera.Video', 'camera.Video(%s)', {} )
        VECTOR    ( 'math.Vector' , 'math.Vector(%s,false)', {} ) 
        PINHOLE   ( 'camera.models.pinhole.Pinhole', 'camera.models.pinhole.Pinhole(%s)' , {} )
        HOMOGRAPHY( 'camera.models.pinhole.Homography', 'camera.models.pinhole.Homography(%s)' , {} )
        PLENOPTIC_HOMOGRAPHY( 'camera.models.plenoptic.standard.Homography' ...
                            , 'camera.models.plenoptic.standard.Homography(%s)' , {} )
        IMAGE_SIZE( 'image.ImageSize'     , 'image.ImageSize(%s)' ...
                  , { 'lightfield.image.MicrolensImageSize', 'lightfield.image.ViewpointImageSize' ...
                    , 'lightfield.epi.EpipolarPlaneImageSize', 'image.epi.EpipolarPlaneImageSize' } )
        IMAGE_RAY ( 'lightfield.ImageRay' , 'lightfield.ImageRay(%s,false)' , {} ) 
        OBJECT_RAY( 'lightfield.ObjectRay', 'lightfield.ObjectRay(%s,false)', {} ) 
    end
end