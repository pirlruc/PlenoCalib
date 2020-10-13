classdef CameraTypes
    %CAMERATYPES
    %   Enumerates the camera types available in a lightfield.
    
    properties
        keyword             % Camera type keyword
        description         % Camera type description
    end
    
    properties (Dependent)
        encoding            % Encoding indices
        decoding            % Decoding indices
    end
    
    methods
        function self = CameraTypes(keyword,description)
            %
            % Camera types instance.
            %
            self.keyword     = keyword;
            self.description = description;
        end
        
        function indices = get.encoding(~)
            indices = [3,4,1,2];
        end
        
        function indices = get.decoding(self)
            [~,indices] = sort(self.encoding);
        end
    end
    
    enumeration
        VIEWPOINT ( 'ij' , 'Viewpoint cameras.' ) 
        MICROLENS ( 'kl' , 'Microlens cameras.' )
    end
end
