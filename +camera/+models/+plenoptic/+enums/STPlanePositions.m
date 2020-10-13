classdef STPlanePositions
    %STPLANEPOSITIONS
    %   Enumerate (s,t)-plane positions that can be considered for the
    %   lytro camera calibration.
    
    properties
        keyword         % (s,t)-plane position keyword
        description     % (s,t)-plane position description
        parameters      % Intrinsic parameters to estimate
    end
    
    methods
        function self = STPlanePositions(keyword, description, parameters)
            %
            % (s,t)-plane positions instance.
            %
            self.keyword     = keyword;
            self.description = description;
            self.parameters  = parameters;
        end
    end
    
    enumeration
        % Do not optimize h_sk and h_tl
        VIEWPOINTS_CENTERS_PROJECTION  ( 'viewpoints' , '(s,t)-plane corresponds to the plane containing the viewpoints centers of projection.' ...
                                       , [1 0 0 0 1;0 1 0 0 1;1 0 1 1 1;0 1 0 1 1;0 0 0 0 0] );
    end
end

