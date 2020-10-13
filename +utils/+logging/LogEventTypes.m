classdef LogEventTypes
    %LOGEVENTTYPES
    %   Enumerate log event types available.
    
    properties
        keyword         % Log event type keyword
        description     % Log event type description
    end
    
    methods
        function self = LogEventTypes(keyword, description)
            %
            % Log event type instance.
            %
            self.keyword     = keyword;
            self.description = description;
        end
    end
    
    enumeration
        ALL   (0,'Log everything.')
        TRACE (1,'Log functions called.')
        DEBUG (2,'Log processing information.')
        INFO  (3,'Information event.')
        WARN  (4,'Warning event.')
        ERROR (5,'Expected event error.')
        FATAL (6,'Unexpected event error.')
        OFF   (7,'No logging.')
    end
end