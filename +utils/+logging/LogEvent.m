classdef LogEvent < event.EventData
    %LOGEVENT
    %   Log event data.
    
    properties
        type    = utils.logging.LogEventTypes.INFO  % Log event type
        message = ''                                % Message to log
    end
    
    methods
        function self = LogEvent(varargin)
            %
            % Create log event.
            %
            % INPUTS:
            %   1. logEventType - log event type.
            %   2. logMessage   - log message.
            % 
            narginchk(0,2);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.type = varargin{1};
                end
                
                if nargin >= 2
                    self.message = varargin{2};
                end
            end
        end
    end
end

