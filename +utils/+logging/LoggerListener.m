classdef LoggerListener < handle
    %LOGGERLISTENER
    %   Methods and utilities for logging in Matlab.
    
    properties (Access = protected)
        listener    % Logger listener
    end
    
    methods
        function self = LoggerListener(loggerObject)
            %
            % Create logger object and listener.
            %
            % INPUTS:
            %   1. loggerObject - logger object data.
            %
            narginchk(0,1);
            
            if nargin <= 0
                loggerObject = utils.logging.Logger.getLogger();
            end
            
            % Add listener to logger
            self.listener = addlistener(loggerObject,'logEvent',@utils.logging.LoggerListener.handleLogEvent);
        end
    end
    
    methods
        function activate(self)
            %
            % Activate listener for logging.
            %
            self.listener.Enabled = true;
        end
        
        function deactivate(self)
            %
            % Deactivate listener for logging.
            %
            self.listener.Enabled = false;
        end
        
        function delete(self)
            %
            % Delete listener for logging.
            %
            delete(self.listener);
        end
    end
    
    methods (Static)
        function handleLogEvent(loggerObject,logEventData)
           %
           % Log event data.
           %
           % INPUTS:
           %   1. loggerObject - logger object with properties for logging.
           %   2. logEventData - event data for logging.
           %
           narginchk(2,2);

           if logEventData.type.keyword < loggerObject.commandWindowLevel.keyword ...
           && logEventData.type.keyword < loggerObject.logLevel.keyword ...
           && logEventData.type.keyword < loggerObject.emailLogLevel.keyword
               return
           end
           
           % Obtain call stack
           callStack = dbstack;
           callStack = callStack(end);
           
           % Display to command window if log event type is greater than
           % the level defined for command window
           if logEventData.type.keyword >= loggerObject.commandWindowLevel.keyword
               fprintf( '%s:%s:%s:%6d - %s\n' ...
                      , logEventData.type.char ...
                      , callStack.file, callStack.name, callStack.line ...
                      , logEventData.message );
           end
           
           % Write to logging file if log event type is greater than
           % the level defined for file
           if logEventData.type.keyword >= loggerObject.logLevel.keyword
               % Append new log to log file
               try
                   fileID = fopen(loggerObject.filepath,'a');
                   fprintf( fileID,'%s,%s,%s,%s,%6d - %s\r\n' ...
                          , datestr(now,'yyyy-mm-dd HH:MM:SS.FFF') ...
                          , logEventData.type.char ...
                          , callStack.file, callStack.name, callStack.line ...
                          , logEventData.message );
                   fclose(fileID);
               catch error
                   display(error);
               end
           end
           
           % Email log event if log event type is greater than the level
           % defined for emaling log
           if ~isempty(loggerObject.email.recipients) ...
           && logEventData.type.keyword >= loggerObject.logLevel.keyword
               % Send email with log record
               try
                   loggerObject.email.subject = sprintf( '%s,%s,%s,%6d' ...
                                              , logEventData.type.char ...
                                              , callStack.file, callStack.name, callStack.line );
                   loggerObject.email.message = sprintf( fileID,'%s - %s\r\n' ...
                                              , datestr(now,'yyyy-mm-dd HH:MM:SS.FFF') ...
                                              , logEventData.message );
                   loggerObject.email.send();
               catch error
                   display(error);
               end
           end
        end
    end
end
