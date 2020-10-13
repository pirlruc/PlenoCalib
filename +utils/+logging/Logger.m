classdef Logger < handle
    %LOGGER
    %   Methods and utilities for logging in Matlab.
    %
    % The code is based on log4m object (Luke Winslow).
    %
    
    properties
        filepath = 'logging.log'                                % Filepath for logging
        email    = utils.Email()                                % Email for sending log
        logLevel = utils.logging.LogEventTypes.INFO             % File level logging
        commandWindowLevel = utils.logging.LogEventTypes.ALL    % Command window level logging
        emailLogLevel      = utils.logging.LogEventTypes.ERROR  % Email level logging
    end
    
    events
        logEvent    % Event for logging
    end
    
    methods
        function self = Logger(varargin)
            %
            % Create logger object.
            %
            % INPUTS:
            %   1. logFilepath - filepath for logging file.
            %   2. commandWindowLevel - log level for command window.
            %   3. logLevel        - log level for logging file.
            %   4. emailLogLevel   - log level for emailing message.
            %   5. emailFilepath   - email configuration filepath.
            %   6. emailRecipients - email recipients.
            %
            narginchk(0,6);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.filepath = varargin{1};
                end
                
                if nargin >= 2
                    self.commandWindowLevel = varargin{2};
                end
                
                if nargin >= 3
                    self.logLevel = varargin{3};
                end
                
                if nargin >= 4
                    self.emailLogLevel = varargin{4};
                end
                
                if nargin >= 5
                    self.email = varargin{5};
                end
                
                if nargin >= 6
                    self.email.recipients = varargin{6};
                end
            end
        end
        
        function set.filepath(self,newFilepath)
            [fileID,message] = fopen(newFilepath, 'a');
            if(fileID < 0)
                error = MException( 'Logger:setFilepath' ...
                                  , ['Problem with supplied filepath: ' message] );
                error.throw();
            end
            fclose(fileID);
            
            % Define filepath
            self.filepath = newFilepath;
        end
        
        function set.email(self,newEmailFilepath)
            [fileID,message] = fopen(newEmailFilepath, 'a');
            if(fileID < 0)
                error = MException( 'Logger:setEmail' ...
                                  , ['Problem with supplied filepath: ' message] );
                error.throw();
            end
            fclose(fileID);
            
            % Define filepath
            self.email = utils.Email.EmailFromFile(newEmailFilepath);
        end
    end
    
    methods (Static)
        function self = getLogger(logFilepath,commandWindowLevel,logLevel,emailLogLevel,emailFilepath,emailRecipients)
            %
            % Obtain logger object.
            %
            % INPUTS:
            %   1. logFilepath - filepath for logging file.
            %   2. commandWindowLevel - log level for command window.
            %   3. logLevel      - log level for logging file.
            %   4. emailLogLevel - log level for emailing message.
            %   5. emailFilepath   - email configuration filepath.
            %   6. emailRecipients - email recipients.
            %
            narginchk(0,6);
            
            persistent localLogger
            
            if nargin < 1
                logFilepath = 'logging.log';
            end
            
            if nargin < 2
                commandWindowLevel = utils.logging.LogEventTypes.INFO;
            end
            
            if nargin < 3
                logLevel = utils.logging.LogEventTypes.ALL;
            end
            
            if nargin < 4
                emailLogLevel = utils.logging.LogEventTypes.ERROR;
            end
            
            if nargin < 5
                emailFilepath = '';
            end
            
            if nargin < 6
                emailRecipients = '';
            end
            
            % Initialize logger
            if isempty(localLogger) || ~isvalid(localLogger)
                localLogger = utils.logging.Logger();
                localLogger.filepath = logFilepath;
                localLogger.logLevel = logLevel;
                if ~isempty(emailFilepath)
                    localLogger.email = utils.Email.EmailFromFile(emailFilepath);
                    localLogger.email.recipients = emailRecipients;
                end
                localLogger.emailLogLevel      = emailLogLevel;
                localLogger.commandWindowLevel = commandWindowLevel;
            end
            self = localLogger;
        end
        
        function testLoggingSpeed()
            %
            % Test logging speed for logger.
            %

            % Obtain logger
            logger   = utils.logging.Logger.getLogger();
            listener = utils.logging.LoggerListener(logger);
            
            % Start tests
            disp('Testing speed logging only to command window - 1e4 logs');
            logger.commandWindowLevel = utils.logging.LogEventTypes.TRACE;
            logger.logLevel           = utils.logging.LogEventTypes.OFF;
            logger.emailLogLevel      = utils.logging.LogEventTypes.OFF;
            tic;
            for iLog = 1:1e4
                logEvent = utils.logging.LogEvent();
                logEvent.type    = utils.logging.LogEventTypes.TRACE;
                logEvent.message = 'test_logging_to_command_window';
                notify(logger,'logEvent',logEvent);
            end
            toc;
            
            
            disp('Testing speed logging only to file - 1e4 logs');
            logger.commandWindowLevel = utils.logging.LogEventTypes.OFF;
            logger.logLevel           = utils.logging.LogEventTypes.TRACE;
            logger.emailLogLevel      = utils.logging.LogEventTypes.OFF;
            tic;
            for iLog = 1:1e4
                logEvent = utils.logging.LogEvent();
                logEvent.type    = utils.logging.LogEventTypes.TRACE;
                logEvent.message = 'test_logging_to_file';
                notify(logger,'logEvent',logEvent);
            end
            toc;
           
 
            disp('Testing speed when logging is off - 1e4 logs');
            logger.commandWindowLevel = utils.logging.LogEventTypes.OFF;
            logger.logLevel           = utils.logging.LogEventTypes.OFF;
            logger.emailLogLevel      = utils.logging.LogEventTypes.OFF;
            tic;
            for iLog = 1:1e4
                logEvent = utils.logging.LogEvent();
                logEvent.type    = utils.logging.LogEventTypes.TRACE;
                logEvent.message = 'test_logging_off';
                notify(logger,'logEvent',logEvent);
            end
            toc;
            
            % Delete listener
            listener.delete();
        end
    end
end
