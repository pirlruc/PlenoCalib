classdef Time
    %TIME
    %   Utility for measuring time and converting time formats.
    
    properties (Constant, Hidden)
        SECONDS_PER_MINUTE = 60;                                 % Number of seconds per minute.
        SECONDS_PER_HOUR   = 60 * utils.Time.SECONDS_PER_MINUTE; % Number of seconds per hour.
        SECONDS_PER_DAY    = 24 * utils.Time.SECONDS_PER_HOUR;   % Number of seconds per day.
    end
    
    properties (~Dependent, Hidden, SetAccess = private, GetAccess = private)
        counterID = uint64(zeros)   % Counter identification associated with the time measurement.
    end
    
    properties (~Dependent)
        time = 0                    % Measured time in seconds.
    end
    
    properties (Dependent)
        days            % Number of days.          
        hours           % Number of hours. Days are not considered.
        minutes         % Number of minutes. Days and hours are not considered.
        seconds         % Number of seconds. Days, hours and minutes are not considered.
        miliseconds     % Number of miliseconds. Days, hours, minutes and seconds are not considered.
    end
    
    methods
        function self = Time(varargin)
            % 
            % Time counter instance.
            %
            % INPUTS:
            %   1. time - time in seconds.
            %
            narginchk(0,1);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.time = varargin{1};
                end
            end
        end
        
        function self = set.time(self,newTime)
            self.time = newTime;
        end
        
        function days = get.days(self)
            days = floor(self.time / utils.Time.SECONDS_PER_DAY);
        end
        
        function hours = get.hours(self)
            remainingSeconds = self.time - self.days * utils.Time.SECONDS_PER_DAY;
            hours            = floor(remainingSeconds / utils.Time.SECONDS_PER_HOUR);
        end
        
        function minutes = get.minutes(self)
            remainingSeconds = self.time - self.days  * utils.Time.SECONDS_PER_DAY ...
                                         - self.hours * utils.Time.SECONDS_PER_HOUR;
            minutes = floor(remainingSeconds / utils.Time.SECONDS_PER_MINUTE);
        end
        
        function seconds = get.seconds(self)
            remainingSeconds = self.time - self.days    * utils.Time.SECONDS_PER_DAY ...
                                         - self.hours   * utils.Time.SECONDS_PER_HOUR ...
                                         - self.minutes * utils.Time.SECONDS_PER_MINUTE;
            seconds = floor(remainingSeconds);
        end
        
        function miliseconds = get.miliseconds(self)
            miliseconds = self.time - self.days    * utils.Time.SECONDS_PER_DAY ...
                                    - self.hours   * utils.Time.SECONDS_PER_HOUR ...
                                    - self.minutes * utils.Time.SECONDS_PER_MINUTE ...
                                    - self.seconds;
        end
    end
    
    methods
        function self = startTiming(self)
            %
            % Start measuring time.
            %
            self.counterID = tic;
        end
        
        function self = stopTiming(self)
            %
            % Stop measuring time.
            %
            self.time = toc(self.counterID);
        end
    end
end