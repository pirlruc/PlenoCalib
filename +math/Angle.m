classdef Angle
    %ANGLE
    %   Angle methods and properties.
    
    properties
        anglesInRadians = []     % List of angles values in radians
    end
    
    properties (Dependent)
        anglesInDegrees          % List of angle values in degrees
    end
    
    methods
        function self = Angle(varargin)
            %
            % Angle instance.
            %
            % INPUTS:
            %   1. angle - list of angle values in radians.
            %
            narginchk(0,1);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.anglesInRadians = varargin{1};
                end
            end
        end
        
        function self = set.anglesInRadians(self, newAngles)
            % Angles must be provided in radians
            self.anglesInRadians = newAngles;
        end
        
        function angles = get.anglesInDegrees(self)
            %
            % Convert angle values from radians to degrees.
            %
            angles = radtodeg(self.anglesInRadians);
        end
    end
end
