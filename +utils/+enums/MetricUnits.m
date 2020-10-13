classdef MetricUnits
    %METRICUNITS
    %   Enumerates the metric units available
    
    properties
        description         % Description of metric unit
        symbol              % Symbol of metric unit
        power               % Metric unit power relatively to meter
    end
    
    methods
        function self = MetricUnits(description,symbol,power)
            %
            % Metric units instance.
            %
            self.description = description;
            self.symbol      = symbol;
            self.power       = power;
        end
    end
    
    methods
        function factor = obtainConversionFactor(self,unit)
            %
            % Obtain conversion factor between two metric units.
            %
            % INPUTS:
            %   1. unit - target metric unit. Default is meter.
            %
            narginchk(1,2);
            
            if nargin <= 1
                unit = utils.enums.MetricUnits.METER();
            end
            
            % Obtain conversion factor
            factor = 10^(unit.power - self.power);
        end
    end
    
    enumeration
        METER     ( 'Meter'     , 'm' , 0 ) 
        CENTIMETER( 'Centimeter', 'cm', 2 )
        MILIMETER ( 'Milimeter' , 'mm', 3 )
    end
end

