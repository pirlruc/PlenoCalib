classdef Colors
    %COLORS
    %   Enumerates the colors available
    
    properties
        color       % Description of color
        symbol      % Matlab symbol for color
    end
    
    methods
        function self = Colors(color,symbol)
            %
            % Colors instance.
            %
            self.color  = color;
            self.symbol = symbol;
        end
    end
    
    enumeration
        RED  ( 'Red'  , 'r' ) 
        BLUE ( 'Blue' , 'b' )
        CYAN ( 'Cyan' , 'c' )
        GREEN( 'Green', 'g' )
        BLACK( 'Black', 'k' )
        WHITE( 'White', 'w' )
        YELLOW ( 'Yellow' , 'y' )
        MAGENTA( 'Magenta', 'm' )
    end
end

