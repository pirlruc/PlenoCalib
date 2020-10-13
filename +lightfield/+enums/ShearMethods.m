classdef ShearMethods
    %SHEARMETHODS
    %   Package to define lightfield shearing methods available.
    
    properties (SetAccess = immutable)
        symbol          % Matching cost symbol
        description     % Description of matching cost method
    end
   
    methods
        function self = ShearMethods(symbol,description)
            %
            % Lightfield shear methods instance.
            %
            self.symbol      = symbol;
            self.description = description;
        end
    end
    
    enumeration
        FREQUENCY('frequency', 'Lightfield shearing using fourier transform properties.' )
        CUBIC    ('cubic'    , [ 'Cubic interpolation - the interpolated value at a query point is based on ' ...
                                 'a cubic interpolation of the values at neighboring grid points in each ' ...
                                 'respective dimension. The interpolation is based on a cubic convolution.' ]);
        NEAREST  ('nearest'  , [ 'Nearest-neighbor interpolation - the output pixel is assigned the value ' ...
                                 'of the pixel that the point falls within. No other pixels are considered.' ]);
        LINEAR   ('linear'   , [ 'Linear interpolation - the interpolated value at a query point is based on ' ...
                                 'a linear interpolation of the values at neighboring grid points in each ' ...
                                 'respective dimension.' ]);
    end
end