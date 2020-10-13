classdef InterpolationMethods
    %INTERPOLATIONMETHODS
    %   Image interpolation methods.
    
    properties (~Dependent, SetAccess = immutable, GetAccess = public)
        method          % Interpolation method
        description     % Description of interpolation method
    end
    
    methods
        function self = InterpolationMethods(method, description)
            %
            % Interpolation method enumeration instance.
            %
            self.method      = method;
            self.description = description;
        end
    end
    
    enumeration
        BICUBIC  ('bicubic' , [ 'Bicubic interpolation - the output pixel value is a weighted average ' ...
                                'of pixels in the nearest 4-by-4 neighborhood.' ]);
        CUBIC    ('cubic'   , [ 'Cubic interpolation - the interpolated value at a query point is based on ' ...
                                'a cubic interpolation of the values at neighboring grid points in each ' ...
                                'respective dimension. The interpolation is based on a cubic convolution.' ]);
        NEAREST  ('nearest' , [ 'Nearest-neighbor interpolation - the output pixel is assigned the value ' ...
                                'of the pixel that the point falls within. No other pixels are considered.' ]);
        BILINEAR ('bilinear', [ 'Bilinear interpolation - the output pixel value is a weighted average ' ...
                                'of pixels in the nearest 2-by-2 neighborhood.' ]);
        LINEAR   ('linear'  , [ 'Linear interpolation - the interpolated value at a query point is based on ' ...
                                'a linear interpolation of the values at neighboring grid points in each ' ...
                                'respective dimension.' ]);
        BOX      ('box'     , 'Box-shaped kernel.');
        LANCZOS2 ('lanczos2', 'Lanczos-2 kernel.');
        LANCZOS3 ('lanczos3', 'Lanczos-3 kernel.');
    end
end