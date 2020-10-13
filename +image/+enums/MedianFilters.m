classdef MedianFilters
    %MEDIANFILTERS
    %   Median filters available.
    
    properties (~Dependent, SetAccess = immutable, GetAccess = public)
        filter          % Median filter
        description     % Description of median filter
    end
    
    methods
        function self = MedianFilters(filter, description)
            %
            % Median filters enumeration instance.
            %
            self.filter      = filter;
            self.description = description;
        end
    end
    
    enumeration
        MEDIAN_2D ('medfilt2(%s,%s)' , [ 'Median filter 2D. This filter can be applied to grayscale images or ' ...
                                       , 'to each color channel independently.' ]);
        MEDIAN_3D ('medfilt3(%s,%s)' , 'Median filter 3D. This filter can be applied to color images.');
    end
end