classdef FilterCollection
    %FILTERCOLLECTION
    %   Methods and properties for a collection of filters.
    
    properties
        filters = []        % List of filters
    end
    
    properties (Dependent)
        mask                % Resulting mask of applying the filters in collection
        numberFilters       % Number of filters in collection
        boundary            % Boundary conditions considered for collection
        type                % Algorithm considered to apply filtering for collection
    end
    
    methods
        function self = FilterCollection(varargin)
            %
            % Create filter collection instance.
            %
            % INPUTS:
            %   1. filters - list of filters. The filters should have the 
            %   same properties regarding the boundary conditions and the
            %   algorithm used to apply the resulting mask.
            %
            narginchk(0,1);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.filters = varargin{1};
                end
            end
        end
        
        function self = set.filters(self,newFilters)
            self.filters = utils.enums.Classes.FILTER().convert(newFilters);
        end
        
        function numberFilters = get.numberFilters(self)
            numberFilters = length(self.filters);
        end
        
        function boundary = get.boundary(self)
            % If there are no filters in collection, consider the default
            % boundary conditions.
            if isempty(self.filters)
                boundary = image.filters.Filter().boundary;
            
            % Otherwise, obtain the boundary conditions of a filter. The
            % filters shoud have the same boundary conditions.
            else
                boundary = self.filters(1).boundary;
            end
        end
        
        function type = get.type(self)
            % If there are no filters in collection, consider the default
            % algorithm to apply the filtering.
            if isempty(self.filters)
                type = image.filters.Filter().type;
            
            % Otherwise, obtain the filtering algorithm of a filter. The
            % filters shoud have the same filtering algorithm.
            else
                type = self.filters(1).type;
            end
        end
                
        function mask = get.mask(self)
            % Since the convolution of filters is faster than the
            % convolution of filters with images, process first the filters
            % with themselves and then apply to image.
            %
            % Notice that the process of creating a mask returns the 
            % central part of the convolution with the same size as the
            % first mask in collection.
            %
            
            % If there are no filters in collection, return the default
            % mask
            if self.numberFilters == 0
                mask = image.filters.Filter().mask;

            % Otherwise, obtain the final mask by sequentially convolving
            % the filters in collection.
            else
                % Convolution needs masks in matlab format
                % Obtain matlab format and convert masks
                matlabFormat = image.enums.Formats.MATLAB_FORMAT();
                mask = matlabFormat.encode(self.filters(1).mask);
                for iFilter = 2:self.numberFilters
                    mask = conv2( mask ...
                                , matlabFormat.encode(self.filters(iFilter).mask) ...
                                , 'same');
                end

                % Set to zero components that have value less that numerical
                % precision
                mask(abs(mask) < eps) = 0;

                % Normalize filter to apply to image using the L1-norm.
                mask = mask ./ sum(abs(mask(:)));
                
                % Convert to internal format
                mask = matlabFormat.decode(mask);
            end
        end
    end
    
    methods
        function self = addFilters(self, newFilters)
            %
            % Add filters to the collection of filters.
            %
            % INPUTS:
            %   1. newFilters - new filters to be added to the collection. 
            %   The filters should have the same properties regarding the 
            %   boundary conditions and the algorithm used to apply the 
            %   resulting mask.
            %
            narginchk(2,2);
            
            % Convert new filters to filters data
            newFilters = utils.enums.Classes.FILTER().convert(newFilters);
            
            % Add filters to collection
            self.filters = cat(1, self.filters, newFilters);
        end
        
        function filteredImage = apply(self,imageData)
            %
            % Apply collection of filters to image.
            %
            % INPUTS:
            %   1. imageData - image data. The image should be formatted in
            %   the image internal format.
            %
            narginchk(2,2);
            
            % Obtain new filter with resulting mask
            filter = image.filters.Filter(self.mask, self.boundary, self.type);
            
            % Apply filter to image
            filteredImage = filter.apply(imageData);
        end
    end
end
