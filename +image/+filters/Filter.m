classdef Filter
    %FILTER
    %   Filter methods and properties.
    
    properties
        mask     = []                                           % Filter mask
        boundary = image.enums.BoundaryConditions.SYMMETRIC     % Boundary conditions
        type     = image.enums.FilteringAlgorithms.CONVOLUTION  % Algorithm to apply filter
    end
    
    methods
        function self = Filter(varargin)
            %
            % Create filter instance.
            %
            % INPUTS:
            %   1. mask - filter mask.
            %   2. boundaryConditions - boundary conditions to be
            %   considered. 
            %   3. filteringAlgorithm - algorithm to apply filter.
            %
            narginchk(0,3)
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.mask = varargin{1};
                end
                
                if nargin >= 2
                    self.boundary = varargin{2};
                end
                
                if nargin >= 3
                    self.type = varargin{3};
                end
            end
        end
        
        function self = set.mask(self,newMask)
            self.mask = newMask;
        end
        
        function self = set.boundary(self, newBoundaryConditions)
            self.boundary = newBoundaryConditions;
        end
        
        function self = set.type(self, newFilteringAlgorithm)
            self.type = newFilteringAlgorithm;
        end
    end
    
    methods
        function filteredImage = apply(self,imageData)
            %
            % Apply filter to image. The method allows to apply filtering
            % to color images. It assumes symmetric boundary conditions.
            %
            % INPUTS:
            %   1. imageData - image data. The image should be formatted in
            %   the image internal format.
            %
            narginchk(2,2);
            
            % If mask is empty, throw error
            if isempty(self.mask)
                error = MException('Filter:apply:noMask', 'Mask not defined...');
                error.throw();
            end
            
            % Convert to image data
            imageInstance = utils.enums.Classes.IMAGE().convert(imageData);
            
            % Obtain matlab format
            matlabFormat = image.enums.Formats.MATLAB_FORMAT();
            
            % Encode matlab format and filter image
            filteredImage = imfilter( matlabFormat.encode(imageInstance.data) ...
                                    , matlabFormat.encode(self.mask) ...
                                    , self.boundary.method ...
                                    , self.type.algorithm ...
                                    , 'same');

            % Set to zero components that have value less that numerical
            % precision
            filteredImage(abs(filteredImage) < eps) = 0;

            % Convert back to internal format and create image instance
            filteredImage = image.Image(matlabFormat.decode(filteredImage));
        end
    end
end