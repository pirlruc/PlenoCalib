classdef CheckerboardCorners
    %CHECKERBOARDCORNERS
    %   Checkerboard corners utility.
    
    properties
        corners           = image.Pixel()   % Corners location
        positionInPattern = math.Vector()   % Position of corner in pattern
    end
    
    properties (Dependent)
        cornersIntoPattern                  % Corners location associated with pattern position
        patternSize                         % Pattern size
        numberCorners                       % Number corners
    end
    
    methods
        function self = CheckerboardCorners(varargin)
            %
            % Create checkerboard corners instance.
            %
            % INPUTS:
            %   1. corners           - corners location on image.
            %   2. positionInPattern - position of corners in pattern.
            %
            narginchk(0,2);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.corners = varargin{1};
                end
                
                if nargin >= 2
                    self.positionInPattern = varargin{2};
                end
            end
        end
        
        function self = set.corners(self,newCorners)
            self.corners = utils.enums.Classes.PIXEL().convert(newCorners);
        end
        
        function self = set.positionInPattern(self,newPositionInPattern)
            self.positionInPattern = utils.enums.Classes.VECTOR().convert(newPositionInPattern);
        end
        
        function number = get.numberCorners(self)
            number = self.corners.numberVectors;
        end
        
        function patternSize = get.patternSize(self)
            if self.numberCorners > 0
                positionInPattern_u = self.positionInPattern.obtainComponents(1);
                positionInPattern_v = self.positionInPattern.obtainComponents(2);
                
                patternSize = image.Pixel( [ length(unique(positionInPattern_u.data)) ...
                                           , length(unique(positionInPattern_v.data)) ], false );
            else
                patternSize = image.Pixel();
            end
        end
        
        function pattern = get.cornersIntoPattern(self)
            % Reorganize corners coordinates if corners are defined
            pattern = [];
            if self.numberCorners > 0
                % Obtain minimum position in pattern. This will allow to
                % translate the positions for indexes starting in (1,1).
                positionInPattern_u = self.positionInPattern.obtainComponents(1);
                positionInPattern_v = self.positionInPattern.obtainComponents(2);
                minimumPosition_u   = min(positionInPattern_u.data);
                minimumPosition_v   = min(positionInPattern_v.data);

                % Define corners location on pattern and define linear indices
                cornersLocation = sub2ind( [self.patternSize.u,self.patternSize.v] ...
                                         , positionInPattern_u.data - minimumPosition_u + 1 ...
                                         , positionInPattern_v.data - minimumPosition_v + 1 );

                % Initialize corners pattern
                pattern = nan(self.corners.numberComponents,self.patternSize.u,self.patternSize.v);

                % Obtain pattern
                pattern(1,cornersLocation) = self.corners.u;
                pattern(2,cornersLocation) = self.corners.v;
            end
        end
    end
    
    methods
        function show(self,imageData)
            %
            % Show corners on image.
            %
            % INPUTS:
            %   1. imageData - image data.
            %
            narginchk(2,2);
            
            % Convert to image object
            imageData = utils.enums.Classes.IMAGE().convert(imageData);
            
            % Obtain corners in pattern
            pattern = self.cornersIntoPattern;
            
            % Obtain corners at the right and at the bottom
            rightCorners  = nan(size(pattern));
            bottomCorners = nan(size(pattern));
            rightCorners (:,1:end-1,:) = pattern(:,2:end,:);
            bottomCorners(:,:,1:end-1) = pattern(:,:,2:end);
            
            % Display image
            imageData.show;
            hold on;
            plot([pattern(1,:);rightCorners(1,:) ],[pattern(2,:);rightCorners(2,:) ],'g-');
            plot([pattern(1,:);bottomCorners(1,:)],[pattern(2,:);bottomCorners(2,:)],'c-');
            plot(self.corners.u,self.corners.v,'r.');
            hold off;
        end
    end
    
    methods (Static)
        function self = CheckerboardCornersFromCheckerboardImage(imageData)
            %
            % Obtain corners from checkerboard image.
            %
            % INPUTS:
            %   1. imageData - checkerboard corner image.
            %
            narginchk(1,1);
            
            % Convert to image instance
            imageData = utils.enums.Classes.IMAGE().convert(imageData);
            
            % Obtain gray image and transform to matlab format
            imageData = imageData.rgb2gray();
            grayImage_matlabFormat = imageData.MATLAB_FORMAT().encode(imageData.data);
            
            % The corner detector requires that the image must be in uint8
            % and with double precision
            grayImage_matlabFormat = double(im2uint8(grayImage_matlabFormat));
            
            % Obtain corners using Bok's corner detector
            corners = CheckerboardCorner(grayImage_matlabFormat);
            
            % Create checkerboard corners object
            self = image.corners.CheckerboardCorners();
            if ~isempty(corners)
                self.corners           = corners(3:4,:);
                self.positionInPattern = corners(1:2,:);
            end
        end
    end
    
    methods (Static)
        function corners = ensureTopLeftOrigin(corners,patternSize)
            %
            % Reorganize corners in order to have the origin of the corners
            % in the top left corner.
            %
            % INPUT: 
            %   1. corners - ordered corners for checkerboard pattern.
            %   2. patternSize - checkerboard pattern size.
            %
            narginchk(2,2);

            % Convert to pixel object
            corners     = utils.enums.Classes.PIXEL().convert(corners);
            patternSize = utils.enums.Classes.PIXEL().convert(patternSize);
            
            % Obtain checkerboard centroid
            patternCentroid = image.Pixel([mean(corners.u);mean(corners.v)],false);

            % Top right corner origin correction
            isTopRight = corners.u(1) > patternCentroid.u && corners.v(1) < patternCentroid.v;
            if isTopRight
                pattern = reshape(corners.data, [2, patternSize.u, patternSize.v]);
                pattern = pattern(:, end:-1:1, :);
                pattern = permute(pattern, [1,3,2]);
                pattern = reshape(pattern, 2, []);
                corners.data = pattern;
            end

            % Bottom left corner origin correction
            isBottomLeft = corners.u(1) < patternCentroid.u && corners.v(1) > patternCentroid.v;
            if isBottomLeft
                pattern = reshape(corners.data, [2, patternSize.u, patternSize.v]);
                pattern = pattern(:, :, end:-1:1);
                pattern = permute(pattern, [1,3,2]);
                pattern = reshape(pattern, 2, []);
                corners.data = pattern;
            end

            % Bottom right corner origin correction
            isBottomRight = corners.u(1) > patternCentroid.u && corners.v(1) > patternCentroid.v;
            if isBottomRight
                pattern = reshape(corners.data, [2, patternSize.u, patternSize.v]);
                pattern = pattern(:, end:-1:1, end:-1:1);
                pattern = permute(pattern, [1,3,2]);
                pattern = reshape(pattern, 2, []);
                corners.data = pattern;
            end

            % If the origin is still not top left, display warning
            origin    = corners.obtainVectors(1);
            isTopLeft = all(origin.data < patternCentroid.data);
            if ~isTopLeft
                warning( 'CheckerboardCorners:incorrectOrigin' ...
                       , 'Origin of checkerboard pattern is not the top left corner.' );
            end
        end
        
        function corners = ensureBottomRightFromRightBottomOrdering(corners,patternSize)
            %
            % Obtain corners ordered from top to bottom and then right.
            %
            % INPUTS:
            %   1. corners     - checkerboard corners information.
            %   2. patternSize - checkerboard pattern size.
            % 
            narginchk(2,2);

            % Convert to pixel object
            corners     = utils.enums.Classes.PIXEL().convert(corners);
            patternSize = utils.enums.Classes.PIXEL().convert(patternSize);
            
            % Ensure the top left corner as the origin for the pattern
            corners = image.corners.CheckerboardCorners.ensureTopLeftOrigin(corners,patternSize);
                
            % Order corners in order to have origin on its top left corner
            % and have the following direction: down and right. This 
            % ordering corresponds to the matlab ordering for the function
            % detectCheckerboardPoints. The current ordering of the corner
            % is: right and down.
            newCorners  = [];
            for iCorner = 1:patternSize.u
                % Obtain indices for sampling and concatenate corners
                indices    = iCorner:patternSize.u:corners.numberVectors;
                newCorners = cat(2,newCorners,corners.obtainVectors(indices).data);
            end
            
            % Convert to pixel object
            corners = image.Pixel(newCorners,false);
        end
        
        function points = obtainGroundTruthCheckerboardCorners(patternSize,patternLength)
            %
            % Obtain ground truth checkerboard corners assuming the top
            % left corner as the origin and that the direction to describe
            % the corners is: bottom and right.
            %
            % INPUTS:
            %   1. patternSize   - number of corners in pattern
            %   2. patternLength - spacing between corners in pattern in
            %   metric units.
            %
            narginchk(2,2);
            
            % Convert to pixel object
            patternSize   = utils.enums.Classes.PIXEL().convert(patternSize);
            patternLength = utils.enums.Classes.PIXEL().convert(patternLength);
            
            % Obtain x and y ground truth coordinates
            pattern_x = patternLength.u .* (0:patternSize.u - 1);
            pattern_y = patternLength.v .* (0:patternSize.v - 1);
            
            % Obtain grid
            [pattern_y, pattern_x] = ndgrid(pattern_y, pattern_x);
            points = math.Point([pattern_x(:),pattern_y(:),zeros(numel(pattern_x),1)]',false);
        end
    end        
end

