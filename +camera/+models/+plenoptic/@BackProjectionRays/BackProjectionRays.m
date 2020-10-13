classdef BackProjectionRays
    %BACKPROJECTIONRAYS
    %   Utility to display back-projection rays from plenoptic cameras.

    properties (Constant)
        NUMBER_POINTS_PER_RAY = 4   % Number of points to represent ray in a lightfield
    end
    
    properties
        rays      = math.Point()    % Back-projection rays considering the center of the pixels
        outerRays = math.Point()    % Back-projection rays considering the outer border of the pixels
        innerRays = math.Point()    % Back-projection rays considering the inner border of the pixels
    end
    
    properties (Dependent)
        numberRays              % Number of rays to be represented
        rays2D                  % 2D representation of the back-projection rays for the center of the pixels
        outerRays2D             % 2D representation of the back-projection rays for the outer border of the pixels
        innerRays2D             % 2D representation of the back-projection rays for the inner border of the pixels
    end
    
    methods
        function self = BackProjectionRays(varargin)
            %
            % Create back-projection rays instance.
            %
            % INPUTS:
            %   1. rays      - back-projection rays.
            %   2. innerRays - back-projection rays considering the outer 
            %      border of the pixels.
            %   3. outerRays - back-projection rays considering the inner 
            %      border of the pixels
            %
            % Each ray point should be given in different columns.
            %
            narginchk(0,3);
            
            if ~isempty(varargin)
                if nargin >= 3
                    self.outerRays = varargin{3};
                end
                
                if nargin >= 2
                    self.innerRays = varargin{2};
                end
                
                if nargin >= 1
                    self.rays = varargin{1};
                end
            end
        end
        
        function self = set.rays(self,newRays)
            self.rays = utils.enums.Classes.POINT().convert(newRays);
        end
        
        function self = set.innerRays(self,newInnerRays)
            self.innerRays = utils.enums.Classes.POINT().convert(newInnerRays);
        end
        
        function self = set.outerRays(self,newOuterRays)
            self.outerRays = utils.enums.Classes.POINT().convert(newOuterRays);
        end
        
        function number = get.numberRays(self)
            number = self.rays.numberVectors / self.NUMBER_POINTS_PER_RAY;
        end
        
        function rays2D = get.rays2D(self)
            if isempty(self.rays.data)
                rays2D = [];
            else
                rays2D = unique([self.rays.x;self.rays.z],'rows');
            end
        end
        
        function rays2D = get.outerRays2D(self)
            if isempty(self.outerRays.data)
                rays2D = [];
            else
                rays2D = unique([self.outerRays.x;self.outerRays.z],'rows');
            end
        end
        
        function rays2D = get.innerRays2D(self)
            if isempty(self.innerRays.data)
                rays2D = [];
            else
                rays2D = unique([self.innerRays.x;self.innerRays.z],'rows');
            end
        end
    end
    
    methods
        function plot(self,color)
            %
            % Display back-projection rays in 2D.
            %
            % INPUTS
            %   1. color - color to display rays. Default is black.
            %
            narginchk(1,2);
            
            if nargin <= 1
                color = utils.enums.Colors.BLACK();
            end
            
            % To display the back-projection rays, we need to have at least
            % one ray defined
            if self.numberRays == 0
                error = MException('BackProjectionRays:plot:noRays', 'Rays not provided...');
                error.throw();
            end

            % Obtain 2D representation of the rays
            oRays2D      = math.Vector(self.rays2D);
            oInnerRays2D = math.Vector(self.innerRays2D);
            oOuterRays2D = math.Vector(self.outerRays2D);
            
            % Display middle pixel rays
            gcf;
            hold on;
            plot( oRays2D.obtainComponents(1).data ...
                , oRays2D.obtainComponents(2).data, 'xr');
            
            % Decide if display for back-projection rays should include
            % inner and outer rays
            % Plot only middle pixels
            if oOuterRays2D.numberVectors == 0 ...
            && oInnerRays2D.numberVectors == 0
                for iRay = 1:self.numberRays
                    % Plot middle pixel
                    plot( oRays2D.obtainComponents(1).data(iRay:self.numberRays:end) ...
                        , oRays2D.obtainComponents(2).data(iRay:self.numberRays:end), color.symbol );
                end
                
            % Plot middle and inner pixels
            elseif oOuterRays2D.numberVectors == 0
                plot( oInnerRays2D.obtainComponents(1).data ...
                    , oInnerRays2D.obtainComponents(2).data, 'xr');
                for iRay = 1:self.numberRays
                    % Plot middle pixel
                    plot( oRays2D.obtainComponents(1).data(iRay:self.numberRays:end) ...
                        , oRays2D.obtainComponents(2).data(iRay:self.numberRays:end), color.symbol );

                    % Plot inner pixel
                    plot( oInnerRays2D.obtainComponents(1).data(iRay:self.numberRays:end) ...
                        , oInnerRays2D.obtainComponents(2).data(iRay:self.numberRays:end), 'b' );
                end
                
            % Plot outer and middle pixels
            elseif oInnerRays2D.numberVectors == 0
                plot( oOuterRays2D.obtainComponents(1).data ...
                    , oOuterRays2D.obtainComponents(2).data, 'xr');
                for iRay = 1:self.numberRays
                    % Plot outer pixel
                    plot( oOuterRays2D.obtainComponents(1).data(iRay:self.numberRays:end) ...
                        , oOuterRays2D.obtainComponents(2).data(iRay:self.numberRays:end), 'g' );

                    % Plot middle pixel
                    plot( oRays2D.obtainComponents(1).data(iRay:self.numberRays:end) ...
                        , oRays2D.obtainComponents(2).data(iRay:self.numberRays:end), color.symbol );
                end
                
            % Plot rays with information from depth of field
            else
                plot( oOuterRays2D.obtainComponents(1).data ...
                    , oOuterRays2D.obtainComponents(2).data, 'xr');
                plot( oInnerRays2D.obtainComponents(1).data ...
                    , oInnerRays2D.obtainComponents(2).data, 'xr');
                for iRay = 1:self.numberRays
                    % Plot outer pixel
                    plot( oOuterRays2D.obtainComponents(1).data(iRay:self.numberRays:end) ...
                        , oOuterRays2D.obtainComponents(2).data(iRay:self.numberRays:end), 'g' );

                    % Plot middle pixel
                    plot( oRays2D.obtainComponents(1).data(iRay:self.numberRays:end) ...
                        , oRays2D.obtainComponents(2).data(iRay:self.numberRays:end), color.symbol );

                    % Plot inner pixel
                    plot( oInnerRays2D.obtainComponents(1).data(iRay:self.numberRays:end) ...
                        , oInnerRays2D.obtainComponents(2).data(iRay:self.numberRays:end), 'b' );
                end
            end
            hold off;
        end
        
        function plot3(self,color)
            %
            % Display back-projection rays in 3D.
            %
            % INPUTS
            %   1. color - color to display rays. Default is black.
            %
            narginchk(1,2);
            
            if nargin <= 1
                color = utils.enums.Colors.BLACK();
            end
            
            % To display the back-projection rays, we need to have at least
            % one ray defined
            if self.numberRays == 0
                error = MException('BackProjectionRays:plot3:noRays', 'Rays not provided...');
                error.throw();
            end
            
            % Display middle pixel rays
            gcf;
            hold on;
            plot3(self.rays.x, self.rays.y, self.rays.z, 'xr');
            
            % Decide if display for back-projection rays should include
            % inner and outer rays
            % Plot only middle pixels
            if self.outerRays.numberVectors == 0 ...
            && self.innerRays.numberVectors == 0
                for iRay = 1:self.numberRays
                    % Plot middle pixel
                    plot3( self.rays.x(iRay:self.numberRays:end) ...
                         , self.rays.y(iRay:self.numberRays:end) ...
                         , self.rays.z(iRay:self.numberRays:end), color.symbol );
                end
                
            % Plot middle and inner pixels
            elseif self.outerRays.numberVectors == 0
                plot3(self.innerRays.x, self.innerRays.y, self.innerRays.z, 'xr');
                for iRay = 1:self.numberRays
                    % Plot middle pixel
                    plot3( self.rays.x(iRay:self.numberRays:end) ...
                         , self.rays.y(iRay:self.numberRays:end) ...
                         , self.rays.z(iRay:self.numberRays:end), color.symbol );

                    % Plot inner pixel
                    plot3( self.innerRays.x(iRay:self.numberRays:end) ...
                         , self.innerRays.y(iRay:self.numberRays:end) ...
                         , self.innerRays.z(iRay:self.numberRays:end), 'b' );
                end
                
            % Plot outer and middle pixels
            elseif self.innerRays.numberVectors == 0
                plot3(self.outerRays.x, self.outerRays.y, self.outerRays.z, 'xr');
                for iRay = 1:self.numberRays
                    % Plot outer pixel
                    plot3( self.outerRays.x(iRay:self.numberRays:end) ...
                         , self.outerRays.y(iRay:self.numberRays:end) ...
                         , self.outerRays.z(iRay:self.numberRays:end), 'g' );

                    % Plot middle pixel
                    plot3( self.rays.x(iRay:self.numberRays:end) ...
                         , self.rays.y(iRay:self.numberRays:end) ...
                         , self.rays.z(iRay:self.numberRays:end), color.symbol );
                end
                
            % Plot rays with information from depth of field
            else
                plot3(self.outerRays.x, self.outerRays.y, self.outerRays.z, 'xr');
                plot3(self.innerRays.x, self.innerRays.y, self.innerRays.z, 'xr');
                for iRay = 1:self.numberRays
                    % Plot outer pixel
                    plot3( self.outerRays.x(iRay:self.numberRays:end) ...
                         , self.outerRays.y(iRay:self.numberRays:end) ...
                         , self.outerRays.z(iRay:self.numberRays:end), 'g' );

                    % Plot middle pixel
                    plot3( self.rays.x(iRay:self.numberRays:end) ...
                         , self.rays.y(iRay:self.numberRays:end) ...
                         , self.rays.z(iRay:self.numberRays:end), color.symbol );

                    % Plot inner pixel
                    plot3( self.innerRays.x(iRay:self.numberRays:end) ...
                         , self.innerRays.y(iRay:self.numberRays:end) ...
                         , self.innerRays.z(iRay:self.numberRays:end), 'b' );
                end
            end
            hold off;
        end
    end
    
    methods (Static)
        function self = BackProjectionRaysFromIntrinsicMatrix(intrinsicMatrixInstance,imageRays_ijkl)
            %
            % Create back-projection rays instance from a plenoptic camera 
            % intrinsic matrix.
            %
            % Obtain the coordinates of the back-projection rays. These
            % coordinates correspond to the image, microlens, main lens and
            % the world plane in metric units.
            %
            % INPUTS:
            %   1. intrinsicMatrixInstace - intrinsic matrix instance for
            %      plenoptic cameras. Default is a standard plenoptic
            %      camera intrinsic matrix.
            %   2. imageRays_ijkl - image space rays in pixels and
            %      microlenses indices. Each ray should be provided in 
            %      different columns.
            %
            narginchk(0,2);
            
            % If intrinsic matrix instance is not provided, obtain example
            % of intrinsic matrix with standard plenoptic camera.
            if nargin <= 0
                intrinsicMatrixInstance = camera.models.plenoptic.standard.IntrinsicMatrix();
            end
            
            % If image rays are not provided, obtain a regular sampling of
            % rays
            if nargin <= 1
                [pixel_i,pixel_j,microlens_k,microlens_l] = ndgrid(1:4:10,1:4:10,1:20:100,1:20:100);
                imageRays_ijkl = [pixel_i(:),pixel_j(:),microlens_k(:),microlens_l(:)]';
            end
            imageRays_ijkl = utils.enums.Classes.IMAGE_RAY().convert(imageRays_ijkl);

            % Represent image rays with homogeneous coordinates
            imageRays_ijkl = imageRays_ijkl.setHomogeneousCoordinates();
            
            % Obtain back-projection rays
            imagePlane     = lightfield.ImageRay(intrinsicMatrixInstance.imagePlaneMatrix         * imageRays_ijkl.data,true);
            microlensPlane = lightfield.ImageRay(intrinsicMatrixInstance.microlensPlaneMatrix     * imageRays_ijkl.data,true);
            mainLensPlane  = lightfield.ImageRay(intrinsicMatrixInstance.mainLensPlaneMatrix      * imageRays_ijkl.data,true);
            worldPlane     = lightfield.ImageRay(intrinsicMatrixInstance.intrinsicMatrixDirection * imageRays_ijkl.data,true);
            imagePlane     = imagePlane.removeHomogeneousCoordinates();
            microlensPlane = microlensPlane.removeHomogeneousCoordinates();
            mainLensPlane  = mainLensPlane.removeHomogeneousCoordinates();
            worldPlane     = worldPlane.removeHomogeneousCoordinates();
            
            % Obtain 3D coordinates of back-projection rays. Only the first
            % 2 components matter since the remaining correspond to
            % directions.
            backProjectionRays3D = [ [imagePlane.i    ; imagePlane.j    ] ...
                                   , [microlensPlane.i; microlensPlane.j] ...
                                   , [mainLensPlane.i ; mainLensPlane.j ] ...
                                   , [worldPlane.i    ; worldPlane.j    ] ...
                                   ; - (   intrinsicMatrixInstance.worldPlaneDistance ...
                                         + intrinsicMatrixInstance.distanceMicrolensToMainLens ...
                                         + intrinsicMatrixInstance.microlensFocalLength ) ...
                                     * ones(1,imagePlane.numberVectors) ...
                                   , - (   intrinsicMatrixInstance.worldPlaneDistance ...
                                         + intrinsicMatrixInstance.distanceMicrolensToMainLens ) ...
                                     * ones(1,microlensPlane.numberVectors) ...
                                   , - intrinsicMatrixInstance.worldPlaneDistance ...
                                     * ones(1,mainLensPlane.numberVectors) ...
                                   , zeros(1,worldPlane.numberVectors) ...
                                   ];
                               
            % Create back-projection rays class to display results
            self = camera.models.plenoptic.BackProjectionRays(backProjectionRays3D);
        end
        
        function self = BackProjectionRaysFromIntrinsicMatrixWithDepthField(intrinsicMatrixInstance, imageRays_ijkl)
            %
            % Create back-projection rays instance with adjacent rays from 
            % a plenoptic camera intrinsic matrix.
            %
            % Obtain the coordinates of the back-projection rays including
            % the back-projection rays of the pixel limits. These 
            % coordinates correspond to the image, microlens, main lens and
            % the world plane in metric units.
            %
            % INPUTS:
            %   1. intrinsicMatrixInstace - intrinsic matrix instance for
            %      plenoptic cameras. Default is a standard plenoptic
            %      camera intrinsic matrix.
            %   2. imageRays_ijkl - image space rays in pixels and
            %      microlenses indices. Each ray should be provided in 
            %      different columns.
            %
            narginchk(0,2);
            
            % If intrinsic matrix instance is not provided, obtain example
            % of intrinsic matrix
            if nargin <= 0
                intrinsicMatrixInstance = camera.models.plenoptic.standard.IntrinsicMatrix();
            end
            
            % If image rays are not provided, obtain a regular sampling of
            % rays
            if nargin <= 1
                [pixel_i,pixel_j,microlens_k,microlens_l] = ndgrid(1:4:10,1:4:10,1:20:100,1:20:100);
                imageRays_ijkl = [pixel_i(:),pixel_j(:),microlens_k(:),microlens_l(:)]';
            end
            imageRays_ijkl = utils.enums.Classes.IMAGE_RAY().convert(imageRays_ijkl);

            % Obtain more points for depth of field
            outerImageRays_ijkl = [imageRays_ijkl.i - 0.5;imageRays_ijkl.j - 0.5;imageRays_ijkl.k;imageRays_ijkl.l];
            innerImageRays_ijkl = [imageRays_ijkl.i + 0.5;imageRays_ijkl.j + 0.5;imageRays_ijkl.k;imageRays_ijkl.l];
            rays_ijkl           = [outerImageRays_ijkl,imageRays_ijkl.data,innerImageRays_ijkl];
            
            % Obtain back-projection for all rays
            backProjectionRays3D = ...
                camera.models.plenoptic.BackProjectionRays.BackProjectionRaysFromIntrinsicMatrix( ...
                        intrinsicMatrixInstance, rays_ijkl ).rays.data;
            
            % Divide back-projection into outer, middle and inner pixels
            % Obtain indices to select blocks of rays
            numberRaysPerBlock    = length(rays_ijkl);
            blockSelectionIndices = 1:imageRays_ijkl.numberVectors;
            blockSelectionIndices = [ blockSelectionIndices ...
                                    , blockSelectionIndices +     numberRaysPerBlock ...
                                    , blockSelectionIndices + 2 * numberRaysPerBlock ...
                                    , blockSelectionIndices + 3 * numberRaysPerBlock ];

            self = camera.models.plenoptic.BackProjectionRays();
            self.outerRays = backProjectionRays3D(:,blockSelectionIndices);
            self.rays      = backProjectionRays3D(:,blockSelectionIndices +     imageRays_ijkl.numberVectors);
            self.innerRays = backProjectionRays3D(:,blockSelectionIndices + 2 * imageRays_ijkl.numberVectors);
        end
    end
end
