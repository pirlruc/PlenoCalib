classdef IntrinsicMatrixPoints
    %INTRINSICMATRIXPOINTS
    %   Intrinsic matrix for plenoptic cameras. The intrinsic matrix
    %   represents the lightfield in the object space with a two-plane
    %   parameterization with two points.
    
    properties
        data = camera.models.plenoptic.standard.IntrinsicMatrix().intrinsicMatrixPoints
                % Intrinsic matrix to convert rays from image to object space.
        distancePlanes_st_uv = 1    % Distance between planes (s,t) and (u,v)
    end
    
    properties (Dependent)
        intrinsicMatrixDirection
            % Back-project image rays to the world plane in metric units.
            % Applying this intrinsic matrix generates a two-plane 
            % parameterization with one point and one direction.
        distance_st_ij  % Distance to new plane (s,t) considering the viewpoint images (i,j) 
        distance_uv_ij  % Distance to new plane (u,v) considering the viewpoint images (i,j) 
        distance_st_kl  % Distance to new plane (s,t) considering the microlens images (k,l) 
        distance_uv_kl  % Distance to new plane (u,v) considering the microlens images (k,l) 
        reparameterizedIntrinsicMatrix_ij
            % Reparameterized intrinsic matrix to have the (s,t) plane
            % coincident with the plane with the centers of projection of
            % the viewpoint images (i,j). Additionally, the (u,v) plane is
            % coincident with the plane with the centers of projection of
            % the microlenses images (k,l).
        reparameterizedIntrinsicMatrix_kl
            % Reparameterized intrinsic matrix to have the (s,t) plane
            % coincident with the plane with the centers of projection of
            % the microlens images (k,l). Additionally, the (u,v) plane is
            % coincident with the plane with the centers of projection of
            % the viewpoint images (i,j).
    end
    
    methods
        function self = IntrinsicMatrixPoints(varargin)
            %
            % Create symbolic intrinsc matrix points instance.
            %
            % INPUTS:
            %   1. intrinsicMatrixPoints - intrinsic matrix that represents
            %   the lightfield in the object space with 2 points.
            %   2. distancePlanes_st_uv  - distance between the planes
            %   (s,t) and (u,v).
            %
            narginchk(0,2);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.data = varargin{1};
                end
                
                if nargin >= 2
                    self.distancePlanes_st_uv = varargin{2};
                end
            end
        end
        
        function self = set.data(self,newIntrinsicMatrix)
            self.data = newIntrinsicMatrix;
        end
        
        function intrinsic = get.intrinsicMatrixDirection(self)
            % Obtain matrix to convert the lightfield parameterization with
            % 2 points for a lightfield with one point and one direction.
            pointsToDirection = [ 1, 0, 0, 0, 0 ...
                                ; 0, 1, 0, 0, 0 ...
                                ; -1/self.distancePlanes_st_uv, 0, 1/self.distancePlanes_st_uv, 0, 0 ...
                                ; 0, -1/self.distancePlanes_st_uv, 0, 1/self.distancePlanes_st_uv, 0 ...
                                ; 0, 0, 0, 0, 1 ];

            % Obtain new intrinsic matrix
            intrinsic = camera.models.plenoptic.IntrinsicMatrixDirection( ...
                                pointsToDirection * self.data );
        end
        
        function distance = get.distance_st_ij(self)
            % Obtain the distance to the projection centers of the
            % viewpoint images. 
            
            % Create plenoptic camera from intrinsic matrix.
            intrinsic = camera.models.plenoptic.standard.Camera(self.data);

            % Two distances are given since we do not know if the intrinsic
            % matrix defines pinholes. These are the distances for each
            % coordinate pair (i,k) and (j,l).
            distance.ik = - ( intrinsic.h_sk * self.distancePlanes_st_uv ) ...
                          / ( intrinsic.h_uk - intrinsic.h_sk );
            distance.jl = - ( intrinsic.h_tl * self.distancePlanes_st_uv ) ...
                          / ( intrinsic.h_vl - intrinsic.h_tl );
        end
        
        function distance = get.distance_uv_ij(self)
            % Obtain the distance to new plane used to define the rays 
            % (u,v).
            
            % Create plenoptic camera from intrinsic matrix.
            intrinsic = camera.models.plenoptic.standard.Camera(self.data);

            % Two distances are given since we do not know if the intrinsic
            % matrix defines pinholes. These are the distances for each
            % coordinate pair (i,k) and (j,l).
            distance.ik = - ( intrinsic.h_ui * self.distancePlanes_st_uv ) ...
                          / ( intrinsic.h_ui - intrinsic.h_si );
            distance.jl = - ( intrinsic.h_vj * self.distancePlanes_st_uv ) ...
                          / ( intrinsic.h_vj - intrinsic.h_tj );
        end
        
        function distance = get.distance_st_kl(self)
            % Obtain the distance to the projection centers of the
            % microlens images. 
            
            % Create plenoptic camera from intrinsic matrix.
            intrinsic = camera.models.plenoptic.standard.Camera(self.data);

            % Two distances are given since we do not know if the intrinsic
            % matrix defines pinholes. These are the distances for each
            % coordinate pair (i,k) and (j,l).
            distance.ik = - ( intrinsic.h_si * self.distancePlanes_st_uv ) ...
                          / ( intrinsic.h_ui - intrinsic.h_si );
            distance.jl = - ( intrinsic.h_tj * self.distancePlanes_st_uv ) ...
                          / ( intrinsic.h_vj - intrinsic.h_tj );
        end
        
        function distance = get.distance_uv_kl(self)
            % Obtain the distance to new plane used to define the rays 
            % (u,v).
            
            % Create plenoptic camera from intrinsic matrix.
            intrinsic = camera.models.plenoptic.standard.Camera(self.data);

            % Two distances are given since we do not know if the intrinsic
            % matrix defines pinholes. These are the distances for each
            % coordinate pair (i,k) and (j,l).
            distance.ik = - ( intrinsic.h_uk * self.distancePlanes_st_uv ) ...
                          / ( intrinsic.h_uk - intrinsic.h_sk );
            distance.jl = - ( intrinsic.h_vl * self.distancePlanes_st_uv ) ...
                          / ( intrinsic.h_vl - intrinsic.h_tl );
        end
        
        function intrinsicMatrix = get.reparameterizedIntrinsicMatrix_ij(self)
            %
            % Obtain reparameterized intrinsic matrix to have the (s,t) 
            % plane coincident with the plane of the centers of projection 
            % of the viewpoint images (i,j). Additionally, the (u,v) plane 
            % is coincident with the plane with the centers of projection 
            % of the microlenses images (k,l).
            %

            % Obtain distances to new plane (s',t')
            distance_st = self.distance_st_ij;
                
            % Obtain distances to new plane (u',v')
            distance_uv = self.distance_uv_ij;
            
            % Obtain the matrix to be used for moving the plane (s,t) to 
            % (s',t') and the plane (u,v) to (u',v').
            movePlane_st   = [ 1, 0, distance_st.ik, 0, 0 ...
                             ; 0, 1, 0, distance_st.jl, 0 ...
                             ; 0, 0, 1, 0, 0 ...
                             ; 0, 0, 0, 1, 0 ...
                             ; 0, 0, 0, 0, 1 ];
            movePlane_uv   = [ 1, 0, 0, 0, 0 ...
                             ; 0, 1, 0, 0, 0 ...
                             ; 1, 0, self.distancePlanes_st_uv - distance_st.ik + distance_uv.ik, 0, 0 ...
                             ; 0, 1, 0, self.distancePlanes_st_uv - distance_st.jl + distance_uv.jl, 0 ...
                             ; 0, 0, 0, 0, 1 ];
            movePlane_stuv = movePlane_uv * movePlane_st;
            
            % Obtain the intrinsic matrix that results from moving the
            % plane (s,t) and (u,v)
            intrinsicMatrix = camera.models.plenoptic.IntrinsicMatrixPoints( ...
                                    movePlane_stuv * self.intrinsicMatrixDirection.data );
                                
            % Set to zero entries that are smaller than precision
            intrinsicMatrix.data(abs(intrinsicMatrix.data) <= eps) = 0;
        end
        
        function intrinsicMatrix = get.reparameterizedIntrinsicMatrix_kl(self)
            %
            % Obtain reparameterized intrinsic matrix to have the (s,t) 
            % plane coincident with the plane with the centers of 
            % projection of the microlens images (k,l). Additionally, the 
            % (u,v) plane is coincident with the plane with the centers of 
            % projection of the viewpoint images (i,j).
            %
            
            % Obtain distances to new plane (s',t')
            distance_st = self.distance_st_kl;
                
            % Obtain distances to new plane (u',v')
            distance_uv = self.distance_uv_kl;
            
            % Obtain the matrix to be used for moving the plane (s,t) to 
            % (s',t') and the plane (u,v) to (u',v').
            movePlane_st   = [ 1, 0, distance_st.ik, 0, 0 ...
                             ; 0, 1, 0, distance_st.jl, 0 ...
                             ; 0, 0, 1, 0, 0 ...
                             ; 0, 0, 0, 1, 0 ...
                             ; 0, 0, 0, 0, 1 ];
            movePlane_uv   = [ 1, 0, 0, 0, 0 ...
                             ; 0, 1, 0, 0, 0 ...
                             ; 1, 0, self.distancePlanes_st_uv - distance_st.ik + distance_uv.ik, 0, 0 ...
                             ; 0, 1, 0, self.distancePlanes_st_uv - distance_st.jl + distance_uv.jl, 0 ...
                             ; 0, 0, 0, 0, 1 ];
            movePlane_stuv = movePlane_uv * movePlane_st;
            
            % Obtain the intrinsic matrix that results from moving the
            % plane (s,t) and (u,v)
            intrinsicMatrix = camera.models.plenoptic.IntrinsicMatrixPoints( ...
                                    movePlane_stuv * self.intrinsicMatrixDirection.data );

            % Set to zero entries that are smaller than precision
            intrinsicMatrix.data(abs(intrinsicMatrix.data) <= eps) = 0;
        end
    end
    
    methods (Static)
        function self = IntrinsicMatrixPointsFromPinholeIntrinsicMatrix( ...
                                    pinholeIntrinsicMatrix ...
                                  , viewpointDistance ...
                                  , lightfieldSize )
            %
            % Create intrinsic matrix considering a lightfield obtained
            % from a given camera defined with a specific pinhole intrinsic
            % matrix.
            %
            narginchk(0,3);
            
            % Obtain default intrinsic matrix
            if nargin <= 0
                pinholeIntrinsicMatrix = camera.models.pinhole.IntrinsicMatrix();
            end
            
            % Obtain default viewpoint distance
            if nargin <= 1
                viewpointDistance = [5,5,0] .* 1e-3;
            end
            
            % Obtain default lightfield size
            if nargin <= 2
                lightfieldSize = lightfield.LightfieldSize([3,10,10,380,380]);
            end
            
            % Convert input data to instances
            viewpointDistance = utils.enums.Classes.POINT().convert(viewpointDistance);
            lightfieldSize    = utils.enums.Classes.LIGHTFIELD_SIZE().convert(lightfieldSize);
            
            % Obtain offsets for each axis of the viewpoint projection
            % centers plane.
            s0_x = -viewpointDistance.x .* ( lightfieldSize.numberPixels_i + 1 ) ./ 2;
            s0_y = -viewpointDistance.y .* ( lightfieldSize.numberPixels_j + 1 ) ./ 2;
            
            % Obtain pinhole intrinsic matrix
            intrinsic = pinholeIntrinsicMatrix.intrinsicMatrix;
            
            % Obtain intrinsic matrix for plenoptic camera considering a
            % parameterization using two points. In this parameterization,
            % the (s,t)-plane corresponds to the plane with the centers of
            % projection and the (u,v)-plane corresponds to the image
            % planes.
            self = camera.models.plenoptic.IntrinsicMatrixPoints( ...
                                    [ viewpointDistance.x, 0, 0, 0, s0_x ...
                                    ; 0, viewpointDistance.y, 0, 0, s0_y ...
                                    ; zeros(3,2), inv(intrinsic) ] ...
                                  , pinholeIntrinsicMatrix.imagePlaneDistance );
        end
    end
end
