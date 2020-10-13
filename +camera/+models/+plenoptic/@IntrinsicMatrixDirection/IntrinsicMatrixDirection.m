classdef IntrinsicMatrixDirection
    %INTRINSICMATRIXDIRECTION
    %   Intrinsic matrix for plenoptic cameras. The intrinsic matrix 
    %   represents the lightfield in the object space with a two-plane
    %   parameterization with one point and one direction.
    
    properties
        data = camera.models.plenoptic.standard.IntrinsicMatrix().intrinsicMatrixDirection
                % Intrinsic matrix to convert rays from image to object space.
    end
    
    properties (Dependent)
        distance_st_ij  % Distance to new plane (s,t) considering the viewpoint images (i,j) 
        distance_st_kl  % Distance to new plane (s,t) considering the microlens images (k,l) 
        reparameterizedIntrinsicMatrix_ij
            % Reparameterized intrinsic matrix to have the (s,t) plane
            % coincident with the plane with the centers of projection of
            % the viewpoint images (i,j).
        reparameterizedIntrinsicMatrix_kl
            % Reparameterized intrinsic matrix to have the (s,t) plane
            % coincident with the plane with the centers of projection of
            % the microlens images (k,l).
        baselineCorrectedIntrinsicMatrix
            % Correct intrinsic matrix to have positive baseline distances
    end
    
    methods
        function self = IntrinsicMatrixDirection(varargin)
            %
            % Create symbolic intrinsc matrix direction instance.
            %
            % INPUTS:
            %   1. intrinsicMatrixDirection - intrinsic matrix that
            %   represents the lightfield in the object space with one
            %   point and one direction.
            %
            narginchk(0,1);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
        
        function self = set.data(self,newIntrinsicMatrix)
            self.data = newIntrinsicMatrix;
        end
        
        function distance = get.distance_st_ij(self)
            % Obtain the distance to the projection centers of the
            % viewpoint images. 
            
            % Create plenoptic camera from intrinsic matrix.
            intrinsic = camera.models.plenoptic.standard.Camera(self.data);

            % Two distances are given since we do not know if the intrinsic
            % matrix defines pinholes. These are the distances for each
            % coordinate pair (i,k) and (j,l).
            distance.ik = - intrinsic.h_sk / intrinsic.h_uk;
            distance.jl = - intrinsic.h_tl / intrinsic.h_vl;
        end
        
        function distance = get.distance_st_kl(self)
            % Obtain the distance to the projection centers of the
            % microlens images. 
            
            % Create plenoptic camera from intrinsic matrix.
            intrinsic = camera.models.plenoptic.standard.Camera(self.data);

            % Two distances are given since we do not know if the intrinsic
            % matrix defines pinholes. These are the distances for each
            % coordinate pair (i,k) and (j,l).
            distance.ik = - intrinsic.h_si / intrinsic.h_ui;
            distance.jl = - intrinsic.h_tj / intrinsic.h_vj;
        end
        
        function intrinsicMatrix = get.reparameterizedIntrinsicMatrix_ij(self)
            %
            % Obtain reparameterized intrinsic matrix to have the (s,t) 
            % plane coincident with the plane of the centers of projection 
            % of the viewpoint images (i,j).
            %

            % Obtain distances for viewpoint projection centers
            distance_st = self.distance_st_ij;
            
            % Obtain the matrix to be used for moving the plane (s,t) to 
            % (s',t')
            movePlane_st = [ 1, 0, distance_st.ik, 0, 0 ...
                           ; 0, 1, 0, distance_st.jl, 0 ...
                           ; 0, 0, 1, 0, 0 ...
                           ; 0, 0, 0, 1, 0 ...
                           ; 0, 0, 0, 0, 1 ];
    
            % Obtain the intrinsic matrix
            intrinsicMatrix = camera.models.plenoptic.IntrinsicMatrixDirection( ...
                                    movePlane_st * self.data );
                                
            % Set to zero entries that are smaller than precision
            intrinsicMatrix.data(abs(intrinsicMatrix.data) <= eps) = 0;
        end
                
        function intrinsicMatrix = get.reparameterizedIntrinsicMatrix_kl(self)
            %
            % Obtain reparameterized intrinsic matrix to have the (s,t) 
            % plane coincident with the plane of the centers of projection
            % of the microlens images (k,l).            
            %
            
            % Obtain distances for microlens projection centers
            distance_st = self.distance_st_kl;
            
            % Obtain the matrix to be used for moving the plane (s,t) to 
            % (s',t')
            movePlane_st = [ 1, 0, distance_st.ik, 0, 0 ...
                           ; 0, 1, 0, distance_st.jl, 0 ...
                           ; 0, 0, 1, 0, 0 ...
                           ; 0, 0, 0, 1, 0 ...
                           ; 0, 0, 0, 0, 1 ];
    
            % Obtain the intrinsic matrix
            intrinsicMatrix = camera.models.plenoptic.IntrinsicMatrixDirection( ...
                                    movePlane_st * self.data );
                                
            % Set to zero entries that are smaller than precision
            intrinsicMatrix.data(abs(intrinsicMatrix.data) <= eps) = 0;
        end
        
        function intrinsicMatrix = get.baselineCorrectedIntrinsicMatrix(self)
            % Obtain corrected intrinsic matrix to have positive baseline
            % distances for viewpoints and microlenses. This corresponds to
            % having the symmetric intrinsic matrix H.
            
            % Obtain intrinsic matrix reparameterized for viewpoints
            viewpointIntrinsicMatrix = self.reparameterizedIntrinsicMatrix_ij.data;
            
            % Intialize intrinsic matrix
            intrinsicMatrix = self;
            
            % Validate if entries corresponding to baselines are positive
            if viewpointIntrinsicMatrix(1,1) < 0 || viewpointIntrinsicMatrix(2,2) < 0
                intrinsicMatrix.data = -intrinsicMatrix.data;
                intrinsicMatrix.data(end,end) = 1;
            end
        end
    end
    
    methods 
        function intrinsic = obtainRecenteredIntrinsicMatrix(self,cameraCoordinateOriginRay)
            % 
            % Obtain recentered intrinsic matrix. This gives an intrinsic
            % matrix whose origin of the coordinate system (x,y)
            % corresponds to the central ray of the lightfield.
            %
            % INPUTS:
            %   1. cameraCoordinateOriginRay - ray coordinates for the
            %   camera coordinate system.
            %
            narginchk(2,2);
            
            cameraCoordinateOriginRay = utils.enums.Classes.IMAGE_RAY().convert(cameraCoordinateOriginRay);
            cameraCoordinateOriginRay = cameraCoordinateOriginRay.setHomogeneousCoordinates();
            
            % Obtain recentered intrinsic matrix
            intrinsicMatrix = self.data;
            intrinsicMatrix(1:4,end) = 0;
            
            % Obtain camera coordinate system origin
            cameraOrigin = lightfield.ImageRay(intrinsicMatrix * cameraCoordinateOriginRay.data,true);
            cameraOrigin = cameraOrigin.removeHomogeneousCoordinates();
            intrinsicMatrix(1:4,end) = -cameraOrigin.data;
            
            intrinsic = camera.models.plenoptic.IntrinsicMatrixDirection( intrinsicMatrix );
        end
        
        function intrinsic = obtainIntrinsicMatrixPoints(self,distancePlanes_st_uv)
            %
            % Obtain intrinsic matrix points considering a distance between
            % the planes (s,t) and (u,v). This allows to map a lightfield
            % in the image sensor to a lightfield in the object space
            % parameterized with 2 points.
            %
            % INPUTS:
            %   1. distancePlanes_st_uv - distance between planes (s,t) and
            %   (u,v).
            %
            narginchk(1,2);
            
            if nargin <= 1
                distancePlanes_st_uv = 1;
            end
                
            % Obtain matrix to propagate directions to points
            worldTo2Points = [ 1, 0,                    0,                    0, 0 ...
                             ; 0, 1,                    0,                    0, 0 ...
                             ; 1, 0, distancePlanes_st_uv,                    0, 0 ...
                             ; 0, 1,                    0, distancePlanes_st_uv, 0 ...
                             ; 0, 0,                    0,                    0, 1 ];
            
            % Obtain intrinsic matrix parameterizing the lightfield with 2
            % points.
            intrinsic = camera.models.plenoptic.IntrinsicMatrixPoints( ...
                            worldTo2Points * self.data );
        end
    end
end

