classdef Distortion
    %DISTORTION
    %   Distortion for camera rays.
    %
    %   For plenoptic cameras, the distortion is applied to the coordinates
    %   (s,t) and (u,v) in the object space. The coordinates (u,v) are
    %   corrected assuming a scaling of the coordinates (s,t). The
    %   coordinates (s,t) are only corrected if we are using the distortion 
    %   model from Bok et al.
    %
    %   Considering the radial distortion parameters k = (k1,k2,k3,k4), the
    %   tangential distortion parameters p = (p1,p2,p3), the scaling
    %   (k5,k6,k7,k8) for s and t, and the center of distortion 
    %   c = (c_u,c_v). The undistorted coordinates m_n are given by:
    %   u_n = (1 + k1 * r^2 + k2 * r^4 + k3 * r^6 + k4 * r^8) (u_d - c_u) 
    %       + [ p1 (r^2 + 2 (u_d - c_u)^2) + 2 p2 (u_d - c_u) (v_d - c_v) ]
    %         (1 + p3 r^2) + k5 * s + k6 * s^2 + c_u
    %
    %   v_n = (1 + k1 * r^2 + k2 * r^4 + k3 * r^6 + k4 * r^8) (v_d - c_v) 
    %       + [ p2 (r^2 + 2 (v_d - c_v)^2) + 2 p1 (u_d - c_u) (v_d - c_v) ]
    %         (1 + p3 r^2) + k7 * t + k8 * t^2 + c_v
    %
    %   with r = sqrt( (u_d - c_u)^2 + (v_d - c_v)^2 ) and m_r are the
    %   distorted coordinates. The same model is considered for the (s,t)
    %   coordinates.
    %
    
    properties
        radial     = [0;0;0;0]              % Parameters for radial distortion
        tangential = [0;0;0]                % Parameters for tangential distortion
        scaling    = image.Pixel([0,0;0,0]) % Parameters for scaling (s,t)
        center     = image.Pixel([0;0])     % Distortion center
        parameters = []
        distortionModel = camera.models.plenoptic.enums.DistortionModels.DANSEREAU_ORIGINAL()
                % Distortion model
    end
    
    properties (Dependent)
        k1      % Radial distortion parameter (r^2)
        k2      % Radial distortion parameter (r^4)
        k3      % Radial distortion parameter (r^6)
        k4      % Radial distortion parameter (r^8)
        k5      % Scaling parameter for s
        k6      % Scaling parameter for s^2
        k7      % Scaling parameter for t
        k8      % Scaling parameter for t^2
        p1      % Tangential distortion parameter (coordinate u)
        p2      % Tangential distortion parameter (coordinate v)
        p3      % Tangential distortion parameter (coordinates u,v)
    end
    
    methods
        function self = Distortion(varargin)
            %
            % Create distortion instance.
            %
            % INPUTS:
            %   1. radialDistortionParameters - radial distortion 
            %   parameters.
            %   2. tangentialDistortionParameters - tangential distortion
            %   parameters.
            %   3. scalingParameters - scaling parameters for s and t.
            %   4. center - center of distortion.
            %   5. distortionModel - distortion model.
            %
            narginchk(0,5);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.radial = varargin{1};
                end
                
                if nargin >= 2
                    self.tangential = varargin{2};
                end
                
                if nargin >= 3
                    self.scaling = varargin{3};
                end
                
                if nargin >= 4
                    self.center = varargin{4};
                end
                
                if nargin >= 5
                    self.distortionModel = varargin{5};
                end
            end
        end
        
        function self = set.radial(self,newRadialDistortionParameters)
            self.radial = newRadialDistortionParameters(:);
        end
        
        function self = set.tangential(self,newTangentialDistortionParameters)
            self.tangential = newTangentialDistortionParameters(:);
        end
        
        function self = set.scaling(self,newScalingParameters)
            self.scaling = utils.enums.Classes.PIXEL().convert(newScalingParameters);
        end
        
        function self = set.center(self,newDistortionCenter)
            self.center = utils.enums.Classes.PIXEL().convert(newDistortionCenter);
        end
        
        function value = get.k1(self)
            value = self.radial(1);
        end
        
        function value = get.k2(self)
            value = self.radial(2);
        end
        
        function value = get.k3(self)
            value = self.radial(3);
        end
        
        function value = get.k4(self)
            value = self.radial(4);
        end
        
        function value = get.k5(self)
            value = self.scaling.u(1);
        end
        
        function value = get.k6(self)
            value = self.scaling.u(2);
        end
        
        function value = get.k7(self)
            value = self.scaling.v(1);
        end
        
        function value = get.k8(self)
            value = self.scaling.v(2);
        end
        
        function value = get.p1(self)
            value = self.tangential(1);
        end
        
        function value = get.p2(self)
            value = self.tangential(2);
        end
        
        function value = get.p3(self)
            value = self.tangential(3);
        end
    end
    
    methods
        function parameters = encode(self)
            %
            % Encode distortion models into a vectorized version.
            %
            parameters = self.distortionModel.encode(self);
        end
        
        function self = decode(self,parameters)
            %
            % Decode encoded parameters into a distortion model.
            %
            % INPUTS:
            %   1. parameters - encoded distortion parameters. 
            %
            narginchk(2,2);
            self = self.distortionModel.decode(parameters);
        end
        
        function undistortedObjectRays_stuv = removeDistortion(self,distortedObjectRays_stuv)
            % 
            % Remove radial and tangential distortion from distorted object 
            % rays.
            %
            % INPUTS:
            %   1. distortedObjectRays_stuv - rays in object space with
            %   distortion.
            %
            narginchk(1,2);

            % If no rays are provided, generate random rays
            if nargin <= 1
                distortedObjectRays_stuv = rand(2,50);
            end
            
            % Convert to object rays
            distortedObjectRays_stuv = utils.enums.Classes.OBJECT_RAY().convert(distortedObjectRays_stuv);
            
            % Remove distortion from directions (u,v)
            % Obtain distance to the center ray
            centerRayDistance    = image.Pixel( [distortedObjectRays_stuv.u;distortedObjectRays_stuv.v] ...
                                              - repmat(self.center.data,1,distortedObjectRays_stuv.numberVectors) ...
                                              , false );
            centerRayDistance_r2 = sum(centerRayDistance.data.^2);
            
            % Remove distortion to object rays
            radial_uv     = centerRayDistance.data ...
                         .* repmat( ( 1 + self.k1 .* centerRayDistance_r2 ...      % (r^2)
                                        + self.k2 .* centerRayDistance_r2.^2 ...   % (r^4)
                                        + self.k3 .* centerRayDistance_r2.^3 ...   % (r^6)
                                        + self.k4 .* centerRayDistance_r2.^4 ) ... % (r^8)
                                  , 2, 1) ...
                          + repmat(self.center.data,1,distortedObjectRays_stuv.numberVectors);
            tangential_uv = (        [self.p1; self.p2] .* ( repmat(centerRayDistance_r2,2,1) + 2 * centerRayDistance.data.^2 ) ...
                              + 2 .* [self.p2; self.p1] .* ( repmat(centerRayDistance.u .* centerRayDistance.v,2,1) )) ...
                         .* ( 1 + self.p3 .* repmat(centerRayDistance_r2,2,1) );
            undistortedObjectRays_uv = radial_uv + tangential_uv ...
                                     + [ self.k5 .* distortedObjectRays_stuv.s ...
                                       + self.k6 .* distortedObjectRays_stuv.s.^2 ...
                                       ; self.k7 .* distortedObjectRays_stuv.t ...
                                       + self.k8 .* distortedObjectRays_stuv.t.^2 ];
            
            % Remove distortion from positions (s,t) if we are dealing with
            % the distortion model from Bok.
            if self.distortionModel == camera.models.plenoptic.enums.DistortionModels.BOK_ORIGINAL ...
            || self.distortionModel == camera.models.plenoptic.enums.DistortionModels.BOK_DISTORTION_CENTER
                % Obtain distance to the center position
                centerPositionDistance    = image.Pixel( [distortedObjectRays_stuv.s;distortedObjectRays_stuv.t] ...
                                                       - repmat(self.center.data,1,distortedObjectRays_stuv.numberVectors) ...
                                                       , false );
                centerPositionDistance_r2 = sum(centerPositionDistance.data.^2);

                % Remove distortion to object positions
                radial_st     = centerPositionDistance.data ...
                             .* repmat( ( 1 + self.k1 .* centerPositionDistance_r2 ...      % (r^2)
                                            + self.k2 .* centerPositionDistance_r2.^2 ...   % (r^4)
                                            + self.k3 .* centerPositionDistance_r2.^3 ...   % (r^6)
                                            + self.k4 .* centerPositionDistance_r2.^4 ) ... % (r^8)
                                      , 2, 1) ...
                              + repmat(self.center.data,1,distortedObjectRays_stuv.numberVectors);
                tangential_st = (        [self.p1; self.p2] .* ( repmat(centerPositionDistance_r2,2,1) + 2 * centerPositionDistance.data.^2 ) ...
                                  + 2 .* [self.p2; self.p1] .* ( repmat(centerPositionDistance.u .* centerPositionDistance.v,2,1) )) ...
                             .* ( 1 + self.p3 .* repmat(centerPositionDistance_r2,2,1) );
                undistortedObjectRays_st = radial_st + tangential_st;

            % Otherwise, maintain the positions unchanged.
            else
                undistortedObjectRays_st = [distortedObjectRays_stuv.s;distortedObjectRays_stuv.t];
            end
            
            % Obtain undistorted object ray
            undistortedObjectRays_stuv = lightfield.ObjectRay( [ undistortedObjectRays_st ...
                                                               ; undistortedObjectRays_uv ] ...
                                                             , false );
        end
        
        function distortedObjectRays_stuv = addDistortion( self ...
                                                         , undistortedObjectRays_stuv ...
                                                         , numberIterations )
            % 
            % Add radial and tangential distortion to undistorted rays.
            %
            % INPUTS:
            %   1. undistortedObjectRays_stuv - rays in object space
            %   without distortion.
            %   2. numberIterations    - number iterations consider to 
            %   invert the distortion model. This model is nonlinear. 
            %   Default is 5.
            %
            narginchk(1,3);
            
            % If no rays are provided, generate random rays
            if nargin <= 1
                undistortedObjectRays_stuv = rand(2,50);
            end
            
            % If no number of iterations is provided, assume default value
            if nargin <= 2
                numberIterations = 5;
            end
            
            % If rays are provided in rows instead of columns, transpose
            % the rays
            undistortedObjectRays_stuv = utils.enums.Classes.OBJECT_RAY().convert(undistortedObjectRays_stuv);
            
            % The inverse model for the complete distortion model is not
            % implemented. The distortion model is only completely
            % characterized when there is no tangential distortion
            if sum(abs(self.tangential)) > 0
                warning( 'addDistortion:tangentialDistortionNotCharacterized' ...
                       , [ 'Inverse distortion model for tangential parameters is not ' ...
                         , 'implemented. The distortion will be added considering that ' ...
                         , 'these parameters are zero.' ]);
            end

            % Add distortion to coordinates (u,v)
            % Obtain distance to the center ray
            centerRayDistance = [undistortedObjectRays_stuv.u;undistortedObjectRays_stuv.v] ...
                              - repmat(self.center.data,1,undistortedObjectRays_stuv.numberVectors) ...
                              - [ self.k5 .* undistortedObjectRays_stuv.s ...
                                + self.k6 .* undistortedObjectRays_stuv.s.^2 ...
                                ; self.k7 .* undistortedObjectRays_stuv.t ...
                                + self.k8 .* undistortedObjectRays_stuv.t.^2 ];
            
            % Obtain distorted rays distance iteratively
            distortedCenterRayDistance = centerRayDistance;
            for iIteration = 1:numberIterations
                centerRayDistance_r2 = sum(distortedCenterRayDistance.^2);

                % Add distortion to rays
                distortedCenterRayDistance = centerRayDistance ...
                                          ./ repmat( ( 1 + self.k1 .* centerRayDistance_r2 ...      % (r^2)
                                                         + self.k2 .* centerRayDistance_r2.^2 ...   % (r^4)
                                                         + self.k3 .* centerRayDistance_r2.^3 ...   % (r^4)
                                                         + self.k4 .* centerRayDistance_r2.^4 ) ... % (r^8)
                                                   , 2, 1 );
            end
        
            % Add center of distortion and scaling to obtain distorted rays
            distortedObjectRays_uv = distortedCenterRayDistance ...
                                   + repmat(self.center.data,1,undistortedObjectRays_stuv.numberVectors);

            % Add distortion to positions (s,t) if we are dealing with the
            % distortion model from Bok.
            if self.distortionModel == camera.models.plenoptic.enums.DistortionModels.BOK_ORIGINAL ...
            || self.distortionModel == camera.models.plenoptic.enums.DistortionModels.BOK_DISTORTION_CENTER
                % Obtain distance to the center position
                centerPositionDistance = [undistortedObjectRays_stuv.s;undistortedObjectRays_stuv.t] ...
                                       - repmat(self.center.data,1,undistortedObjectRays_stuv.numberVectors);

                % Obtain undistorted rays distance iteratively
                distortedCenterPositionDistance = centerPositionDistance;
                for iIteration = 1:numberIterations
                    centerPositionDistance_r2 = sum(distortedCenterPositionDistance.^2);

                    % Add distortion to positions
                    distortedCenterPositionDistance = centerPositionDistance ...
                                                   ./ repmat( ( 1 + self.k1 .* centerPositionDistance_r2 ...      % (r^2)
                                                                  + self.k2 .* centerPositionDistance_r2.^2 ...   % (r^4)
                                                                  + self.k3 .* centerPositionDistance_r2.^3 ...   % (r^4)
                                                                  + self.k4 .* centerPositionDistance_r2.^4 ) ... % (r^8)
                                                            , 2, 1 );
                end
        
                % Add center of distortion to obtain distorted positions
                distortedObjectRays_st = distortedCenterPositionDistance ...
                                       + repmat(self.center.data,1,undistortedObjectRays_stuv.numberVectors);

            % Otherwise, maintain the positions unchanged.
            else
                distortedObjectRays_st = [undistortedObjectRays_stuv.s;undistortedObjectRays_stuv.t];
            end
            
            % Obtain distorted object ray
            distortedObjectRays_stuv = lightfield.ObjectRay( [ distortedObjectRays_st ...
                                                             ; distortedObjectRays_uv ] ...
                                                           , false );
        end
    end
end
