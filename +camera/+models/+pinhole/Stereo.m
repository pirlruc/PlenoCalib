classdef Stereo
    %STEREO
    %   Utility for stereo reconstruction using pinhole camera models.
    
    properties
        camera_A = camera.models.pinhole.Pinhole() 	% Pinhole camera instance A
        camera_B = camera.models.pinhole.Pinhole() 	% Pinhole camera instance B
    end
    
    methods
        function self = Stereo(varargin)
            %
            % Create stereo reconstruction instance.
            %
            % INPUTS: 
            %   1. pinholeCamera_A - pinhole camera model.
            %   2. pinholeCamera_B - pinhole camera model.
            %
            narginchk(0,2);
            
            if ~isempty(varargin)
                if nargin >= 2
                    self.camera_B = varargin{2};
                end
                
                if nargin >= 1
                    self.camera_A = varargin{1};
                end
            end
        end
        
        function self = set.camera_A(self,newPinholeCamera_A)
            self.camera_A = newPinholeCamera_A;
        end
        
        function self = set.camera_B(self,newPinholeCamera_B)
            self.camera_B = newPinholeCamera_B;
        end
    end
    
    methods
        function [worldPoint,rmse,sse,residuals] = reconstruct(self,pixels_A, pixels_B)
            %
            % Obtain world points by stereo reconstruction from image
            % correspondences.
            %
            % INPUTS:
            %   1. pixels_A - image pixel correspondences for camera A.
            %   2. pixels_B - image pixel correspondences for camera B.
            %
            narginchk(3,3);
            
            % If pixels are provided in rows instead of columns, transpose
            % the pixels
            pixels_A = utils.enums.Classes.PIXEL().convert(pixels_A);
            pixels_B = utils.enums.Classes.PIXEL().convert(pixels_B);
            
            % Select entries of projection matrices
            p_uv_lambda_A     = self.camera_A.projectionMatrix(1:2,4);
            p_uv_lambda_B     = self.camera_B.projectionMatrix(1:2,4);
            p_lambda_lambda_A = self.camera_A.projectionMatrix(3,4);
            p_lambda_lambda_B = self.camera_B.projectionMatrix(3,4);
            p_u_A             = self.camera_A.projectionMatrix(1,1:3);
            p_u_B             = self.camera_B.projectionMatrix(1,1:3);
            p_v_A             = self.camera_A.projectionMatrix(2,1:3);
            p_v_B             = self.camera_B.projectionMatrix(2,1:3);
            p_lambda_A        = self.camera_A.projectionMatrix(3,1:3);
            p_lambda_B        = self.camera_B.projectionMatrix(3,1:3);
            
            % Construct vector and matrix for least squares
            vector = [ p_uv_lambda_A - pixels_A.data .* p_lambda_lambda_A ...
                     ; p_uv_lambda_B - pixels_B.data .* p_lambda_lambda_B ];
            matrix = [ repmat(pixels_A.data,1,length(p_lambda_A)) .* repmat(p_lambda_A,pixels_A.numberComponents,1) ...
                       - [p_u_A;p_v_A]  ... 
                     ; repmat(pixels_B.data,1,length(p_lambda_B)) .* repmat(p_lambda_B,pixels_B.numberComponents,1) ...
                       - [p_u_B;p_v_B] ];

            % Solve least squares problem
            worldPoint = math.Point(mldivide(matrix,vector),false);
            
            % Obtain residuals and errors for least squares problem
            if nargout > 1
                residuals = abs(matrix * worldPoint.data - vector);
                sse  = sum(residuals.^2);               % Sum of squared errors
                rmse = sqrt(sse / length(residuals));   % Root mean squared error
            end
        end
    end
    
    methods (Static)
        function worldPoints = obtainWorldPointsForDepth(depthValues,numberPointsPerDepth)
            %
            % Obtain world points for different depth values.
            %
            % INPUTS:
            %   1. depthValues          - depth values to be considered.
            %   2. numberPointsPerDepth - number of points generated per
            %      depth.
            %
            narginchk(0,2);
            
            if nargin <= 1
                numberPointsPerDepth = 20;
            end
            
            if nargin <= 0
                depthValues = 0.01:0.05:1.0;
            end
            
            % Obtain world points
            worldPoints = [];
            for depth = depthValues
                % Generate random points for depth
                xy = depth .* (-1 + 2.* rand(numberPointsPerDepth,2));
                z  = depth .* ones(numberPointsPerDepth,1);
                
                % Concatenate points
                worldPoints = cat(1,worldPoints,[xy,z]);
            end
            
            % Ensure that world points are given in the correct order in
            % terms of dimensions.
            worldPoints = math.Point(worldPoints,false);
        end
    end
end

