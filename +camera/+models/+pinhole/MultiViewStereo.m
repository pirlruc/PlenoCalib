classdef MultiViewStereo
    %MULTIVIEWSTEREO
    %   Utility for multi-view stereo reconstruction using pinhole camera 
    %   models.
    
    properties
        cameras = []        % Pinhole camera collection
    end
    
    properties (Dependent)
        numberCameras
    end
    
    methods
        function self = MultiViewStereo(varargin)
            %
            % Create multi-view stereo reconstruction instance.
            %
            % INPUTS: 
            %   1. pinholeCameras - pinhole camera collection.
            %
            narginchk(0,1);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.cameras = varargin{1};
                end
            end
        end
        
        function self = set.cameras(self,newPinholeCameras)
            self.cameras = newPinholeCameras;
        end
        
        function number = get.numberCameras(self)
            number = length(self.cameras);
        end
    end
    
    methods
        function [worldPoint,rmse,sse,residuals] = reconstruct(self,pixels)
            %
            % Obtain world points by multi-view stereo reconstruction from
            % image correspondences.
            %
            % INPUTS:
            %   1. pixels - image pixel correspondences for each camera.
            %      The pixels should be provided in a cell in which each
            %      element corresponds to the number of the camera.
            %
            narginchk(2,2);
            
            % We need at least two cameras defined
            if self.numberCameras <= 1
                error = MException( 'MultiViewStereo:reconstruct:noCameras' ...
                                  , 'Define at least two cameras to reconstruct...' );
                error.throw();
            end
            
            % We need at least two sets of pixels
            if length(pixels) <= 1
                error = MException( 'MultiViewStereo:reconstruct:noPixels' ...
                                  , 'Define at least two sets of pixels to reconstruct...' );
                error.throw();
            end
            
            % Obtain vector and matrix for least squares. Consider only
            % cameras for which we have pixels
            vector = [];
            matrix = [];
            % Stop if there are no more sets of pixels
            for iCamera = 1:min(self.numberCameras,length(pixels)) 
                % Obtain image correspondence pixel for camera and
                % corresponding projection projection matrix
                cameraPixel      = utils.enums.Classes.PIXEL().convert(pixels{iCamera});
                projectionMatrix = self.cameras(iCamera).projectionMatrix;
                
                % If camera has no correspondences skip this camera.
                if cameraPixel.numberVectors == 0
                    continue
                end
                
                % Select entries of projection matrix for camera
                p_uv_lambda     = projectionMatrix(1:2,4);
                p_lambda_lambda = projectionMatrix(3,4);
                p_u             = projectionMatrix(1,1:3);
                p_v             = projectionMatrix(2,1:3);
                p_lambda        = projectionMatrix(3,1:3);
                
                % Construct vector and matrix for least squares
                vector = cat(1,vector, p_uv_lambda - cameraPixel.data .* p_lambda_lambda);
                matrix = cat(1,matrix,    repmat(cameraPixel.data,1,length(p_lambda)) ...
                                       .* repmat(p_lambda,cameraPixel.numberComponents,1) ...
                                        - [p_u;p_v] );
            end
            
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
end

