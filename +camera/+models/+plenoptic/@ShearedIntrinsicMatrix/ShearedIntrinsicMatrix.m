classdef ShearedIntrinsicMatrix
    %SHEAREDINTRINSICMATRIX
    %   Sheared intrinsic matrix for plenoptic camera.
    %   Since the sheared intrinsic matrix works with the coordinates 
    %   (i,j,k,l) and these coordinates define the lightfield using the 
    %   two-plane parameterization with two points, we can use any 
    %   intrinsic matrix.    
    %
    %   Remember that the indices (i,j,k,l) are one-based indices.
    
    properties
        data = camera.models.plenoptic.standard.IntrinsicMatrix().intrinsicMatrixDirection
                % Intrinsic matrix to convert rays from image to object space
    end
    
    methods 
        function self = ShearedIntrinsicMatrix(varargin)
            %
            % Create sheared intrinsic matrix instance.
            % 
            % INPUTS: 
            %   1. intrinsicMatrix  - intrinsic matrix data.
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
    end
    
    methods
        function shearedMatrix = obtainIntrisicMatrixChangingWorldFocalPlane_kl( self ...
                                                                               , disparity_dk_di ...
                                                                               , referenceViewpoint )
            %
            % Obtain intrinsic matrix by changing the world focal plane.
            %
            % INPUTS:
            %   1. disparity_dk_di - disparity on viewpoints.
            %   2. referenceViewpoint - reference viewpoint for shearing
            %   lightfield.
            %
            narginchk(1,3);
            
            % If disparity is not provided, assume default value.
            if nargin <= 1
                disparity_dk_di = 1;
            end
            
            % If reference viewpoint is not provided, assume default value.
            if nargin <= 1
                referenceViewpoint = [5,5,nan,nan];
            end

            % Transform reference viewpoint to image ray instance
            referenceViewpoint = utils.enums.Classes.IMAGE_RAY().convert(referenceViewpoint);
            
            % Obtain matrix for resampling the acquired lightfield
            resamplingMatrix = [ 1, 0, 0, 0, 0 ...
                               ; 0, 1, 0, 0, 0 ...
                               ; disparity_dk_di, 0, 1, 0, -referenceViewpoint.i * disparity_dk_di ...
                               ; 0, disparity_dk_di, 0, 1, -referenceViewpoint.j * disparity_dk_di ...
                               ; 0, 0, 0, 0, 1];
        
            % Obtain intrinsic matrix after shearing
            shearedMatrix = camera.models.plenoptic.IntrinsicMatrixDirection( ...
                                    self.data * resamplingMatrix );
        end
        
        function shearedMatrix = obtainIntrisicMatrixChangingWorldFocalPlane_ij( self ...
                                                                               , disparity_dk_di ...
                                                                               , referenceMicrolens )
            %
            % Obtain intrinsic matrix by changing the world focal plane.
            %
            % INPUTS:
            %   1. disparity_dk_di - disparity on viewpoints.
            %   2. referenceMicrolens - reference microlens for shearing
            %   lightfield.
            %
            narginchk(1,3);
            
            % If disparity is not provided, assume default value.
            if nargin <= 1
                disparity_dk_di = 1;
            end
            
            % If reference microlens is not provided, assume default value.
            if nargin <= 1
                referenceMicrolens = [nan,nan,180,180];
            end

            % Transform reference microlens to image ray instance
            referenceMicrolens = utils.enums.Classes.IMAGE_RAY().convert(referenceMicrolens);
            
            % Obtain matrix for resampling the acquired lightfield
            resamplingMatrix = [ 1, 0, 1/disparity_dk_di, 0, -referenceMicrolens.k / disparity_dk_di ...
                               ; 0, 1, 0, 1/disparity_dk_di, -referenceMicrolens.l / disparity_dk_di ...
                               ; 0, 0, 1, 0, 0 ...
                               ; 0, 0, 0, 1, 0 ...
                               ; 0, 0, 0, 0, 1];
        
            % Obtain intrinsic matrix after shearing
            shearedMatrix = camera.models.plenoptic.IntrinsicMatrixDirection( ...
                                    self.data * resamplingMatrix );
        end
    end
end
