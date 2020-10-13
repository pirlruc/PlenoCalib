classdef PointReconstruction
    %POINTRECONSTRUCTION
    %   Utility to evaluate reconstruction for lytro camera.
    
    properties
        plenoptic = camera.models.plenoptic.standard.Camera()
                                                % Plenoptic camera instance
        correspondences        = []             % Correspondences identified in lightfield
        groundTruthWorldPoints = math.Point()   % Ground truth points in world coordinate system
        reconstructionMethod   = camera.models.plenoptic.enums.ReconstructionTypes.POINT_RECONSTRUCTION() 
                % Reconstruction method to be used
    end
    
    properties (Dependent)
        numberPoints                % Number of points considered for the reconstruction
        reconstructionMethodHandle  % Reconstruction method handle
        groundTruthCameraPoints     % Ground truth camera points
        estimatedCameraPoints       % Estimated camera points from stereo reconstruction
    end
    
    methods
        function self = PointReconstruction(varargin)
            %
            % Create reconstruction instance. In this object, each item in
            % correspondences corresponds to a different point.
            %
            % INPUTS:
            %   1. plenoptic   - plenoptic camera instance.
            %   2. correspondences   - correspondences identified in 
            %   ligthfield.
            %   3. groundTruthWorldPoints - ground truth points in world
            %   coordinate system.
            %   4. reconstructionMethod   - reconstruction method.
            %
            narginchk(0,4);
            
            if ~isempty(varargin)
                if nargin >= 4
                    self.reconstructionMethod = varargin{4};
                end
                
                if nargin >= 3
                    self.groundTruthWorldPoints = varargin{3};
                end
                
                if nargin >= 2
                    self.correspondences = varargin{2};
                end
                
                if nargin >= 1
                    self.plenoptic = varargin{1};
                end
            end
        end
        
        function self = set.correspondences(self,newImagePoints)
            self.correspondences = newImagePoints;
        end
        
        function self = set.groundTruthWorldPoints(self,newGroundTruthWorldPoints)
            self.groundTruthWorldPoints = utils.enums.Classes.POINT().convert(newGroundTruthWorldPoints);
        end

        function self = set.reconstructionMethod(self,newReconstructionMethod)
            self.reconstructionMethod = newReconstructionMethod;
        end
        
        function self = set.plenoptic(self,newPlenoptic)
            self.plenoptic = utils.enums.Classes.STANDARD_PLENOPTIC_CAMERA().convert(newPlenoptic);
        end
        
        function number = get.numberPoints(self)
            number = length(self.correspondences);
        end
        
        function handle = get.reconstructionMethodHandle(self)
            % This instance is used with the goal of representing image ray
            % coordinates obtained from the lightfield. Thus, the
            % coordinates have distortion included and should be removed
            % for the projection and reconstruction.
            FIT_EPI           = false;
            REMOVE_DISTORTION = true;
            % Define method handle
            if self.reconstructionMethod == camera.models.plenoptic.enums.ReconstructionTypes.POINT_RECONSTRUCTION()
                handle = @(imageRays_ijkl) reconstruct(self.plenoptic, imageRays_ijkl, REMOVE_DISTORTION);
            elseif self.reconstructionMethod == camera.models.plenoptic.enums.ReconstructionTypes.POINT_RECONSTRUCTION_FITTING_LINES()
                handle = @(imageRays_ijkl) reconstructByImposingLines(self.plenoptic, imageRays_ijkl, REMOVE_DISTORTION, FIT_EPI);
            elseif self.reconstructionMethod == camera.models.plenoptic.enums.ReconstructionTypes.LINE_PARAMETER_RECONSTRUCTION()
                handle = @(imageRays_ijkl) reconstructUsingLineParameters(self.plenoptic, imageRays_ijkl, REMOVE_DISTORTION, FIT_EPI);
            end
        end
        
        function points = get.estimatedCameraPoints(self)
            % Estimate the camera point
            points = [];
            for iPoint = 1:self.numberPoints
                % Obtain image ray coordinates for correspondences
                imageRays_ijkl = self.correspondences(iPoint);
                
                % Obtain the camera point estimate by stereo reconstruction
                pointEstimate = self.reconstructionMethodHandle( imageRays_ijkl );

                % Concatenate point to estimates collection
                points = cat(2,points,pointEstimate.data);
            end
            points = math.Point(points,false);
        end
        
        function points = get.groundTruthCameraPoints(self)
            % Obtain ground truth points in camera coordinate system
            if self.groundTruthWorldPoints.numberVectors > 0
                points = self.plenoptic.extrinsic.obtainCameraPoints(self.groundTruthWorldPoints);
            else
                points = math.Point();
            end
        end
    end
end
