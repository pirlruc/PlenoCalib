classdef Camera < abstract.TemplatePlenopticCamera
    %CAMERA
    %   Lytro camera utility to decode lightfields and calibrate camera.
    
    properties (Dependent)
        numberMicrolensesTypes      % Number of different microlenses types
    end
    
    methods
        function self = Camera(varargin)
            %
            % Lytro camera instance.
            %
            narginchk(0,1);
            
            % Create super class
            self = self@abstract.TemplatePlenopticCamera();

            if ~isempty(varargin)
                if nargin >= 1
                    self.whiteImagesFilepath = varargin{1};
                end
            end
        end
        
        function number = get.numberMicrolensesTypes(~)
            number = 1;
        end
    end
    
    methods
        function [intrinsic,extrinsic,distortion,projectionErrors] = ...
                        calibrate( self, calibrationFolderPath, pattern ...
                                 , linearSolutionParameters ...
                                 , additionalOptions )
            %
            % Calibrate plenoptic camera.
            % 
            % INPUTS:
            %   01. calibrationFolderPath - local folder path to the 
            %   calibration images for the camera. The images are not kept
            %   in memory to avoid a huge increase in memory usage.
            %   02. pattern - information regarding calibration pattern.
            %   03. linearSolutionParameters - parameters for linear
            %   initialization for calibration of lytro cameras.
            %   04. additionalOptions - additional options for calibration.
            %
            narginchk(3,5);
            
            if nargin <= 3
                linearSolutionParameters = struct( 'observationMatrixMethod', camera.models.pinhole.enums.ObservationMatrixMethods.CROSS_PRODUCT() );
            end
            
            if nargin <= 4
                additionalOptions = struct( 'cornersTemplate'  ,'%s__CheckerCorners.mat' ...
                                          , 'calibrationFile'  ,'CalInfo.json' ...
                                          , 'lensletBorderSize', 1 ...
                                          , 'optimizationTol'  , 5e-5 );
            end
            if ~isfield(additionalOptions,'cornersTemplate')
                additionalOptions.cornersTemplate = '%s__CheckerCorners.mat';
            end
            if ~isfield(additionalOptions,'calibrationFile')
                additionalOptions.calibrationFile = 'CalInfo.json';
            end
            if ~isfield(additionalOptions,'optimizationTol')
                additionalOptions.optimizationTol = 5e-5;
            end
            if ~isfield(additionalOptions,'lensletBorderSize')
                additionalOptions.lensletBorderSize = 1;
            end
            
            LIGHTFIELD_FILENAME_TEMPLATE = '%s__Decoded.mat';
            CORNERS_CALIBRATION_FILENAME = 'CheckerboardCorners.mat';
            DANSEREAU_DISTORTION         = camera.models.plenoptic.enums.DistortionModels.DANSEREAU_ORIGINAL;
            
            % Define path for white images database and calibration options
            calibrationOptionsTemp = struct( 'ExpectedCheckerSize'       , pattern.numberRectangles ...
                                           , 'ExpectedCheckerSpacing_m'  , pattern.size ...
                                           , 'CheckerCornersFnamePattern', additionalOptions.cornersTemplate ...
                                           , 'Phase'            , 'Init' ...
                                           , 'CalInfoFname'     , additionalOptions.calibrationFile ...
                                           , 'OptTolX'          , additionalOptions.optimizationTol ...
                                           , 'LFFnamePattern'   , LIGHTFIELD_FILENAME_TEMPLATE ...
                                           , 'CheckerInfoFname' , CORNERS_CALIBRATION_FILENAME ... 
                                           , 'SaveResult'       , true ...
                                           , 'ForceRedoInit'    , true ...
                                           , 'Refine'           , true ...
                                           , 'ShowDisplay'      , false ...
                                           , 'LensletBorderSize', additionalOptions.lensletBorderSize ...
                                           , 'DistortModel'     , DANSEREAU_DISTORTION.keyword );

            % Perform calibration procedure
            projectionErrors    = [];
            calibrationComplete = false;
            firstDistortion     = false;
            while calibrationComplete == false
                % Obtain iteration phase in calibration procedure
                iteration.phase = calibrationOptionsTemp.Phase;
                
                % Initialize calibration with linear solution
                if strcmp(calibrationOptionsTemp.Phase,'Init') > 0
                    % Create file with corners from each pose
                    [corners,worldPoints,lightfieldSize] = self.loadCornersFromFolderPath( calibrationFolderPath ...
                                                                                         , calibrationOptionsTemp );

                    % Obtain the number of poses
                    numberPoses = size(corners,3);

                    % Obtain the homography for each pose
                    homographies = [];
                    for iPose = 1:numberPoses
                        % Obtain pose corners
                        poseCorners    = corners(:,:,iPose);

                        % Obtain correspondences
                        imageRaysForPoints_ijkl = self.obtainCorrespondencesFromCorners( poseCorners ...
                                                                                       , prod(pattern.numberRectangles) );

                        % Associate each ray to the corresponding
                        % world point
                        imageRays_ijkl     = [];
                        worldPointsForRays = [];
                        for iPoint = 1:length(imageRaysForPoints_ijkl)
                            % Obtain rays and world point for corner
                            cornerPoint = worldPoints.obtainVectors(iPoint);
                            cornerRays  = imageRaysForPoints_ijkl(iPoint);

                            % Concatenate rays and points
                            imageRays_ijkl     = cat(2,imageRays_ijkl,cornerRays.data);
                            worldPointsForRays = cat(2,worldPointsForRays,repmat(cornerPoint.data,1,cornerRays.numberVectors));
                        end

                        % Obtain viewpoint camera array homography matrices
                        viewpointHomography = camera.models.plenoptic.standard.ViewpointHomography.HomographyFromPointCorrespondences( ...
                                                    worldPointsForRays, imageRays_ijkl ...
                                                  , linearSolutionParameters.observationMatrixMethod );

                        % Add homography
                        homographies = cat(2,homographies,viewpointHomography);
                    end

                    % Obtain viewpoint camera array projection matrices
                    viewpointArrays  = camera.models.plenoptic.standard.ViewpointProjectionMatrix.ViewpointProjectionMatrixFromHomographies(homographies);

                    % Correct camera coordinate system origin
                    origin_uv = viewpointArrays(1).intrinsicMatrix(1:2,3);
                    cameraCoordinateOriginRay = lightfield.ImageRay([ 0 ...
                                                                    ; 0 ...
                                                                    ; origin_uv ],false);
                    calibrationOptionsTemp.cameraCoordinateOriginRay = cameraCoordinateOriginRay.data;

                    % Obtain intrinsic matrix H for lytro camera
                    plenopticCameras = camera.models.plenoptic.standard.Camera.CameraFromViewpointProjectionMatrix( ...
                                                viewpointArrays, lightfieldSize, cameraCoordinateOriginRay );

                    % Obtain execution information
                    generatedByInfo = struct( 'mfilename' , mfilename ...
                                            , 'time'      , datestr(now,'ddmmmyyyy_HHMMSS') ...
                                            , 'VersionStr', 1.0 );

                    % Report results in Dansereau's format
                    self.reportInDansereauFormat( calibrationFolderPath ...
                                                , plenopticCameras ...
                                                , calibrationOptionsTemp ...
                                                , generatedByInfo )
                    calibrationOptionsTemp.Phase = 'NoDistort';
                
                % Nonlinear optimization to refine linear solution without
                % distortion (only if distortion is selected).
                elseif strcmp(calibrationOptionsTemp.Phase,'NoDistort') > 0
                    calibrationOptionsTemp = camera.lytro.utils.LFCalRefineModified( calibrationFolderPath ...
                                                                                   , calibrationOptionsTemp );
                    close(gcf);
                    calibrationOptionsTemp.Phase = 'WithDistort';
                    
                % Nonlinear optimization with distortion (only if
                % distortion is selected).
                elseif strcmp(calibrationOptionsTemp.Phase,'WithDistort') > 0
                    calibrationOptionsTemp     = camera.lytro.utils.LFCalRefineModified( calibrationFolderPath ...
                                                                                       , calibrationOptionsTemp );
                    close(gcf);

                    if firstDistortion == false
                        if isempty(DANSEREAU_DISTORTION.secondParameters) == true
                            calibrationOptionsTemp.Phase = 'Refine';
                        else
                            calibrationOptionsTemp.Phase = 'WithDistort';
                        end
                        firstDistortion = true;
                    else
                        calibrationOptionsTemp.Phase = 'Refine';
                    end

                % Nonlinear optimization to refine solution
                else
                    calibrationOptionsTemp = camera.lytro.utils.LFCalRefineModified( calibrationFolderPath ...
                                                                                   , calibrationOptionsTemp );
                    close(gcf);
                    calibrationComplete = true;
                end
                
                % If projection errors are required
                if nargout > 3
                    % Create projection error collection
                    iteration.errors = camera.lytro.ProjectionErrorCollection.ProjectionErrorCollectionFromCalibrationFile( ...
                                                [calibrationFolderPath filesep calibrationOptionsTemp.CalInfoFname] ...
                                              , camera.models.plenoptic.enums.ReconstructionTypes.POINT_RECONSTRUCTION );

                    % Add iteration results to projection errors
                    projectionErrors = cat(1,projectionErrors,iteration);
                    fprintf(' ###### PHASE = %s ###### \n',iteration.phase);
                    fprintf('Dansereau Point to Ray Distance: %f \n', mean(iteration.errors.dansereauPointToRayDistance));
                    fprintf('Projection Error               : %f \n', mean(iteration.errors.projectionError)            );
                    fprintf('Reconstruction Error           : %f \n', mean(iteration.errors.reconstructionError)        );
                end
            end
            
            % Obtain final camera model parameters
            [intrinsic,extrinsic,distortion] = ...
                    camera.lytro.Camera.obtainCalibrationParametersFromCalibrationFile( ...
                            [calibrationFolderPath filesep calibrationOptionsTemp.CalInfoFname] );
        end
    end
        
    methods (Static)
        function [intrinsic,extrinsic,distortion,lightfieldSize,reprojectionError] = ...
                        obtainCalibrationParametersFromCalibrationFile(calibrationFilepath)
            %
            % Obtain calibration parameters from calibration file.
            %
            % INPUTS:
            %   1. calibrationFilepath - filepath to the calibration file.
            %
            narginchk(1,1);
            
            try
                [ rawExtrinsic, intrinsic, rawDistortion, calibrationOptions, reprojectionError] = ...
                    LFStruct2Var( LFReadMetadata(calibrationFilepath) ...
                                , 'EstCamPosesV' ...
                                , 'EstCamIntrinsicsH' ...
                                , 'EstCamDistortionV' ...
                                , 'CalOptions' ...
                                , 'ReprojectionError' );
            catch
                [ rawExtrinsic, intrinsic, rawDistortion, calibrationOptions ] = ...
                    LFStruct2Var( LFReadMetadata(calibrationFilepath) ...
                                , 'EstCamPosesV' ...
                                , 'EstCamIntrinsicsH' ...
                                , 'EstCamDistortionV' ...
                                , 'CalOptions' );
                reprojectionError = nan;
            end                        
            % Obtain lightfield size from calibration options
            if ~isfield(calibrationOptions,'LFSize')
                rootPath = fileparts(calibrationFilepath);
                load([rootPath filesep calibrationOptions.CheckerInfoFname],'LFSize');
                calibrationOptions.LFSize = LFSize;
            end
            toolboxFormat      = lightfield.enums.Formats.DANSEREAU_FORMAT();
            internalFormatSize = calibrationOptions.LFSize(toolboxFormat.decodingOrder);
            lightfieldSize     = lightfield.LightfieldSize(internalFormatSize);
            
            % Divide extrinsic parameters into rotation vector and
            % translation vector
            extrinsic.rotationVector    = rawExtrinsic(:,4:6);
            extrinsic.translationVector = rawExtrinsic(:,1:3);

            % Segregate distortion parameters according to distortion model
            distortion = camera.models.plenoptic.Distortion();
            if ~isempty(rawDistortion)
                switch calibrationOptions.DistortModel
                    case camera.models.plenoptic.enums.DistortionModels.DANSEREAU_ORIGINAL().keyword
                        distortion = camera.models.plenoptic.enums.DistortionModels.DANSEREAU_ORIGINAL().decode(rawDistortion);
                    case camera.models.plenoptic.enums.DistortionModels.DANSEREAU_SIMPLIFIED().keyword
                        distortion = camera.models.plenoptic.enums.DistortionModels.DANSEREAU_SIMPLIFIED().decode(rawDistortion);
                    case camera.models.plenoptic.enums.DistortionModels.BOK_ORIGINAL().keyword
                        distortion = camera.models.plenoptic.enums.DistortionModels.BOK_ORIGINAL().decode(rawDistortion);
                    case camera.models.plenoptic.enums.DistortionModels.BOK_DISTORTION_CENTER().keyword
                        distortion = camera.models.plenoptic.enums.DistortionModels.BOK_DISTORTION_CENTER().decode(rawDistortion);
                    case camera.models.plenoptic.enums.DistortionModels.JOHANNSEN_ORIGINAL().keyword
                        distortion = camera.models.plenoptic.enums.DistortionModels.JOHANNSEN_ORIGINAL().decode(rawDistortion);
                    case camera.models.plenoptic.enums.DistortionModels.ZELLER_ORIGINAL().keyword
                        distortion = camera.models.plenoptic.enums.DistortionModels.ZELLER_ORIGINAL().decode(rawDistortion);
                    case camera.models.plenoptic.enums.DistortionModels.ZHANG_ORIGINAL().keyword
                        distortion = camera.models.plenoptic.enums.DistortionModels.ZHANG_ORIGINAL().decode(rawDistortion);
                    case camera.models.plenoptic.enums.DistortionModels.ZHANG_DISTORTION_CENTER().keyword
                        distortion = camera.models.plenoptic.enums.DistortionModels.ZHANG_DISTORTION_CENTER().decode(rawDistortion);
                    case camera.models.plenoptic.enums.DistortionModels.MONTEIRO_ORIGINAL().keyword
                        distortion = camera.models.plenoptic.enums.DistortionModels.MONTEIRO_ORIGINAL().decode(rawDistortion);
                end
            end
        end
        
        function [imageCorners, worldCornersGroundTruth, lightfieldSize] = ...
                        obtainCornersFromCornersFile(cornersFilepath,patternSize)
            % 
            % Obtain corner points from corners information file.
            %
            % INPUTS:
            %   1. cornersFilepath - filepath to the corners information 
            %   file.
            %   2. patternSize     - number of corners in checkerboard
            %   pattern. This is used to correct the origin of the pattern.
            %
            narginchk(1,2);
            
            if nargin <= 1
                patternSize = image.Pixel();
            end
            
            % Convert to pixel object
            patternSize = utils.enums.Classes.PIXEL().convert(patternSize);
            
            % Obtain information from file
            fileData = load(cornersFilepath);
            
            % During the calibration procedure, the pixels i and j were
            % switched. Thus, transpose the information of the observed
            % corners to have the correct order. Pose appears in the last
            % dimension
            %
            % Corners from one ligthfield image
            if isfield(fileData, 'CheckerCorners')
                % In this file, the components are in different dimensions.
                % The number of vectors is in the first component instead
                % of the second component, thus, we must switch the
                % dimensions. Additionally, the cell dimensions are defined
                % in (j,i) instead of (i,j).
                imageCorners = cellfun( @(x) permute(x,[2,1]) ...
                                      , fileData.CheckerCorners, 'UniformOutput', false );
                imageCorners = permute(imageCorners,[2,1]);
                worldCornersGroundTruth = math.Point();
                
                % If pattern size is given, correct origin of pattern only
                % for viewpoints with that number of corners
                if patternSize.numberVectors > 0
                    for viewpoint_i = 1:size(imageCorners,1)
                        for viewpoint_j = 1:size(imageCorners,2)
                            % Obtain viewpoint corners
                            viewpointCorners = image.Pixel(imageCorners{viewpoint_i,viewpoint_j},false);
                            
                            % Correct origin to top left corner
                            if viewpointCorners.numberVectors == prod(patternSize.data)
                                viewpointCorners = image.corners.CheckerboardCorners.ensureTopLeftOrigin( viewpointCorners ...
                                                                                                        , patternSize );

                                % Set corners in cell structure
                                imageCorners{viewpoint_i,viewpoint_j} = viewpointCorners.data;
                            end
                        end
                    end
                end
            
            % Corners from several lightfield images
            else
                % The variable CheckerObs is defined as (j,i,poses) instead
                % of (poses,i,j).
                imageCorners = permute(fileData.CheckerObs,[3,2,1]);
                worldCornersGroundTruth = math.Point(fileData.IdealChecker,false);
            end
            
            % Obtain lightfield size
            toolboxFormat      = lightfield.enums.Formats.DANSEREAU_FORMAT();
            internalFormatSize = fileData.LFSize(toolboxFormat.decodingOrder);
            lightfieldSize     = lightfield.LightfieldSize(internalFormatSize);
        end
        
        function correspondences = obtainCorrespondencesFromCorners( cornerPoints_kl ...
                                                                   , numberCorners )
            % 
            % Obtain the image rays corresponding to the corner image 
            % points obtained from the viewpoints images. 
            % 
            % INPUTS: 
            %   1. cornerPoints_kl - interest corner points identified in 
            %   the viewpoint images. The structure should be given as a 
            %   cell matrix with dimensions:
            %           numberPixels_i x numberPixels_j
            %   2. numberCorners - number of corners to be found. This 
            %   assumes that the corner points are ordered in each cell
            %   viewpoint entry. All viewpoint entries with a different
            %   number of corner points will be discarded
            %
            narginchk(2,2);
            
            % Obtain projection rays for all interest points
            rays_ijkl = [];
            for pixel_i = 1:size(cornerPoints_kl,1)
                for pixel_j = 1:size(cornerPoints_kl,2)
                    % Consider only viewpoint images with all interest 
                    % points to ensure the ordering of the interest points.
                    % The number of correspondences is given by the 
                    % different columns.
                    imagePixels = image.Pixel(cornerPoints_kl{pixel_i,pixel_j}, false);
                    if imagePixels.numberVectors ~= numberCorners ...
                    || numberCorners == 0
                        continue
                    end
                    
                    rays_ijkl = cat(2,rays_ijkl, [ repmat([pixel_i;pixel_j],1,numberCorners) ... % (i,j)
                                                 ; imagePixels.data ]);                          % (k,l)
                end
            end
            rays_ijkl = lightfield.ImageRay(rays_ijkl,false);
            
            % Segregate each projection ray to the corresponding point
            if numberCorners == 0 || rays_ijkl.numberVectors == 0
                correspondences = lightfield.ImageRay();
            else
                correspondences = [];
                for iPoint = 1:numberCorners
                    % Obtain rays for corner
                    cornerRays      = rays_ijkl.obtainVectors(iPoint:numberCorners:rays_ijkl.numberVectors);
                    correspondences = cat( 1, correspondences ...
                                         , lightfield.ImageRay( cornerRays.data, false ) );
                end
            end
        end
        
        function [imageCorners, worldCornersGroundTruth, lightfieldSize] = loadCornersFromFolderPath( imagesPath ...
                                                                                                    , calibrationOptions )
            %
            % Load corners from folder path.
            %
            % INPUTS:
            %   1. imagesPath  - corners folder path.
            %   2. patternSize - number of corners in checkerboard pattern.
            %   This is used to correct the origin of the pattern.
            %   3. patternLength - length between consecutive corners of
            %   checkerboard pattern.
            %
            narginchk(2,2);
            
            % Convert to pixel object
            patternSize   = utils.enums.Classes.PIXEL().convert(calibrationOptions.ExpectedCheckerSize);
            patternLength = utils.enums.Classes.PIXEL().convert(calibrationOptions.ExpectedCheckerSpacing_m);

            % Corners filename template
            DANSEREAU_CORNERS_STRING   = '__CheckerCorners';
            DANSEREAU_CORNERS_FILENAME = 'CheckerboardCorners.mat';

            % Obtain corner files
            cornerFiles    = dir([ imagesPath filesep '**' filesep '*' DANSEREAU_CORNERS_STRING '.mat' ]);
            cornerFilepath = [imagesPath filesep DANSEREAU_CORNERS_FILENAME];

            % Load first file just to get the number of viewpoints
            filepath             = [cornerFiles(1).folder filesep cornerFiles(1).name];
            load(filepath,'LFSize','CamInfo','LensletGridModel','DecodeOptions');
            [~,~,lightfieldSize] = camera.lytro.Camera.obtainCornersFromCornersFile(filepath);
            numberPoses          = length(cornerFiles);
            numberViewpoints_i   = lightfieldSize.numberPixels_i;
            numberViewpoints_j   = lightfieldSize.numberPixels_j;
            
            % Obtain and validate corners of each file and pose
            imageCorners = cell(numberViewpoints_i,numberViewpoints_j,numberPoses);
            for iFile = 1:numberPoses
                % Obtain filename
                filepath = [cornerFiles(iFile).folder filesep cornerFiles(iFile).name];
                
                % Obtain corners from file
                cornerPoses = camera.lytro.Camera.obtainCornersFromCornersFile( filepath ...
                                                                              , patternSize );
                
                % Remove viewpoints that do not have the full grid and that
                % do not have the top left corner as the origin
                for viewpoint_i = 1:numberViewpoints_i
                    for viewpoint_j = 1:numberViewpoints_j
                        % Obtain viewpoint corners
                        viewpointCorners = image.Pixel( cornerPoses{viewpoint_i,viewpoint_j} ...
                                                      , false );

                        % Obtain centroid and pattern origin
                        centroid = mean(viewpointCorners.data,2);
                        origin   = viewpointCorners.obtainVectors(1);

                        % Grid must be fully observable
                        if viewpointCorners.numberVectors ~= prod(patternSize.data)
                            viewpointCorners.data = [];

                        % Pattern must have the top left corner as the
                        % origin
                        elseif ~(all(origin.data < centroid))
                            viewpointCorners.data = [];
                        end

                        % Set corners in cell structure
                        imageCorners{viewpoint_i,viewpoint_j,iFile} = viewpointCorners.data;
                    end
                end
            end

            % Transform corners of each pose from (i,j,poses) to 
            % (poses,j,i)
            CheckerObs = permute(imageCorners,[3,2,1]);

            % Obtain world pattern points
            worldCornersGroundTruth = image.corners.CheckerboardCorners.obtainGroundTruthCheckerboardCorners( patternSize ...
                                                                                                            , patternLength );
            IdealChecker            = worldCornersGroundTruth.data;

            % Save corners from all poses used in calibration
            GeneratedByInfo = struct( 'mfilename' , mfilename ...
                                    , 'time'      , datestr(now,'ddmmmyyyy_HHMMSS') ...
                                    , 'VersionStr', 1.0 );
            CalOptions      = calibrationOptions;
            save( cornerFilepath ...
                , 'LFSize', 'CamInfo', 'LensletGridModel', 'DecodeOptions' ...
                , 'GeneratedByInfo', 'CalOptions', 'CheckerObs', 'IdealChecker' );
        end
        
        function reportInDansereauFormat( calibrationFolderPath ...
                                        , plenopticCameras ...
                                        , calibrationOptions ...
                                        , generatedByInfo )
            %
            % Report results in Dansereau's format.
            %
            % INPUTS:
            %   1. calibrationFolderPath - folder path to save results.
            %   2. plenopticCameras      - calibration result.
            %   3. calibrationOptions    - input for calibration.
            %   4. generatedByInfo       - script information.
            %
            narginchk(3,4);

            if nargin <= 3
                % Obtain execution information
                generatedByInfo = struct( 'mfilename' , mfilename ...
                                        , 'time'      , datestr(now,'ddmmmyyyy_HHMMSS') ...
                                        , 'VersionStr', 1.0 );
            end                    
            
            % Corners filename template
            DANSEREAU_CORNERS_FILENAME     = 'CheckerboardCorners.mat';
            DANSEREAU_CALIBRATION_FILENAME = 'CalInfo.json';

            % Obtain corner and calibration files
            cornerFilepath      = [calibrationFolderPath filesep DANSEREAU_CORNERS_FILENAME];
            calibrationFilepath = [calibrationFolderPath filesep DANSEREAU_CALIBRATION_FILENAME];

            % Load information from lightfield
            load(cornerFilepath,'LFSize','CamInfo','LensletGridModel','DecodeOptions');

            % Save calibration initialization 
            GeneratedByInfo   = generatedByInfo;
            CalOptions        = calibrationOptions;
            CalOptions.LFSize = LFSize;
            EstCamIntrinsicsH = plenopticCameras(1).intrinsicMatrix;
            EstCamDistortionV = [];
            ReprojectionError = [];
            
            % Obtain extrinsic parameters for each pose
            numberPoses  = length(plenopticCameras);
            EstCamPosesV = zeros(numberPoses,6);
            for iPose = 1:numberPoses
                % Obtain extrinsic parameters
                extrinsic = plenopticCameras(iPose).extrinsic;
                
                % Obtain rotation matrix and translation vector
                EstCamPosesV(iPose,1:3) = extrinsic.translationVector.data';
                EstCamPosesV(iPose,4:6) = extrinsic.rotationMatrix.vector';
            end
            LFWriteMetadata( calibrationFilepath, LFVar2Struct( GeneratedByInfo ...
                                                              , LensletGridModel ...
                                                              , EstCamIntrinsicsH ...
                                                              , EstCamDistortionV ...
                                                              , EstCamPosesV ...
                                                              , CamInfo ...
                                                              , CalOptions ...
                                                              , DecodeOptions ...
                                                              , ReprojectionError ));
            fprintf(' ---Finished calibration initialization---\n');
            fprintf('Estimate of camera intrinsics: \n');
            disp(EstCamIntrinsicsH);
        end
        
        function [corners_ijkl,worldPoints_xyz,centers_uv,calibrationInformation] = loadCornersInformationForPoses( calibrationFolderPath ...
                                                                                                                  , calibrationOptions ...
                                                                                                                  , calibrationInformation )
            %
            % Obtain corners information.
            %
            % INPUTS:
            %   1. calibrationFolderPath  - folder path containing corners
            %   information.
            %   2. calibrationOptions     - options to consider during 
            %   calibration.
            %   3. calibrationInformation - geenral information regarding 
            %   the calibration procedure.
            %
            narginchk(3,3);

            % Create file with corners from each pose
            calibrationOptions.ExpectedCheckerSize      = calibrationInformation.patternInformation.numberRectangles.data;
            calibrationOptions.ExpectedCheckerSpacing_m = calibrationInformation.patternInformation.sizeRectangles.data;
            [corners,worldPoints,lightfieldSize] = camera.lytro.Camera.loadCornersFromFolderPath( ...
                                                                calibrationFolderPath ...
                                                              , calibrationOptions );
                        
            % Obtain the number of poses
            calibrationInformation.numberPoses = size(corners,3);
                        
            % For each pose, concatenate information of corners
            numberRays      = 0;
            corners_ijkl    = cell(calibrationInformation.numberMicrolensesTypes,calibrationInformation.numberPoses);
            worldPoints_xyz = cell(calibrationInformation.numberMicrolensesTypes,calibrationInformation.numberPoses);
            centers_uv      = cell(calibrationInformation.numberMicrolensesTypes,calibrationInformation.numberPoses);
            for iLens = 1:calibrationInformation.numberMicrolensesTypes
                for iPose = 1:calibrationInformation.numberPoses
                    % Obtain pose corners
                    poseCorners = corners(:,:,iPose);

                    % Obtain correspondences
                    imageRaysForPoints_ijkl = camera.lytro.Camera.obtainCorrespondencesFromCorners( ...
                                                    poseCorners ...
                                                  , prod(calibrationInformation.patternInformation.numberRectangles.data) );

                    % Associate each ray to the corresponding
                    % world point
                    imageRays_ijkl     = [];
                    worldPointsForRays = [];
                    for iPoint = 1:length(imageRaysForPoints_ijkl)
                        % Obtain rays and world point for corner
                        cornerPoint = worldPoints.obtainVectors(iPoint);
                        cornerRays  = imageRaysForPoints_ijkl(iPoint);

                        % Concatenate rays and points
                        imageRays_ijkl     = cat(2,imageRays_ijkl,cornerRays.data);
                        worldPointsForRays = cat(2,worldPointsForRays,repmat(cornerPoint.data,1,cornerRays.numberVectors));
                    end
                    
                    % Define corners and world points for pose
                    corners_ijkl{iLens,iPose}    = lightfield.ImageRay(imageRays_ijkl,false);
                    worldPoints_xyz{iLens,iPose} = math.Point(worldPointsForRays,false);
                    
                    % Update number of rays detected
                    numberRays = numberRays + corners_ijkl{iLens,iPose}.numberVectors;
                end
            end
            calibrationInformation.numberRays     = numberRays;
            calibrationInformation.lightfieldSize = lightfieldSize;
        end
    end
end
