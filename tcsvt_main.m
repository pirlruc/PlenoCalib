% Setup environment
ins = setup.Setup();

% If there is an error while installing Bok toolbox, download directly and
% put it in toolboxes folder with the zipped file and rerun script
fprintf('Installing dependencies...\n');
ins.install();      

% Load sample datasets
[datasets,patterns,whiteImages] = loadDansereauCalibrationDatasets();

% Pre-process datasets
preProcessDatasets(datasets,whiteImages);

% Calibrate datasets
calibrateDatasets(datasets,whiteImages,patterns);



% Auxiliar Functions
function [datasets,patterns,whiteImages] = loadDansereauCalibrationDatasets()
    %     
    % Load Dansereau calibration datasets and extract white images
    % database.
    % 
    
    % Load Dansereau Calibration Datasets
    whiteImages = setup.enums.datasets.LightfieldDatasets.LFTOOLBOX_SAMPLES;
    fprintf('Downloading white images...\n');
    whiteImages.download();
    datasetA    = setup.enums.datasets.LightfieldDatasets.CALIBRATION_DATASET_CVPR2013_A;
    datasetB    = setup.enums.datasets.LightfieldDatasets.CALIBRATION_DATASET_CVPR2013_B;
    datasetC    = setup.enums.datasets.LightfieldDatasets.CALIBRATION_DATASET_CVPR2013_C;
    datasetD    = setup.enums.datasets.LightfieldDatasets.CALIBRATION_DATASET_CVPR2013_D;
    datasetE    = setup.enums.datasets.LightfieldDatasets.CALIBRATION_DATASET_CVPR2013_E;
    datasets    = [datasetA;datasetB;datasetC;datasetD;datasetE];
    fprintf('Downloading Dataset A from Dansereau...\n');
    datasetA.download();    
    fprintf('Downloading Dataset B from Dansereau...\n');
    datasetB.download();
    fprintf('Downloading Dataset C from Dansereau...\n');
    datasetC.download();
    fprintf('Downloading Dataset D from Dansereau...\n');
    datasetD.download();
    fprintf('Downloading Dataset E from Dansereau...\n');
    datasetE.download();

    % Obtain lytro archive folder path and white images path
    lytroArchivePath = whiteImages.dataset.folderPath;
    whiteImagesPath  = [lytroArchivePath filesep 'WhiteImageDatabase.mat'];

    % Obtain lytro camera and white images database
    if ~exist(whiteImagesPath,'file')
        LFUtilUnpackLytroArchive(lytroArchivePath);
        LFUtilProcessWhiteImages(lytroArchivePath);
    end

    % Obtain pattern information
    patternA = struct('numberRectangles',[19,19],'size',[3.61,3.61] .* 1e-3);
    patternB = struct('numberRectangles',[19,19],'size',[3.61,3.61] .* 1e-3);
    patternC = struct('numberRectangles',[19,19],'size',[7.22,7.22] .* 1e-3);
    patternD = struct('numberRectangles',[19,19],'size',[7.22,7.22] .* 1e-3);
    patternE = struct('numberRectangles',[ 8, 6],'size',[35.1,35.0] .* 1e-3);
    patterns = [patternA;patternB;patternC;patternD;patternE];
end

function preProcessDatasets(datasets,whiteImages)
    %
    % Pre-process datasets for calibration.
    %
    % INPUTS:
    %   1. datasets     - calibration datasets information.
    %   2. whiteImages  - white images information.
    %
    narginchk(2,2);

    % Obtain lytro archive folder path and white images path
    lytroArchivePath = whiteImages.dataset.folderPath;
    whiteImagesPath  = [lytroArchivePath filesep 'WhiteImageDatabase.mat'];

    % Obtain lytro camera and white images database
    if ~exist(whiteImagesPath,'file')
        LFUtilUnpackLytroArchive(lytroArchivePath);
        LFUtilProcessWhiteImages(lytroArchivePath);
    end

    % Define decoding properties and corner extraction options
    decodeOptions = struct( 'WhiteImageDatabasePath', whiteImagesPath ...
                          , 'ResampMethod'          , 'triangulation' ...
                          , 'DoDehex'               , true ...
                          , 'DoSquareST'            , true ...
                          , 'OptionalTasks'         , 'ColourCorrect' );
    cornerOptions = struct('ShowDisplay', false);

    % Process each dataset
    for iDataset = 1:length(datasets)
        % Decode lightfields from lytro image
        LFUtilDecodeLytroFolder(datasets(iDataset).dataset.folderPath, [], decodeOptions);

        % Extract corners from lytro images
        LFCalFindCheckerCorners(datasets(iDataset).dataset.folderPath, cornerOptions);
    end
end

function calibrateDatasets(datasets,whiteImages,patterns)
    %
    % Calibrate datasets.
    %
    % INPUTS:
    %   1. datasets     - calibration datasets information.
    %   2. whiteImages  - white images information.
    %   3. patterns     - calibration patterns information.
    %
    narginchk(3,3);
    
    % Obtain white image path
    lytroArchivePath = whiteImages.dataset.folderPath;
    whiteImagesPath  = [lytroArchivePath filesep 'WhiteImageDatabase.mat'];
    lytro            = camera.lytro.Camera();
    lytro.whiteImagesFilepath = whiteImagesPath;

    MONTEIRO_PARAMETERS = struct('observationMatrixMethod',camera.models.pinhole.enums.ObservationMatrixMethods.CROSS_PRODUCT());
    for iDataset = 1:length(datasets)
        [intrinsic,extrinsic,distortion,projectionErrors] = ...
                           lytro.calibrate( datasets(iDataset).dataset.folderPath ...
                                          , patterns(iDataset) ...
                                          , MONTEIRO_PARAMETERS );
        fprintf(' ###### FINAL CALIBRATION DATASET %s ######\n', char(datasets(iDataset)));
        fprintf('LFIM:\n');
        disp(intrinsic);
        fprintf('Errors:\n');
        fprintf('Dansereau Point to Ray Distance: %f \n', mean(projectionErrors(end).errors.dansereauPointToRayDistance));
        fprintf('Projection Error               : %f \n', mean(projectionErrors(end).errors.projectionError)            );
        fprintf('Reconstruction Error           : %f \n', mean(projectionErrors(end).errors.reconstructionError)        );
    end
end
