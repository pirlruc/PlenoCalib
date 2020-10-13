DISPARITIES = [0.1,0.5];

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

% Shear datasets
shearDatasets(datasets,DISPARITIES);

% Calibrate datasets
calibrateDatasets(datasets,whiteImages,patterns,DISPARITIES);

% Show results of calibration
for iDataset = 1:length(datasets)
    plotModelVersusEstimated([pwd filesep 'calibration.mat'],DISPARITIES);
end







% Auxiliar Functions
function [datasets,patterns,whiteImages] = loadDansereauCalibrationDatasets()
    %     
    % Load Dansereau calibration dataset and extract white images
    % database.
    % 
    
    % Load Dansereau Calibration Datasets
    whiteImages = setup.enums.datasets.LightfieldDatasets.LFTOOLBOX_SAMPLES;
    fprintf('Downloading white images...\n');
    whiteImages.download();
    datasetA    = setup.enums.datasets.LightfieldDatasets.CALIBRATION_DATASET_CVPR2013_A;
    datasets    = datasetA;
    fprintf('Downloading Dataset A from Dansereau...\n');
    datasetA.download();    

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
    patterns = patternA;
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

function shearDatasets(datasets,disparities)
    %
    % Shear datasets for calibration.
    %
    % INPUTS:
    %   1. datasets    - calibration datasets information.
    %   2. disparities - disparity values to consider for shearing.
    %
    narginchk(1,2);

    if nargin <= 1
        disparities = 0.1:0.1:2;
    end

    % Process each dataset
    for iDataset = 1:length(datasets)
        calibrationFolderPath = fileparts(datasets(iDataset).dataset.folderPath);

        % Obtain files in dataset
        [datasetFiles, rootFolder] = LFFindFilesRecursive( datasets(iDataset).dataset.folderPath ...
                                                         , 'IMG_*__Decoded.mat' );

        % Create disparity folders
        for iDisparity = 1:length(disparities)
            disparity       = disparities(iDisparity);
            disparityFolder = [calibrationFolderPath filesep sprintf('disparity_%02d',round(disparity * 10))];
            mkdir(disparityFolder);
        end
        
        % Create sheared lightfields
        for iFile = 1:length(datasetFiles)
            % Obtain lightfield filename
            filepath = [rootFolder filesep datasetFiles{iFile}];
            load(filepath);
            [~,filename,extension] = fileparts(datasetFiles{iFile});
            filename = [filename,extension];

            % Obtain lightfield
            lfield = lightfield.lytro.Lightfield.LightfieldFromDansereauFile(filepath);

            % Shear and save lightfield
            for iDisparity = 1:length(disparities)
                disparity       = disparities(iDisparity);
                disparityFolder = [calibrationFolderPath filesep sprintf('disparity_%02d',round(disparity * 10))];

                shearedLfield = lfield.shear(disparity,lightfield.enums.ShearMethods.CUBIC,[5,5,nan,nan],true);
                data       = shearedLfield.TOOLBOX_FORMAT.encode(shearedLfield.data);
                confidence = shearedLfield.TOOLBOX_FORMAT.encode(shearedLfield.confidence);
                LF         = cat(5,data,confidence);

                save([ disparityFolder filesep filename ] ...
                     , 'DecodeOptions','GeneratedByInfo','LF','LFMetadata','LensletGridModel','RectOptions','WhiteImageMetadata' );
            end
        end

        % Extract corners in disparity folder path
        cornerOptions = struct('ShowDisplay', false);
        for iDisparity = 1:length(disparities)
            disparity       = disparities(iDisparity);
            disparityFolder = [calibrationFolderPath filesep sprintf('disparity_%02d',round(disparity * 10))];
            % Extract corners from lytro images
            LFCalFindCheckerCorners(disparityFolder, cornerOptions);
        end        
    end
end
        
function calibrateDatasets(datasets,whiteImages,patterns,disparities)
    %
    % Calibrate datasets.
    %
    % INPUTS:
    %   1. datasets     - calibration datasets information.
    %   2. whiteImages  - white images information.
    %   3. patterns     - calibration patterns information.
    %   4. disparities - disparity values to consider for shearing.
    %
    narginchk(3,4);
    
    if nargin <= 3
        disparities = 0.1:0.1:2;
    end
    
    % Obtain white image path
    lytroArchivePath = whiteImages.dataset.folderPath;
    whiteImagesPath  = [lytroArchivePath filesep 'WhiteImageDatabase.mat'];
    lytro            = camera.lytro.Camera();
    lytro.whiteImagesFilepath = whiteImagesPath;

    MONTEIRO_PARAMETERS = struct('observationMatrixMethod',camera.models.pinhole.enums.ObservationMatrixMethods.CROSS_PRODUCT());
    for iDataset = 1:length(datasets)
        calibrationFolderPath    = fileparts(datasets(iDataset).dataset.folderPath);
        [~,~,~,originalErrors] = ...
                           lytro.calibrate( datasets(iDataset).dataset.folderPath ...
                                          , patterns(iDataset) ...
                                          , MONTEIRO_PARAMETERS );
        fprintf('Errors:\n');
        fprintf('Dansereau Point to Ray Distance: %f \n', mean(originalErrors(end).errors.dansereauPointToRayDistance));
        fprintf('Projection Error               : %f \n', mean(originalErrors(end).errors.projectionError)            );
        fprintf('Reconstruction Error           : %f \n', mean(originalErrors(end).errors.reconstructionError)        );
        
        % Process calibration for each disparity dataset
        disparityErrors = cell(length(disparities),1);
        for iDisparity = 1:length(disparities)
            disparity       = disparities(iDisparity);
            disparityFolder = [calibrationFolderPath filesep sprintf('disparity_%02d',round(disparity * 10))];
            [~,~,~,disparityErrors{iDisparity}] = ...
                lytro.calibrate( disparityFolder ...
                               , patterns(iDataset) ...
                               , MONTEIRO_PARAMETERS );
            fprintf('Errors: %02d\n',round(disparity * 10));
            fprintf('Dansereau Point to Ray Distance: %f \n', mean(disparityErrors{iDisparity}(end).errors.dansereauPointToRayDistance));
            fprintf('Projection Error               : %f \n', mean(disparityErrors{iDisparity}(end).errors.projectionError)            );
            fprintf('Reconstruction Error           : %f \n', mean(disparityErrors{iDisparity}(end).errors.reconstructionError)        );
        end
        save('calibration.mat','originalErrors','disparityErrors');
    end
end

function plotModelVersusEstimated(resultsFilepath,disparities)
    %
    % Plot model for shearing vs. estimated parameters.
    %
    % INPUT:
    %   1. resultsFilepath - filepath to results file.
    %   2. disparities     - disparity values to consider for shearing.
    %
    narginchk(1,2);

    if nargin <= 1
        disparities = 0.1:0.1:2;
    end

    % Load results data
    resultsData = load(resultsFilepath);

    % Obtain reference viewpoint
    referenceViewpoint = round([ resultsData.originalErrors(end).errors.plenoptic.lightfieldSize.numberPixels_i ...
                               , resultsData.originalErrors(end).errors.plenoptic.lightfieldSize.numberPixels_j ] / 2 );

    % Obtain intrinsic matrix for original calibration dataset 
    originalPlenoptic = resultsData.originalErrors(end).errors.plenoptic;

    % Obtain plenoptic camera for each disparity value
    numberDisparities   = length(resultsData.disparityErrors);
    modelIntrinsic      = nan(numberDisparities,12);
    estimatedIntrinsic  = nan(numberDisparities,12);
    modelProjection     = nan(numberDisparities,30);
    estimatedProjection = nan(numberDisparities,30);
    for iDisparity = 1:numberDisparities
        % Obtain disparity value and estimate corresponding entry
        disparity = disparities(iDisparity);

        % Intrinsic Matrix H
        % Obtain shearing model intrinsic matrix
        samplingMatrix = eye(5,5);
        samplingMatrix(3,1)   =  disparity;
        samplingMatrix(4,2)   =  disparity;
        samplingMatrix(3,end) = -disparity * referenceViewpoint(1);
        samplingMatrix(4,end) = -disparity * referenceViewpoint(2);
        intrinsicModel    = originalPlenoptic.intrinsicMatrix * samplingMatrix;
        intrinsicModel    = camera.models.plenoptic.standard.Camera(intrinsicModel);
        originalViewpoint = camera.models.plenoptic.standard.ViewpointProjectionMatrix.ViewpointProjectionMatrixFromStandardPlenopticCamera(intrinsicModel);
        originalMicrolens = camera.models.plenoptic.standard.MicrolensProjectionMatrix.MicrolensProjectionMatrixFromStandardPlenopticCamera(intrinsicModel);

        % Get values for intrinsic matrix model
        modelIntrinsic(iDisparity,01) = intrinsicModel.h_si;
        modelIntrinsic(iDisparity,02) = intrinsicModel.h_sk;
        modelIntrinsic(iDisparity,03) = intrinsicModel.h_s;
        modelIntrinsic(iDisparity,04) = intrinsicModel.h_tj;
        modelIntrinsic(iDisparity,05) = intrinsicModel.h_tl;
        modelIntrinsic(iDisparity,06) = intrinsicModel.h_t;
        modelIntrinsic(iDisparity,07) = intrinsicModel.h_ui;
        modelIntrinsic(iDisparity,08) = intrinsicModel.h_uk;
        modelIntrinsic(iDisparity,09) = intrinsicModel.h_u;
        modelIntrinsic(iDisparity,10) = intrinsicModel.h_vj;
        modelIntrinsic(iDisparity,11) = intrinsicModel.h_vl;
        modelIntrinsic(iDisparity,12) = intrinsicModel.h_v;

        % Obtain intrinsic matrix for sheared calibration dataset 
        % and corresponding viewpoint and microlens arrays
        shearedPlenoptic = resultsData.disparityErrors{iDisparity}(end).errors.plenoptic;
        shearedPlenoptic.extrinsic = camera.models.ExtrinsicMatrix();
        shearedViewpoint = camera.models.plenoptic.standard.ViewpointProjectionMatrix.ViewpointProjectionMatrixFromStandardPlenopticCamera(shearedPlenoptic);
        shearedMicrolens = camera.models.plenoptic.standard.MicrolensProjectionMatrix.MicrolensProjectionMatrixFromStandardPlenopticCamera(shearedPlenoptic);

        % Get values for intrinsic matrix entries
        estimatedIntrinsic(iDisparity,01) = shearedPlenoptic.h_si;
        estimatedIntrinsic(iDisparity,02) = shearedPlenoptic.h_sk;
        estimatedIntrinsic(iDisparity,03) = shearedPlenoptic.h_s;
        estimatedIntrinsic(iDisparity,04) = shearedPlenoptic.h_tj;
        estimatedIntrinsic(iDisparity,05) = shearedPlenoptic.h_tl;
        estimatedIntrinsic(iDisparity,06) = shearedPlenoptic.h_t;
        estimatedIntrinsic(iDisparity,07) = shearedPlenoptic.h_ui;
        estimatedIntrinsic(iDisparity,08) = shearedPlenoptic.h_uk;
        estimatedIntrinsic(iDisparity,09) = shearedPlenoptic.h_u;
        estimatedIntrinsic(iDisparity,10) = shearedPlenoptic.h_vj;
        estimatedIntrinsic(iDisparity,11) = shearedPlenoptic.h_vl;
        estimatedIntrinsic(iDisparity,12) = shearedPlenoptic.h_v;

        % Viewpoints projection matrix
        % Obtain viewpoint projection model
        modelProjection(iDisparity,01) = originalViewpoint.intrinsicMatrix(1,1);
        modelProjection(iDisparity,02) = originalViewpoint.intrinsicMatrix(2,2);
        modelProjection(iDisparity,03) = originalViewpoint.intrinsicMatrix(1,3);
        modelProjection(iDisparity,04) = originalViewpoint.incrementIntrinsicMatrix_i(1,3);
        modelProjection(iDisparity,05) = originalViewpoint.intrinsicMatrix(2,3);
        modelProjection(iDisparity,06) = originalViewpoint.incrementIntrinsicMatrix_j(2,3);
        modelProjection(iDisparity,07) = originalViewpoint.extrinsicMatrix(1,4);
        modelProjection(iDisparity,08) = originalViewpoint.incrementExtrinsicMatrix_i(1,4);
        modelProjection(iDisparity,09) = originalViewpoint.extrinsicMatrix(2,4);
        modelProjection(iDisparity,10) = originalViewpoint.incrementExtrinsicMatrix_j(2,4);
        modelProjection(iDisparity,11) = originalViewpoint.extrinsicMatrix(3,4);
        modelProjection(iDisparity,12) = modelProjection(iDisparity,03) + modelProjection(iDisparity,04);
        modelProjection(iDisparity,13) = modelProjection(iDisparity,05) + modelProjection(iDisparity,06);
        modelProjection(iDisparity,14) = modelProjection(iDisparity,07) + modelProjection(iDisparity,08);
        modelProjection(iDisparity,15) = modelProjection(iDisparity,09) + modelProjection(iDisparity,10);

        % Obtain viewpoint projection estimate
        estimatedProjection(iDisparity,01) = shearedViewpoint.intrinsicMatrix(1,1);
        estimatedProjection(iDisparity,02) = shearedViewpoint.intrinsicMatrix(2,2);
        estimatedProjection(iDisparity,03) = shearedViewpoint.intrinsicMatrix(1,3);
        estimatedProjection(iDisparity,04) = shearedViewpoint.incrementIntrinsicMatrix_i(1,3);
        estimatedProjection(iDisparity,05) = shearedViewpoint.intrinsicMatrix(2,3);
        estimatedProjection(iDisparity,06) = shearedViewpoint.incrementIntrinsicMatrix_j(2,3);
        estimatedProjection(iDisparity,07) = shearedViewpoint.extrinsicMatrix(1,4);
        estimatedProjection(iDisparity,08) = shearedViewpoint.incrementExtrinsicMatrix_i(1,4);
        estimatedProjection(iDisparity,09) = shearedViewpoint.extrinsicMatrix(2,4);
        estimatedProjection(iDisparity,10) = shearedViewpoint.incrementExtrinsicMatrix_j(2,4);
        estimatedProjection(iDisparity,11) = shearedViewpoint.extrinsicMatrix(3,4);
        estimatedProjection(iDisparity,12) = estimatedProjection(iDisparity,03) + estimatedProjection(iDisparity,04);
        estimatedProjection(iDisparity,13) = estimatedProjection(iDisparity,05) + estimatedProjection(iDisparity,06);
        estimatedProjection(iDisparity,14) = estimatedProjection(iDisparity,07) + estimatedProjection(iDisparity,08);
        estimatedProjection(iDisparity,15) = estimatedProjection(iDisparity,09) + estimatedProjection(iDisparity,10);

        % Microlens projection matrix
        % Obtain microlens projection model
        modelProjection(iDisparity,16) = originalMicrolens.intrinsicMatrix(1,1);
        modelProjection(iDisparity,17) = originalMicrolens.intrinsicMatrix(2,2);
        modelProjection(iDisparity,18) = originalMicrolens.intrinsicMatrix(1,3);
        modelProjection(iDisparity,19) = originalMicrolens.incrementIntrinsicMatrix_k(1,3);
        modelProjection(iDisparity,20) = originalMicrolens.intrinsicMatrix(2,3);
        modelProjection(iDisparity,21) = originalMicrolens.incrementIntrinsicMatrix_l(2,3);
        modelProjection(iDisparity,22) = originalMicrolens.extrinsicMatrix(1,4);
        modelProjection(iDisparity,23) = originalMicrolens.incrementExtrinsicMatrix_k(1,4);
        modelProjection(iDisparity,24) = originalMicrolens.extrinsicMatrix(2,4);
        modelProjection(iDisparity,25) = originalMicrolens.incrementExtrinsicMatrix_l(2,4);
        modelProjection(iDisparity,26) = originalMicrolens.extrinsicMatrix(3,4);
        modelProjection(iDisparity,27) = modelProjection(iDisparity,18) + modelProjection(iDisparity,19);
        modelProjection(iDisparity,28) = modelProjection(iDisparity,20) + modelProjection(iDisparity,21);
        modelProjection(iDisparity,29) = modelProjection(iDisparity,22) + modelProjection(iDisparity,23);
        modelProjection(iDisparity,30) = modelProjection(iDisparity,24) + modelProjection(iDisparity,25);

        % Obtain microlens projection estimate
        estimatedProjection(iDisparity,16) = shearedMicrolens.intrinsicMatrix(1,1);
        estimatedProjection(iDisparity,17) = shearedMicrolens.intrinsicMatrix(2,2);
        estimatedProjection(iDisparity,18) = shearedMicrolens.intrinsicMatrix(1,3);
        estimatedProjection(iDisparity,19) = shearedMicrolens.incrementIntrinsicMatrix_k(1,3);
        estimatedProjection(iDisparity,20) = shearedMicrolens.intrinsicMatrix(2,3);
        estimatedProjection(iDisparity,21) = shearedMicrolens.incrementIntrinsicMatrix_l(2,3);
        estimatedProjection(iDisparity,22) = shearedMicrolens.extrinsicMatrix(1,4);
        estimatedProjection(iDisparity,23) = shearedMicrolens.incrementExtrinsicMatrix_k(1,4);
        estimatedProjection(iDisparity,24) = shearedMicrolens.extrinsicMatrix(2,4);
        estimatedProjection(iDisparity,25) = shearedMicrolens.incrementExtrinsicMatrix_l(2,4);
        estimatedProjection(iDisparity,26) = shearedMicrolens.extrinsicMatrix(3,4);
        estimatedProjection(iDisparity,27) = estimatedProjection(iDisparity,18) + estimatedProjection(iDisparity,19);
        estimatedProjection(iDisparity,28) = estimatedProjection(iDisparity,20) + estimatedProjection(iDisparity,21);
        estimatedProjection(iDisparity,29) = estimatedProjection(iDisparity,22) + estimatedProjection(iDisparity,23);
        estimatedProjection(iDisparity,30) = estimatedProjection(iDisparity,24) + estimatedProjection(iDisparity,25);
    end
    
    % Plot model against estimate for intrinsic matrix
    intrinsicLabels  = { 'h_{si}','h_{sk}','h_s','h_{tj}','h_{tl}','h_t' ...
                       , 'h_{ui}','h_{uk}','h_u','h_{vj}','h_{vl}','h_v' };
    location         = { '','','northeast','','','northeast' ...
                       , '','','northeast','','','' };
    for iEntry = 1:size(modelIntrinsic,2)
        % Obtain minimum and maximum for entry
        minEntry = min([modelIntrinsic(:,iEntry);estimatedIntrinsic(:,iEntry)]);
        maxEntry = max([modelIntrinsic(:,iEntry);estimatedIntrinsic(:,iEntry)]);
        if abs(minEntry - maxEntry) <= eps && minEntry == 0
            continue
        end
        minPower = ceil(log10(abs(minEntry)) - 1) - 1;
        maxPower = ceil(log10(abs(maxEntry)) - 1) - 1;
        minAxisY = floor(minEntry * 10^(-minPower)) * 10^minPower;
        maxAxisY = ceil(maxEntry * 10^(-maxPower))  * 10^maxPower;

        figure();
        hold on;
        plot(disparities,modelIntrinsic(:,iEntry),'-ok');
        plot(disparities,estimatedIntrinsic(:,iEntry),'-xr');
        if isempty(location{iEntry})
            labelLocation = 'southeast'; 
        else
            labelLocation = location{iEntry}; 
        end
        legend({'Model','Estimate'},'Location',labelLocation);
        ylim([minAxisY,maxAxisY]);
        xlabel('Disparity [pix]');
        ylabel(intrinsicLabels{iEntry})
        axis square;
        grid on;
        utils.vislab.toprint;
        saveas(gcf,[fileparts(resultsFilepath) filesep intrinsicLabels{iEntry} '.png']);
        close(gcf);
    end

    % Fuse entries for display
    intrinsicFusion  = [ 1,4 ; 2,5; 3,6; 7,10; 8,11; 9,12 ];
    location         = {'','','northeast','','','southwest'};
    for iFusion = 1:size(intrinsicFusion,1)
        model     = modelIntrinsic(:,intrinsicFusion(iFusion,:));
        estimated = estimatedIntrinsic(:,intrinsicFusion(iFusion,:));

        % Obtain minimum and maximum for entry
        minEntry = min([model(:);estimated(:)]);
        maxEntry = max([model(:);estimated(:)]);
        if abs(minEntry - maxEntry) <= eps && minEntry == 0
            continue
        end
        minPower = ceil(log10(abs(minEntry)) - 1) - 1;
        maxPower = ceil(log10(abs(maxEntry)) - 1) - 1;
        minAxisY = floor(minEntry * 10^(-minPower)) * 10^minPower;
        maxAxisY = ceil(maxEntry * 10^(-maxPower))  * 10^maxPower;

        figure();
        hold on;
        plot(disparities,model(:,1),'-om');
        plot(disparities,model(:,2),'-oc');
        plot(disparities,estimated(:,1),'-xr');
        plot(disparities,estimated(:,2),'-xb');
        if isempty(location{iFusion})
            labelLocation = 'southeast'; 
        else
            labelLocation = location{iFusion}; 
        end
        legend( { ['Model '    intrinsicLabels{intrinsicFusion(iFusion,1)}] ...
                , ['Model '    intrinsicLabels{intrinsicFusion(iFusion,2)}] ...
                , ['Estimate ' intrinsicLabels{intrinsicFusion(iFusion,1)}] ...
                , ['Estimate ' intrinsicLabels{intrinsicFusion(iFusion,2)}] } ...
              , 'Location',labelLocation );
        ylim([minAxisY,maxAxisY]);
        xlabel('Disparity [pix]');
        axis square;
        grid on;
        utils.vislab.toprint;
        saveas(gcf,[ fileparts(resultsFilepath) filesep ...
                   , intrinsicLabels{intrinsicFusion(iFusion,1)} '_' ...
                   , intrinsicLabels{intrinsicFusion(iFusion,2)} '.png' ]);
        close(gcf);
    end

    % Plot model against estimate for projection matrix
    projectionLabels = { 'k_{11}','k_{22}','k_{13}^0','k_{13}^i','k_{23}^0','k_{23}^j' ...
                       , 't_1^0 [mm]','t_1^i [mm]','t_2^0 [mm]','t_2^j [mm]','t_3 [mm]','k_{13}','k_{23}','t_1 [mm]','t_2 [mm]' };
    location = { '','','','northeast','','northeast' ...
               , '','northeast','','northeast','','','','northeast','northeast' ...
               , 'northeast','northeast','','northeast','northeast','' ...
               , '','northeast','northeast','','','','northeast','','northeast' };
    for iEntry = 1:size(modelProjection,2)
        iLabel = mod(iEntry,15);
        if iLabel == 0; iLabel = 15; end;
        if isempty(strfind(projectionLabels{iLabel},'[mm]'))
            model     = modelProjection(:,iEntry);
            estimated = estimatedProjection(:,iEntry);
            suffix    = '';
        else
            model     = modelProjection(:,iEntry)     * 1000;
            estimated = estimatedProjection(:,iEntry) * 1000;
            suffix    = ' [mm]';
        end

        % Obtain minimum and maximum for entry
        minEntry = min([model;estimated]);
        maxEntry = max([model;estimated]);
        if abs(minEntry - maxEntry) <= eps && minEntry == 0
            continue
        end
        minPower = ceil(log10(abs(minEntry)) - 1) - 1;
        maxPower = ceil(log10(abs(maxEntry)) - 1) - 1;
        if minPower < 0 && maxPower < 0
            minAxisY = floor(minEntry * 10^(-minPower)) * 10^minPower;
            maxAxisY = ceil(maxEntry * 10^(-maxPower))  * 10^maxPower;
        else
            minAxisY = floor(minEntry);
            maxAxisY = ceil(maxEntry);
        end

        figure();
        hold on;
        plot(disparities,model,'-ok');
        plot(disparities,estimated,'-xr');
        if isempty(location{iEntry})
            labelLocation = 'southeast'; 
        else
            labelLocation = location{iEntry}; 
        end
        legend({'Model','Estimate'},'Location',labelLocation);
        ylim([minAxisY,maxAxisY]);
        xlabel('Disparity [pix]');
        if iEntry <= 15
            ylabel([strrep(projectionLabels{iLabel},' [mm]','') '^{ij}' suffix])
        else
            ylabel([strrep(projectionLabels{iLabel},' [mm]','') '^{kl}' suffix])
        end
        axis square;
        grid on;
        if iEntry <= 15; prefix = 'viewpoints_'; else; prefix = 'microlens_'; end;
        utils.vislab.toprint;
        saveas(gcf,[fileparts(resultsFilepath) filesep prefix projectionLabels{iLabel} '.png']);
        close(gcf);
    end

    % Fuse entries for display
    projectionFusion = [ 1,2 ; 3,5; 4,6; 7,9; 8,10; 12,13; 14,15 ];
    projectionFusion = [projectionFusion; 15 + projectionFusion];
    location         = { 'east','','southwest','','northeast','','northeast' ...
                       , 'northeast','northeast','','northeast','','northeast',''};
    for iFusion = 1:size(projectionFusion,1)
        iLabel_1 = mod(projectionFusion(iFusion,1),15);
        iLabel_2 = mod(projectionFusion(iFusion,2),15);
        if iLabel_1 == 0; iLabel_1 = 15; end;
        if iLabel_2 == 0; iLabel_2 = 15; end;

        if projectionFusion(iFusion,1) > 15
            label_1 = [projectionLabels{iLabel_1} '^{kl}'];
            label_2 = [projectionLabels{iLabel_2} '^{kl}'];
        else
            label_1 = [projectionLabels{iLabel_1} '^{ij}'];
            label_2 = [projectionLabels{iLabel_2} '^{ij}'];
        end

        if isempty(strfind(projectionLabels{iLabel_1},'[mm]'))
            labelY    = '';
            model     = modelProjection(:,projectionFusion(iFusion,:));
            estimated = estimatedProjection(:,projectionFusion(iFusion,:));
        else
            label_1   = strrep(label_1,' [mm]','');
            label_2   = strrep(label_2,' [mm]','');
            labelY    = '[mm]';
            model     = modelProjection(:,projectionFusion(iFusion,:))     * 1000;
            estimated = estimatedProjection(:,projectionFusion(iFusion,:)) * 1000;
        end

        % Obtain minimum and maximum for entry
        minEntry = min([model(:);estimated(:)]);
        maxEntry = max([model(:);estimated(:)]);
        if abs(minEntry - maxEntry) <= eps && minEntry == 0
            continue
        end
        minPower = ceil(log10(abs(minEntry)) - 1) - 1;
        maxPower = ceil(log10(abs(maxEntry)) - 1) - 1;
        if minPower < 0 && maxPower < 0
            minAxisY = floor(minEntry * 10^(-minPower)) * 10^minPower;
            maxAxisY = ceil(maxEntry * 10^(-maxPower))  * 10^maxPower;
        else
            minAxisY = floor(minEntry);
            maxAxisY = ceil(maxEntry);
        end

        figure();
        hold on;
        plot(disparities,model(:,1),'-om');
        plot(disparities,model(:,2),'-oc');
        plot(disparities,estimated(:,1),'-xr');
        plot(disparities,estimated(:,2),'-xb');
        if isempty(location{iFusion})
            labelLocation = 'southeast'; 
        else
            labelLocation = location{iFusion}; 
        end
        legend( { ['Model '    label_1] ...
                , ['Model '    label_2] ...
                , ['Estimate ' label_1] ...
                , ['Estimate ' label_2] } ...
              , 'Location',labelLocation);
        ylim([minAxisY,maxAxisY]);
        xlabel('Disparity [pix]');
        if ~isempty(labelY); ylabel(labelY); end;
        axis square;
        grid on;
        if projectionFusion(iFusion,1) <= 15; prefix = 'viewpoints_'; else; prefix = 'microlens_'; end;
        utils.vislab.toprint;
        saveas(gcf,[ fileparts(resultsFilepath) filesep prefix ...
                   , label_1 '_', label_2 '.png' ]);
        close(gcf);
    end
end
