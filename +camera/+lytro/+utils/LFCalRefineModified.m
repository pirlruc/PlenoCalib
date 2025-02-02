% LFCalRefine - refine calibration by minimizing point/ray reprojection error, called by LFUtilCalLensletCam
%
% Usage: 
%     CalOptions = LFCalRefine( InputPath, CalOptions )
% 
% This function is called by LFUtilCalLensletCam to refine an initial camera model and pose
% estimates through optimization. This follows the calibration procedure described in:
%
% D. G. Dansereau, O. Pizarro, and S. B. Williams, "Decoding, calibration and rectification for
% lenslet-based plenoptic cameras," in Computer Vision and Pattern Recognition (CVPR), IEEE
% Conference on. IEEE, Jun 2013.
%
% Minor differences from the paper: camera parameters are automatically initialized, so no prior
% knowledge of the camera's parameters are required; the free intrinsics parameters have been
% reduced by two: H(3:4,5) were previously redundant with the camera's extrinsics, and are now
% automatically centered; and the light field indices [i,j,k,l] are 1-based in this implementation,
% and not 0-based as described in the paper.
% 
% Inputs:
% 
%     InputPath : Path to folder containing decoded checkerboard images. Checkerboard corners must
%                 be identified prior to calling this function, by running LFCalFindCheckerCorners
%                 for example. An initial estiamte must be provided in a CalInfo file, as generated
%                 by LFCalInit. LFUtilCalLensletCam demonstrates the complete procedure.
% 
%     CalOptions struct controls calibration parameters :
%                        .Phase : 'NoDistort' excludes distortion parameters from the optimization
%                                 process; for any other value, distortion parameters are included
%             .CheckerInfoFname : Name of the file containing the summarized checkerboard
%                                 information, as generated by LFCalFindCheckerCorners. Note that
%                                 this parameter is automatically set in the CalOptions struct
%                                 returned by LFCalFindCheckerCorners.
%                 .CalInfoFname : Name of the file containing an initial estimate, to be refined.
%                                 Note that this parameter is automatically set in the CalOptions
%                                 struct returned by LFCalInit.
%          .ExpectedCheckerSize : Number of checkerboard corners, as recognized by the automatic
%                                 corner detector; edge corners are not recognized, so a standard
%                                 8x8-square chess board yields 7x7 corners
%            .LensletBorderSize : Number of pixels to skip around the edges of lenslets, a low
%                                 value of 1 or 0 is generally appropriate
%                   .SaveResult : Set to false to perform a "dry run"
%        [optional]    .OptTolX : Determines when the optimization process terminates. When the
%                                 estimted parameter values change by less than this amount, the
%                                 optimization terminates. See the Matlab documentation on lsqnonlin,
%                                 option `TolX' for more information. The default value of 5e-5 is set
%                                 within the LFCalRefine function; a value of 0 means the optimization
%                                 never terminates based on this criterion.
%      [optional]    .OptTolFun : Similar to OptTolX, except this tolerance deals with the error value.
%                                 This corresponds to Matlab's lsqnonlin option `TolFun'. The default
%                                 value of 0 is set within the LFCalRefine function, and means the
%                                 optimization never terminates based on this criterion.
%
% Outputs :
% 
%     CalOptions struct maintains the fields of the input CalOptions, and adds the fields:
% 
%                    .LFSize : Size of the light field, in samples
%            .IJVecToOptOver : Which samples in i and j were included in the optimization
%           .IntrinsicsToOpt : Which intrinsics were optimized, these are indices into the 5x5
%                              lenslet camera intrinsic matrix
%     .DistortionParamsToOpt : Which distortion params were optimized
%     .PreviousCamIntrinsics : Previous estimate of the camera's intrinsics
%     .PreviousCamDistortion : Previous estimate of the camera's distortion parameters
%                    .NPoses : Number of poses in the dataset
% 
% 
% See also:  LFUtilCalLensletCam, LFCalFindCheckerCorners, LFCalInit, LFUtilDecodeLytroFolder

% Part of LF Toolbox v0.4 released 12-Feb-2015
% Copyright (c) 2013-2015 Donald G. Dansereau

function CalOptions = LFCalRefineModified( InputPath, CalOptions ) 

%---Defaults---
CalOptions = LFDefaultField( 'CalOptions', 'OptTolX', 5e-5 );
CalOptions = LFDefaultField( 'CalOptions', 'OptTolFun', 0 );
CalOptions = LFDefaultField( 'CalOptions', 'DistortModel', 'directions_radialdistort_center' );

%---Load checkerboard corners and previous cal state---
CheckerInfoFname = fullfile(InputPath, CalOptions.CheckerInfoFname);
CalInfoFname = fullfile(InputPath, CalOptions.CalInfoFname);

load(CheckerInfoFname, 'CheckerObs', 'IdealChecker', 'LFSize');
[EstCamPosesV, EstCamIntrinsicsH, EstCamDistortionV, CamInfo, LensletGridModel, DecodeOptions] = ...
    LFStruct2Var( LFReadMetadata(CalInfoFname), 'EstCamPosesV', 'EstCamIntrinsicsH', 'EstCamDistortionV', 'CamInfo', 'LensletGridModel', 'DecodeOptions' );
CalOptions.LFSize = LFSize;

%---Set up optimization variables---
CalOptions.IJVecToOptOver = CalOptions.LensletBorderSize+1:LFSize(1)-CalOptions.LensletBorderSize;
% If (s,t)-plane coincides with viewpoints centers of projection, the 
% pinhole for the viewpoints is already enforced.
% do not optimize h_sk and h_tl
CalOptions.IntrinsicsToOpt = sub2ind([5,5], [1,3, 2,4, 3, 4], [1,1, 2,2, 3, 4]); 				
additionalIntrinsicsToOpt  = sub2ind([5,5], [3,4],[5,5] );
CalOptions.IntrinsicsToOpt = sort([CalOptions.IntrinsicsToOpt, additionalIntrinsicsToOpt]);

switch( lower(CalOptions.Phase) )
    case 'nodistort'
        CalOptions.DistortionParamsToOpt = [];
    case 'withdistort'
        switch( lower(CalOptions.DistortModel) )
            case 'directions_radialdistort_nocenter'                    % Dansereau et al. simplified
                CalOptions.DistortionParamsToOpt = 1:3;
            case 'directions_radialdistort_center'                      % Dansereau et al. original
                CalOptions.DistortionParamsToOpt = 1:5;
            case 'positions_directions_radialdistort_nocenter'          % Bok et al. original
                CalOptions.DistortionParamsToOpt = 1:3;
            case 'positions_directions_radialdistort_center'            % Bok et al. complex
                CalOptions.DistortionParamsToOpt = 1:5;
            case 'directions_distort_center'                            % Johannsen et al.
                if isempty(EstCamDistortionV)
                    CalOptions.DistortionParamsToOpt = [1:3,7,8];
                else
                    CalOptions.DistortionParamsToOpt = 4:6;
                end
            case 'directions_distort_nocenter'                          % Zeller et al.
                CalOptions.DistortionParamsToOpt = 1:6;
            case 'positions_directions_mixture_radialdistort_nocenter'  % Zhang et al. original
                CalOptions.DistortionParamsToOpt = 1:5;
            case 'positions_directions_mixture_radialdistort_center'    % Zhang et al. complex
                CalOptions.DistortionParamsToOpt = 1:7;
            case 'positions_directions_mixture_radialdistort_center_dansereau_tip'    % Monteiro et al.
                CalOptions.DistortionParamsToOpt = 1:9;
                
        end
    otherwise
        switch( lower(CalOptions.DistortModel) )
            case 'directions_radialdistort_nocenter'                    % Dansereau et al. simplified
                CalOptions.DistortionParamsToOpt = 1:3;
            case 'directions_radialdistort_center'                      % Dansereau et al. original
                CalOptions.DistortionParamsToOpt = 1:5;
            case 'positions_directions_radialdistort_nocenter'          % Bok et al. original
                CalOptions.DistortionParamsToOpt = 1:3;
            case 'positions_directions_radialdistort_center'            % Bok et al. complex
                CalOptions.DistortionParamsToOpt = 1:5;
            case 'directions_distort_center'                            % Johannsen et al.
                CalOptions.DistortionParamsToOpt = 1:8; 
            case 'directions_distort_nocenter'                          % Zeller et al.
                CalOptions.DistortionParamsToOpt = 1:6; 
            case 'positions_directions_mixture_radialdistort_nocenter'  % Zhang et al. original
                CalOptions.DistortionParamsToOpt = 1:5;
            case 'positions_directions_mixture_radialdistort_center'    % Zhang et al. complex
                CalOptions.DistortionParamsToOpt = 1:7; 
            case 'positions_directions_mixture_radialdistort_center_dansereau_tip' % Monteiro et al.
                CalOptions.DistortionParamsToOpt = 1:9; 
        end
end

% Initialize distortion parameters
if( isempty(EstCamDistortionV) && ~isempty(CalOptions.DistortionParamsToOpt) )
    switch( lower(CalOptions.DistortModel) )
        case 'directions_radialdistort_nocenter'                    % Dansereau et al. simplified
            EstCamDistortionV = zeros(1,3);
        case 'directions_radialdistort_center'                      % Dansereau et al. original
            EstCamDistortionV = zeros(1,5);
        case 'positions_directions_radialdistort_nocenter'          % Bok et al. original
            EstCamDistortionV = zeros(1,6);
        case 'positions_directions_radialdistort_center'            % Bok et al. complex
            EstCamDistortionV = zeros(1,10);
        case 'directions_distort_center'                            % Johannsen et al.
            EstCamDistortionV = zeros(1,8);
        case 'directions_distort_nocenter'                          % Zeller et al.
            EstCamDistortionV = zeros(1,6);
        case 'positions_directions_mixture_radialdistort_nocenter'  % Zhang et al. original
            EstCamDistortionV = zeros(1,5);
        case 'positions_directions_mixture_radialdistort_center'    % Zhang et al. complex
            EstCamDistortionV = zeros(1,7);
        case 'positions_directions_mixture_radialdistort_center_dansereau_tip'    % Monteiro et al.
            EstCamDistortionV = zeros(1,9);
    end
end
CalOptions.PreviousCamIntrinsics = EstCamIntrinsicsH;
CalOptions.PreviousCamDistortion = EstCamDistortionV;

fprintf('\n===Calibration refinement step, optimizing:===\n');
fprintf('    Intrinsics: ');
disp(CalOptions.IntrinsicsToOpt);
if( ~isempty(CalOptions.DistortionParamsToOpt) )
    fprintf('    Distortion: ');
    disp(CalOptions.DistortionParamsToOpt);
end

%---Compute initial error between projected and measured corner positions---
IdealChecker = [IdealChecker; ones(1,size(IdealChecker,2))]; % homogeneous coord

%---Encode params and grab info required to build Jacobian sparsity matrix---
CalOptions.NPoses = size(EstCamPosesV,1);
[Params0, ParamsInfo, JacobSensitivity] = EncodeParams( EstCamPosesV, EstCamIntrinsicsH, EstCamDistortionV, CalOptions );

[PtPlaneDist0,JacobPattern] = FindError( Params0, CheckerObs, IdealChecker, CalOptions, ParamsInfo, JacobSensitivity );
if( numel(PtPlaneDist0) == 0 )
    error('No valid grid points found -- possible grid parameter mismatch');
end

fprintf('\n    Start SSE: %g m^2, RMSE: %g m\n', sum((PtPlaneDist0).^2), sqrt(mean((PtPlaneDist0).^2)));

%---Start the optimization---
ObjectiveFunc = @(Params) FindError(Params, CheckerObs, IdealChecker, CalOptions, ParamsInfo, JacobSensitivity );
OptimOptions = optimset('Display','iter', ...
    'TolX', CalOptions.OptTolX, ...
    'TolFun',CalOptions.OptTolFun, ...
    'JacobPattern', JacobPattern, ...
    'PlotFcns', @optimplotfirstorderopt, ...
    'UseParallel', 'Always' );
[OptParams, ~, FinalDist] = lsqnonlin(ObjectiveFunc, Params0, [],[], OptimOptions);

%---Decode the resulting parameters and check the final error---
[EstCamPosesV, EstCamIntrinsicsH, EstCamDistortionV] = DecodeParams(OptParams, CalOptions, ParamsInfo);
fprintf(' ---Finished calibration refinement---\n');

fprintf('Estimate of camera intrinsics: \n');
disp(EstCamIntrinsicsH);
if( ~isempty( EstCamDistortionV ) )
    fprintf('Estimate of camera distortion: \n');
    disp(EstCamDistortionV);
end

ReprojectionError = struct( 'SSE', sum(FinalDist.^2), 'RMSE', sqrt(mean(FinalDist.^2)) );
fprintf('\n    Start SSE: %g m^2, RMSE: %g m\n    Finish SSE: %g m^2, RMSE: %g m\n', ...
    sum((PtPlaneDist0).^2), sqrt(mean((PtPlaneDist0).^2)), ...
    ReprojectionError.SSE, ReprojectionError.RMSE );

if( CalOptions.SaveResult )
    TimeStamp = datestr(now,'ddmmmyyyy_HHMMSS');
    GeneratedByInfo = struct('mfilename', mfilename, 'time', TimeStamp, 'VersionStr', LFToolboxVersion);

    SaveFname = fullfile(InputPath, CalOptions.CalInfoFname);
    fprintf('\nSaving to %s\n', SaveFname);
  
    LFWriteMetadata(SaveFname, LFVar2Struct(GeneratedByInfo, LensletGridModel, EstCamIntrinsicsH, EstCamDistortionV, EstCamPosesV, CamInfo, CalOptions, DecodeOptions, ReprojectionError));
end

end

%---------------------------------------------------------------------------------------------------
function [Params0, ParamsInfo, JacobSensitivity] = EncodeParams( EstCamPosesV, EstCamIntrinsicsH, EstCamDistortionV, CalOptions )
% This makes use of FlattenStruct to reversibly flatten all params into a single array.
% It also applies the same process to a sensitivity list, to facilitate building a Jacobian
% Sparisty matrix.

% The 'P' structure contains all the parameters to encode, and the 'J' structure mirrors it exactly
% with a sensitivity list. Each entry in 'J' lists those poses that are senstitive to the
% corresponding parameter. e.g. The first estimated camera pose affects only observations made
% within the first pose, and so the sensitivity list for that parameter lists only the first pose. A
% `J' value of 0 means all poses are sensitive to that variable -- as in the case of the intrinsics,
% which affect all observations.
P.EstCamPosesV = EstCamPosesV;
J.EstCamPosesV = zeros(size(EstCamPosesV));
for( i=1:CalOptions.NPoses )
    J.EstCamPosesV(i,:) = i;
end

P.IntrinParams = EstCamIntrinsicsH(CalOptions.IntrinsicsToOpt);
J.IntrinParams = zeros(size(CalOptions.IntrinsicsToOpt));

P.DistortParams = EstCamDistortionV(CalOptions.DistortionParamsToOpt);
J.DistortParams = zeros(size(CalOptions.DistortionParamsToOpt));

[Params0, ParamsInfo] = FlattenStruct(P);
JacobSensitivity = FlattenStruct(J);
end
%---------------------------------------------------------------------------------------------------
function [EstCamPosesV, EstCamIntrinsicsH, EstCamDistortionV] = DecodeParams( Params, CalOptions, ParamsInfo )
P = UnflattenStruct(Params, ParamsInfo);
EstCamPosesV = P.EstCamPosesV;

EstCamIntrinsicsH = CalOptions.PreviousCamIntrinsics;
EstCamIntrinsicsH(CalOptions.IntrinsicsToOpt) = P.IntrinParams;

EstCamDistortionV = CalOptions.PreviousCamDistortion;
EstCamDistortionV(CalOptions.DistortionParamsToOpt) = P.DistortParams;

intrinsics = false(5,5);
intrinsics(CalOptions.IntrinsicsToOpt) = true;
CalOptions.intrinsicVariables          = intrinsics;
EstCamIntrinsicsH = camera.lytro.utils.LFRecenterIntrinsicsModified(EstCamIntrinsicsH, CalOptions.LFSize, CalOptions);
end

%---------------------------------------------------------------------------------------------------
function [Params, ParamInfo] = FlattenStruct(P)
Params = [];
ParamInfo.FieldNames = fieldnames(P);
for( i=1:length( ParamInfo.FieldNames ) )
    CurFieldName = ParamInfo.FieldNames{i};
    CurField = P.(CurFieldName);
    ParamInfo.SizeInfo{i} = size(CurField);
    Params = [Params; CurField(:)];
end
end
%---------------------------------------------------------------------------------------------------
function [P] = UnflattenStruct(Params, ParamInfo)
CurIdx = 1;
for( i=1:length( ParamInfo.FieldNames ) )
    CurFieldName = ParamInfo.FieldNames{i};
    CurSize = ParamInfo.SizeInfo{i};
    CurField = Params(CurIdx + (0:prod(CurSize)-1));
    CurIdx = CurIdx + prod(CurSize);
    CurField = reshape(CurField, CurSize);
    P.(CurFieldName) = CurField;
end
end

%---------------------------------------------------------------------------------------------------
function [PtPlaneDists, JacobPattern] = FindError(Params, CheckerObs, IdealChecker, CalOptions, ParamsInfo, JacobSensitivity )
    %---Decode optim params---
    [EstCamPosesV, EstCamIntrinsicsH, EstCamDistortionV] = DecodeParams(Params, CalOptions, ParamsInfo);
    
    %---Tally up the total number of observations---
    TotCornerObs = size( [CheckerObs{:,CalOptions.IJVecToOptOver,CalOptions.IJVecToOptOver}], 2 );
    CheckCornerObs = 0;

    %---Preallocate JacobPattern if it's requested---
    if( nargout >= 2 )
        JacobPattern = zeros(TotCornerObs, length(Params));
    end

    %---Preallocate point-plane distances---
    PtPlaneDists     = zeros(1, TotCornerObs);
    
    %---Compute point-plane distances---
    OutputIdx = 0;
    for( PoseIdx = 1:CalOptions.NPoses )
        %---Convert the pertinent camera pose to a homogeneous transform---
        CurEstCamPoseV = squeeze(EstCamPosesV(PoseIdx, :));
        CurEstCamPoseH = eye(4);
        CurEstCamPoseH(1:3,1:3) = rodrigues(CurEstCamPoseV(4:6));
        CurEstCamPoseH(1:3,4) = CurEstCamPoseV(1:3);

        %---Iterate through the corners---
        for( TIdx = CalOptions.IJVecToOptOver )
            for( SIdx = CalOptions.IJVecToOptOver )
                CurCheckerObs = CheckerObs{PoseIdx, TIdx,SIdx};
                NCornerObs = size(CurCheckerObs,2);
                if( NCornerObs ~= prod(CalOptions.ExpectedCheckerSize) )
                    continue; % this implementation skips incomplete observations
                end
                CheckCornerObs = CheckCornerObs + NCornerObs;
                
                %---Assemble observed corner positions into complete 4D [i,j,k,l] indices---
                CurCheckerObs_Idx = [repmat([SIdx;TIdx], 1, NCornerObs); CurCheckerObs; ones(1, NCornerObs)];
                
                %---Transform ideal 3D corner coords into camera's reference frame---
                IdealChecker_CamFrame = CurEstCamPoseH * IdealChecker;
                IdealChecker_CamFrame = IdealChecker_CamFrame(1:3,:); % won't be needing homogeneous points
                
                %---Project observed corner indices to [s,t,u,v] rays---
                CurCheckerObs_Ray = EstCamIntrinsicsH * CurCheckerObs_Idx;
                
                %---Apply direction-dependent distortion model---
                if( ~isempty(EstCamDistortionV) && any(EstCamDistortionV(:)~=0))
                    switch lower(CalOptions.DistortModel)
                        case 'directions_radialdistort_nocenter'    % Dansereau et al. simplified
                            k1 = EstCamDistortionV(1);
                            k2 = EstCamDistortionV(2);
                            k3 = EstCamDistortionV(3);
                            
                            Direction = CurCheckerObs_Ray(3:4,:);
                            DirectionR2 = sum(Direction.^2);
                            Direction = Direction .* repmat((1 + k1.*DirectionR2    + k2.*DirectionR2.^2 ...
                                                               + k3.*DirectionR2.^3 ),2,1);
                            CurCheckerObs_Ray(3:4,:) = Direction;
                            
                        case 'directions_radialdistort_center'      % Dansereau et al. original
                            k1 = EstCamDistortionV(1);
                            k2 = EstCamDistortionV(2);
                            k3 = EstCamDistortionV(3);
                            b1dir = EstCamDistortionV(4);
                            b2dir = EstCamDistortionV(5);
                            
                            Direction = CurCheckerObs_Ray(3:4,:);
                            Direction = bsxfun(@minus, Direction, [b1dir;b2dir]);
                            DirectionR2 = sum(Direction.^2);
                            Direction = Direction .* repmat((1 + k1.*DirectionR2    + k2.*DirectionR2.^2 ...
                                                               + k3.*DirectionR2.^3 ),2,1);
                            Direction = bsxfun(@plus, Direction, [b1dir;b2dir]);
                            CurCheckerObs_Ray(3:4,:) = Direction;
                            
                        case 'positions_directions_radialdistort_nocenter'  % Bok et al. original
                            k1 = EstCamDistortionV(1);
                            k2 = EstCamDistortionV(2);
                            k3 = EstCamDistortionV(3);

                            Direction = CurCheckerObs_Ray(3:4,:);
                            DirectionR2 = sum(Direction.^2);
                            Direction = Direction .* repmat((1 + k1.*DirectionR2    + k2.*DirectionR2.^2 ...
                                                               + k3.*DirectionR2.^3 ),2,1);
                            CurCheckerObs_Ray(3:4,:) = Direction;

                            Position = CurCheckerObs_Ray(1:2,:);
                            PositionR2 = sum(Position.^2);
                            Position = Position .* repmat((1 + k1.*PositionR2    + k2.*PositionR2.^2 ...
                                                             + k3.*PositionR2.^3 ),2,1);
                            CurCheckerObs_Ray(1:2,:) = Position;

                        case 'positions_directions_radialdistort_center'    % Bok et al. complex
                            k1 = EstCamDistortionV(1);
                            k2 = EstCamDistortionV(2);
                            k3 = EstCamDistortionV(3);
                            b1 = EstCamDistortionV(4);
                            b2 = EstCamDistortionV(5);

                            Direction = CurCheckerObs_Ray(3:4,:);
                            Direction = bsxfun(@minus, Direction, [b1;b2]);
                            DirectionR2 = sum(Direction.^2);
                            Direction = Direction .* repmat((1 + k1.*DirectionR2    + k2.*DirectionR2.^2 ...
                                                               + k3.*DirectionR2.^3 ),2,1);
                            Direction = bsxfun(@plus, Direction, [b1;b2]);
                            CurCheckerObs_Ray(3:4,:) = Direction;

                            Position = CurCheckerObs_Ray(1:2,:);
                            Position = bsxfun(@minus, Position, [b1;b2]);
                            PositionR2 = sum(Position.^2);
                            Position = Position .* repmat((1 + k1.*PositionR2    + k2.*PositionR2.^2 ...
                                                             + k3.*PositionR2.^3 ),2,1);
                            Position = bsxfun(@plus, Position, [b1;b2]);
                            CurCheckerObs_Ray(1:2,:) = Position;
                            
                        case 'directions_distort_center'    % Johannsen et al.
                            k1 = EstCamDistortionV(1);
                            k2 = EstCamDistortionV(2);
                            k3 = EstCamDistortionV(3);
                            p1 = EstCamDistortionV(4);
                            p2 = EstCamDistortionV(5);
                            p3 = EstCamDistortionV(6);
                            b1dir = EstCamDistortionV(7);
                            b2dir = EstCamDistortionV(8);
                            
                            Direction = CurCheckerObs_Ray(3:4,:);
                            Direction0  = bsxfun(@minus, Direction, [b1dir;b2dir]);
                            DirectionR2 = sum(Direction0.^2);
                            Direction_radial     = Direction0 .* repmat((1 + k1.*DirectionR2    + k2.*DirectionR2.^2 ...
                                                                           + k3.*DirectionR2.^3 ),2,1);
                            Direction_tangential = (      [p1;p2] .* (repmat(DirectionR2,2,1) + 2 * Direction0.^2) ...
                                                   + 2 .* [p2;p1] .* repmat(Direction0(1,:) .* Direction(2,:),2,1) ) ...
                                                .* ( 1 + p3 .* repmat(DirectionR2,2,1) );
                            Direction = Direction_radial + Direction_tangential;
                            Direction = bsxfun(@plus, Direction, [b1dir;b2dir]);
                            CurCheckerObs_Ray(3:4,:) = Direction;
                        
                        case 'directions_distort_nocenter'  % Zeller et al.
                            k1 = EstCamDistortionV(1);
                            k2 = EstCamDistortionV(2);
                            k3 = EstCamDistortionV(3);
                            p1 = EstCamDistortionV(4);
                            p2 = EstCamDistortionV(5);
                            p3 = EstCamDistortionV(6);
                            
                            Direction0  = CurCheckerObs_Ray(3:4,:);
                            DirectionR2 = sum(Direction0.^2);
                            Direction_radial     = Direction0 .* repmat((1 + k1.*DirectionR2    + k2.*DirectionR2.^2 ...
                                                                           + k3.*DirectionR2.^3 ),2,1);
                            Direction_tangential = (      [p1;p2] .* (repmat(DirectionR2,2,1) + 2 * Direction0.^2) ...
                                                   + 2 .* [p2;p1] .* repmat(Direction0(1,:) .* Direction0(2,:),2,1) ) ...
                                                .* ( 1 + p3 .* repmat(DirectionR2,2,1) );
                            Direction = Direction_radial + Direction_tangential;
                            CurCheckerObs_Ray(3:4,:) = Direction;

                        case 'positions_directions_mixture_radialdistort_nocenter'  % Zhang et al. original
                            k1 = EstCamDistortionV(1);
                            k2 = EstCamDistortionV(2);
                            k3 = EstCamDistortionV(3);
                            s1 = EstCamDistortionV(4);
                            t1 = EstCamDistortionV(5);

                            Position    = CurCheckerObs_Ray(1:2,:);
                            Direction   = CurCheckerObs_Ray(3:4,:);
                            DirectionR2 = sum(Direction.^2);
                            Direction   = Direction .* repmat((1 + k1.*DirectionR2    + k2.*DirectionR2.^2 ...
                                                                 + k3.*DirectionR2.^3 ),2,1) ...
                                        + [s1;t1] .* Position;
                            CurCheckerObs_Ray(3:4,:) = Direction;
                            
                        case 'positions_directions_mixture_radialdistort_center'    % Zhang et al. complex
                            k1 = EstCamDistortionV(1);
                            k2 = EstCamDistortionV(2);
                            k3 = EstCamDistortionV(3);
                            s1 = EstCamDistortionV(4);
                            t1 = EstCamDistortionV(5);
                            b1dir = EstCamDistortionV(6);
                            b2dir = EstCamDistortionV(7);
                            
                            Position    = CurCheckerObs_Ray(1:2,:);
                            Direction   = CurCheckerObs_Ray(3:4,:);
                            Direction   = bsxfun(@minus, Direction, [b1dir;b2dir]);
                            DirectionR2 = sum(Direction.^2);
                            Direction   = Direction .* repmat((1 + k1.*DirectionR2    + k2.*DirectionR2.^2 ...
                                                                 + k3.*DirectionR2.^3 ),2,1) ...
                                        + [s1;t1] .* Position;
                            Direction   = bsxfun(@plus, Direction, [b1dir;b2dir]);
                            CurCheckerObs_Ray(3:4,:) = Direction;

						case 'positions_directions_mixture_radialdistort_center_dansereau_tip'    % Monteiro et al.
                            k1 = EstCamDistortionV(1);
                            k2 = EstCamDistortionV(2);
                            k3 = EstCamDistortionV(3);
                            b1dir = EstCamDistortionV(4);
                            b2dir = EstCamDistortionV(5);
                            s1 = EstCamDistortionV(6);
                            t1 = EstCamDistortionV(7);
                            s2 = EstCamDistortionV(8);
                            t2 = EstCamDistortionV(9);
                            
                            Position    = CurCheckerObs_Ray(1:2,:);
                            Direction   = CurCheckerObs_Ray(3:4,:);
                            Direction   = bsxfun(@minus, Direction, [b1dir;b2dir]);
                            DirectionR2 = sum(Direction.^2);
                            Direction   = Direction .* repmat((1 + k1.*DirectionR2    + k2.*DirectionR2.^2 ...
                                                                 + k3.*DirectionR2.^3 ),2,1) ...
                                        + [s1;t1] .* Position + [s2;t2] .* Position.^2;
                            Direction   = bsxfun(@plus, Direction, [b1dir;b2dir]);
                            CurCheckerObs_Ray(3:4,:) = Direction;
                    end
                end
                
                %---Find 3D point-ray distance---
                STPlaneIntersect = [CurCheckerObs_Ray(1:2,:); zeros(1,NCornerObs)];
                RayDir = [CurCheckerObs_Ray(3:4,:); ones(1,NCornerObs)];
                
				CurDist3D = LFFind3DPtRayDist( STPlaneIntersect, RayDir, IdealChecker_CamFrame );
				PtPlaneDists(OutputIdx + (1:NCornerObs))     = CurDist3D;
                if( nargout >=2 )
                    % Build the Jacobian pattern. First we enumerate those observations related to
                    % the current pose, then find all parameters to which those observations are
                    % sensitive. This relies on the JacobSensitivity list constructed by the
                    % FlattenStruct function.
                    CurObservationList = OutputIdx + (1:NCornerObs);
                    CurSensitivityList = (JacobSensitivity==PoseIdx | JacobSensitivity==0);
                    JacobPattern(CurObservationList, CurSensitivityList) = 1;
                end
                OutputIdx = OutputIdx + NCornerObs;
            end
        end
    end
    
    %---Check that the expected number of observations have gone by---
    if( CheckCornerObs ~= TotCornerObs )
        error(['Mismatch between expected (%d) and observed (%d) number of corners' ...
            ' -- possibly caused by a grid parameter mismatch'], TotCornerObs, CheckCornerObs);
    end
end

%---Compute distances from 3D rays to a 3D points---
function [Dist] = LFFind3DPtRayDist( PtOnRay, RayDir, Pt3D )

RayDir = RayDir ./ repmat(sqrt(sum(RayDir.^2)), 3,1); % normalize ray
Pt3D = Pt3D - PtOnRay;    % Vector to point

PD1 = dot(Pt3D, RayDir);  
PD1 = repmat(PD1,3,1).*RayDir; % Project point vec onto ray vec

Pt3D = Pt3D - PD1; 
Dist = sqrt(sum(Pt3D.^2, 1)); % Distance from point to projected point

end
