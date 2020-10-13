% LFRecenterIntrinsics - Recenter a light field intrinsic matrix
%
% Usage:
%     H = LFRecenterIntrinsics( H, LFSize )
%
% The recentering works by forcing the central sample in a light field of LFSize samples to
% correspond to the ray [s,t,u,v] = 0. Note that 1-based indexing is assumed in [i,j,k,l].

% Part of LF Toolbox v0.4 released 12-Feb-2015
% Copyright (c) 2013-2015 Donald G. Dansereau

function H = LFRecenterIntrinsicsModified( H, LFSize, CalOptions )
CalOptions = LFDefaultField( 'CalOptions', 'cameraCoordinateOriginRay', (LFSize([2,1,4,3])-1)/2 + 1);
if size(CalOptions.cameraCoordinateOriginRay,1) ~= 1
    CalOptions.cameraCoordinateOriginRay = CalOptions.cameraCoordinateOriginRay';
end

% Store previous H to check changed indices
previousH = H;

% If (s,t)-plane coincides with viewpoints centers of projection, the 
% pinhole for the viewpoints is already enforced.
% No h_sk and h_tl
H(1,3) = 0;
H(2,4) = 0;

% Check if all changed indices are not being optimized
indices = abs(previousH - H) > 0;
if any(CalOptions.intrinsicVariables(indices) == true)
    error = MException( 'LFRecenterIntrinsicsModified:changedEntriesBeingOptimized' ...
                      , 'Entry being optimized changed during optimization.' );
    error.throw();
end

CenterRay = [CalOptions.cameraCoordinateOriginRay, 1]'; % note indices start at 1
H(1:2,5) = 0;
Decentering = H * CenterRay;
H(1:2,5) = -Decentering(1:2);
end
