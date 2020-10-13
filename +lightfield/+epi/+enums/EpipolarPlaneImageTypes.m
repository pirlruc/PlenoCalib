classdef EpipolarPlaneImageTypes
    %EPIPOLARPLANEIMAGETYPES
    %   Epipolar plane image types available for lightfield.
    
    properties (~Dependent, SetAccess = immutable, GetAccess = public)
        symbol              % Epipolar plane image symbol
        fixedCoordinates    % Image ray fixed coordinates
        description         % Description of the epipolar plane images
        encodingOrder       % Ordering to convert lightfield to epipolar plane image
        decodingOrder       % Ordering to convert epipolar plane image to lightfield
    end
    
    methods
        function self = EpipolarPlaneImageTypes( symbol, coordinates, description ...
                                               , encodingOrder, decodingOrder )
            %
            % Derivative masks enumeration instance. Notice that the masks
            % are normalized using L1-norm.
            %
            self.symbol           = symbol;
            self.fixedCoordinates = coordinates;
            self.description      = description;
            self.encodingOrder    = encodingOrder;
            self.decodingOrder    = decodingOrder;
        end
    end
    
    methods 
        function decodedData = decode(self,encodedData)
            %
            % Decode epipolar plane image data to internal format.
            %
            % INPUTS:
            %   1. encodedData - data with epipolar plane image format.
            %
            narginchk(2,2);
            
            % If the decoding order does not have the same dimension,
            % assume that the remaining dimensions will remain in the same
            % position
            if numel(size(encodedData)) > numel(self.decodingOrder)
                decoding = [self.decodingOrder,numel(self.decodingOrder) + 1:numel(size(encodedData))];
            else
                decoding = self.decodingOrder;
            end
            
            decodedData = permute(encodedData,decoding);
        end
        
        function encodedData = encode(self,decodedData)
            %
            % Encode internal format to epipolar plane image format.
            %
            % INPUTS:
            %   1. encodedData - data with lightfield internal format.
            %
            narginchk(2,2);
            
            % If the encoding order does not have the same dimension,
            % assume that the remaining dimensions will remain in the same
            % position
            if numel(size(decodedData)) > numel(self.encodingOrder)
                encoding = [self.encodingOrder,numel(self.encodingOrder) + 1:numel(size(decodedData))];
            else
                encoding = self.encodingOrder;
            end
            
            encodedData = permute(decodedData,encoding);
        end
    end
    
    enumeration
        EPI_IK ( 'ik', '(j,l)' ...
               , [ '(i,k) epipolar plane image fixing coordinates (j,l). ' ...
                 , 'Format: channels x microlens_k x pixel_i x pixel j x microlens_l.' ] ...
               , [1,4,2,3,5], [1,3,4,2,5] )
        EPI_JL ('jl', '(i,k)' ...
               , [ '(j,l) epipolar plane image fixing coordinates (i,k). ' ...
                 , 'Format: channels x microlens_l x pixel_j x pixel i x microlens_k.' ] ...
               , [1,5,3,2,4], [1,4,3,5,2] )
    end
end