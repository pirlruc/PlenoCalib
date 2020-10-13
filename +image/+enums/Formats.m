classdef Formats
    %FORMATS
    %   Package to define image formats available.
    
    properties (SetAccess = immutable)
        formattingType      % Formatting type for files in dataset
        description         % Description of formatting type for files in dataset
        decodingOrder       % Order to decode data from image datasets
        encodingOrder       % Order to encode data from image datasets
        pixelDimension_u    % Dimension of the pixel u information
        pixelDimension_v    % Dimension of the pixel v information
        channelDimension    % Dimension of the channel information
    end
   
    methods
        function self = Formats( formattingType,description ...
                               , decodingOrder,encodingOrder ...
                               , pixelDimension_u,pixelDimension_v,channelDimension)
            %
            % Image formats instance.
            %
            self.formattingType   = formattingType;
            self.description      = description;
            self.decodingOrder    = decodingOrder;
            self.encodingOrder    = encodingOrder;
            self.pixelDimension_u = pixelDimension_u;
            self.pixelDimension_v = pixelDimension_v;
            self.channelDimension = channelDimension;
        end
    end
    
    methods
        function decodedData = decode(self,encodedData)
            %
            % Decode image data to internal structure for image data.
            %
            % INPUTS:
            %   1. encodedData - data encoded with a given format.
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
            % Encode image data from internal structure to a specific
            % image format.
            %
            % INPUTS:
            %   1. decodedData - data encoded in internal format.
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
        MATLAB_FORMAT    ( 'matlab' ...
                         , 'pixel v x pixel u x channels'...
                         , [3,2,1] ...
                         , [3,2,1] ...
                         , 2,1,3 )
        INTERNAL_FORMAT  ( 'internal' ...
                         , 'channels x pixel u x pixel v' ...
                         , [1,2,3] ...
                         , [1,2,3] ...
                         , 2,3,1 )
    end
end