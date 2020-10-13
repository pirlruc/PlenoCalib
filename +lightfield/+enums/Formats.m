classdef Formats
    %FORMATS
    %   Package to define lightfield formats available.
    
    properties (SetAccess = immutable)
        formattingType          % Formatting type for files in dataset
        description             % Description of formatting type for files in dataset
        decodingOrder           % Order to decode data from lightfield datasets
        encodingOrder           % Order to encode data from lightfield datasets
        pixelDimension_i        % Dimension of the microlens pixel i information
        pixelDimension_j        % Dimension of the microlens pixel j information
        microlensDimension_k    % Dimension of the microlens k information
        microlensDimension_l    % Dimension of the microlens l information
        channelDimension        % Dimension of the channel information
    end
   
    methods
        function self = Formats( formattingType,description ...
                               , decodingOrder,encodingOrder ...
                               , pixelDimension_i,pixelDimension_j,microlensDimension_k,microlensDimension_l ...
                               , channelDimension)
            %
            % Lightfield formats instance.
            %
            self.formattingType   = formattingType;
            self.description      = description;
            self.decodingOrder    = decodingOrder;
            self.encodingOrder    = encodingOrder;
            self.pixelDimension_i = pixelDimension_i;
            self.pixelDimension_j = pixelDimension_j;
            self.channelDimension = channelDimension;
            self.microlensDimension_k = microlensDimension_k;
            self.microlensDimension_l = microlensDimension_l;
        end
    end
    
    methods
        function decodedData = decode(self,encodedData)
            %
            % Decode lightfield data to internal structure for lightfield
            % data.
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
            % Encode lightfield data from internal structure to a specific
            % lightfield format.
            %
            % INPUTS:
            %   1. encodedData - data encoded with a given format.
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
        HEIDELBERG_BENCHMARK_FORMAT( 'heidelberg_benchmark' ...
                                   , 'pixel j x pixel i x microlens l x microlens k x channels'...
                                   , [5,2,1,4,3] ...
                                   , [3,2,5,4,1] ...
                                   , 1,2,4,3,5 )
        DANSEREAU_FORMAT ( 'dansereau' ...
                         , 'pixel j x pixel i x microlens l x microlens k x channels'...
                         , [5,2,1,4,3] ...
                         , [3,2,5,4,1] ...
                         , 2,1,4,3,5 )
        HEIDELBERG_FORMAT( 'heidelberg' ...
                         , 'channels x microlens k x microlens l x pixel i x pixel j'...
                         , [1,4,5,2,3] ...
                         , [1,4,5,2,3] ...
                         , 4,5,2,3,1 )
        INTERNAL_FORMAT  ( 'internal' ...
                         , 'channels x pixel i x pixel j x microlens k x microlens l' ...
                         , [1,2,3,4,5] ...
                         , [1,2,3,4,5] ...
                         , 2,3,4,5,1 )
    end
end