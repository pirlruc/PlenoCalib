classdef MicrolensImage < abstract.TemplateImage
    %MICROLENSIMAGE
    %   Microlens image methods and properties.

    properties (Constant)
        % The microlens image should have the format:
        %       channels x pixel i x pixel j
        ENCODING_ORDER = [1,2,3,4,5]    % Order to encode lightfield into microlens image format
        DECODING_ORDER = [1,2,3,4,5]    % Order to decode microlens image format into lightfield
    end
    
    properties (Dependent)
        size            % Microlens image size
    end

    methods
        function self = MicrolensImage(varargin)
            %
            % Microlens image instance.
            %
            % INPUTS:
            %   1. microlensImageData - microlens image data.
            %
            narginchk(0,1);
            
            % Create super class instance
            self = self@abstract.TemplateImage();
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
        
        function microlensSize = get.size(self)
            % If microlens image data is empty, consider the default
            % microlens image size values.
            if isempty(self.data)
                microlensSize = lightfield.image.MicrolensImageSize();
            else
                microlensSize = lightfield.image.MicrolensImageSize(size(self.data));
            end
        end
    end
    
    methods (Static)
        function self = MicrolensImageFromLightfield(lightfieldInstance, microlens_kl)
            %
            % Create microlens image instance from lightfield data.
            %
            % INPUTS:
            %   1. lightfieldInstance - lightfield instance.
            %   2. microlens_kl       - microlens (k,l) in lightfield.
            %   Corresponds to the pixels (i,j) under microlens (k,l).
            %
            narginchk(1,2);
            
            % Convert input data to instances
            lightfieldInstance = utils.enums.Classes.LIGHTFIELD().convert(lightfieldInstance);

            % Define microlens position
            if nargin <= 1
                microlens_kl = [ 1; 1 ...
                               ; round(lightfieldInstance.size.numberMicrolenses_k / 2) ...
                               ; round(lightfieldInstance.size.numberMicrolenses_l / 2) ];
            end
            
            % Convert microlens (k,l) to image ray
            microlens_kl = utils.enums.Classes.IMAGE_RAY().convert(microlens_kl);
            
            % Obtain microlens image from lightfield. The microlens image
            % should have the format:
            %       channels x pixel i x pixel j
            microlensImageData = permute( lightfieldInstance.data(:,:,:,microlens_kl.k,microlens_kl.l) ...
                                        , lightfield.image.MicrolensImage.ENCODING_ORDER );
            
            % Create image instance
            self = lightfield.image.MicrolensImage(microlensImageData);
        end
        
        function MicrolensImageMontageFromLightfield(lightfieldInstance)
            %
            % Create montage of lightfield using microlens images.
            %
            % INPUTS:
            %   1. lightfieldInstance - lightfield instance.
            %
            narginchk(1,1);
            
            % Convert input data to instances
            lightfieldInstance = utils.enums.Classes.LIGHTFIELD().convert(lightfieldInstance);

            % Obtain matlab internal format for images
            matlabFormat = image.enums.Formats.MATLAB_FORMAT();

            % Create image sequence with microlens images from lightfield.
            % The sequence will be represented as an array of k x l images
            % with pixels j x i.
            imageSequence = [];
            for microlens_l = 1:lightfieldInstance.size.numberMicrolenses_l
                microlensImages = lightfieldInstance.data(:,:,:,:,microlens_l);
                imageSequence   = cat(4,imageSequence,matlabFormat.encode(microlensImages));
            end
                
            % Display image sequence using montage
            montage(imageSequence,'Size', [ lightfieldInstance.size.numberMicrolenses_l ...
                                          , lightfieldInstance.size.numberMicrolenses_k ]);
            set(gca,'position',[0 0 1 1],'units','normalized')
            axis tight; 
            axis off;
        end
    end
end