classdef ViewpointImage < abstract.TemplateImage
    %VIEWPOINTIMAGE
    %   Viewpoint image methods and properties.

    properties (Constant)
        % The viewpoint image should have the format:
        %       channels x microlens k x microlens l
        ENCODING_ORDER = [1,4,5,2,3]    % Order to encode lightfield into viewpoint image format
        DECODING_ORDER = [1,4,5,2,3]    % Order to decode viewpoint image format into lightfield
    end
    
    properties (Dependent)
        size            % Viewpoint image size
    end

    methods
        function self = ViewpointImage(varargin)
            %
            % Viewpoint image instance
            %
            % INPUTS:
            %   1. viewpointImageData - viewpoint image data.
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
        
        function viewpointSize = get.size(self)
            % If microlens image data is empty, consider the default
            % microlens image size values.
            if isempty(self.data)
                viewpointSize = lightfield.image.ViewpointImageSize();
            else
                viewpointSize = lightfield.image.ViewpointImageSize(size(self.data));
            end
        end
    end
    
    methods (Static)
        function self = ViewpointImageFromLightfield(lightfieldInstance, viewpoint_ij)
            %
            % Create viewpoint image instance from lightfield data.
            %
            % INPUTS:
            %   1. lightfieldInstance - lightfield instance.
            %   2. viewpoint_ij       - viewpoint (i,j) in lightfield.
            %   Corresponds to the pixels (i,j) under each microlens (k,l).
            %
            narginchk(1,2);
            
            % Convert input data to instances
            lightfieldInstance = utils.enums.Classes.LIGHTFIELD().convert(lightfieldInstance);

            % Define viewpoint position
            if nargin <= 1
                viewpoint_ij = [ round(lightfieldInstance.size.numberPixels_i / 2) ...
                               ; round(lightfieldInstance.size.numberPixels_j / 2) ...
                               ; 1; 1];
            end
            
            % Convert viewpoint (i,j) to image ray
            viewpoint_ij = utils.enums.Classes.IMAGE_RAY().convert(viewpoint_ij);
            
            % Obtain viewpoint image from lightfield. The viewpoint image
            % should have the format:
            %       channels x microlens k x microlens l
            viewpointImageData = lightfieldInstance.data(:,viewpoint_ij.i,viewpoint_ij.j,:,:);
            viewpointImageData = permute(viewpointImageData,lightfield.image.ViewpointImage.ENCODING_ORDER);
            
            % Create image instance
            self = lightfield.image.ViewpointImage(viewpointImageData);
        end
        
        function ViewpointImageMontageFromLightfield(lightfieldInstance)
            %
            % Create montage of lightfield using viewpoint images.
            %
            % INPUTS:
            %   1. lightfieldInstance - lightfield instance.
            %
            narginchk(1,1);
            
            % Convert input data to instances
            lightfieldInstance = utils.enums.Classes.LIGHTFIELD().convert(lightfieldInstance);

            % Obtain matlab internal format for images
            matlabFormat = image.enums.Formats.MATLAB_FORMAT();

            % Create image sequence with viewpoint images from lightfield
            % The sequence will be represented as an array of i x j images
            % with pixels l x k.
            imageSequence = [];
            for viewpoint_i = 1:lightfieldInstance.size.numberPixels_i
                viewpointImages = permute(lightfieldInstance.data(:,viewpoint_i,:,:,:),[1,4,5,3,2]);
                imageSequence = cat(4,imageSequence,matlabFormat.encode(viewpointImages));
            end

            % Display image sequence using montage
            montage(imageSequence,'Size', [ lightfieldInstance.size.numberPixels_i ...
                                          , lightfieldInstance.size.numberPixels_j ]);
            set(gca,'position',[0 0 1 1],'units','normalized')
            axis tight;
            axis off;
        end
    end
end