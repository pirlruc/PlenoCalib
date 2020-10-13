classdef EpipolarPlaneImage < abstract.TemplateEpipolarPlaneImage
    %EPIPOLARPLANEIMAGE
    %   Epipolar plane image for pairs (i,k) or (j,l) coordinates of
    %   lightfield.
    
    properties (Dependent)
        size            % Epipolar plane image size for (i,k) or (j,l) pair.
    end

    methods
        function self = EpipolarPlaneImage(varargin)
            %
            % Create epipolar plane image instance for (i,k) or (j,l) pair.
            %
            % INPUTS:
            %   1. epiData - epipolar plane image data.
            %
            narginchk(0,1);
            
            % Create super class instance.
            self = self@abstract.TemplateEpipolarPlaneImage();
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
        
        function epiSize = get.size(self)
            % If epipolar plane image data is empty, consider the default 
            % epipolar plane image size values.
            if isempty(self.data)
                epiSize = lightfield.epi.EpipolarPlaneImageSize();
            else
                epiSize = lightfield.epi.EpipolarPlaneImageSize(size(self.data));
            end
        end
    end

    methods (Static)
        function self = EpipolarPlaneImageFromLightfield( lightfieldInstance ...
                                                        , epipolarPlaneImageType ...
                                                        , imageRay_ijkl )
            %
            % Create epipolar plane image instance from lightfield data.
            %
            % INPUTS:
            %   1. lightfieldInstance - lightfield instance.
            %   2. epipolarPlaneImageType - type of epipolar plane image.
            %   The types available are (i,k) or (j,l). The default is
            %   (i,k).
            %   3. imageRay_ijkl - image ray coordinates to be considered
            %   for retrieving the epipolar plane image data. Use either
            %   (i,k) or (j,l) pair of coordinates.
            %
            narginchk(1,3);

            % Define default epipolar plane image type
            if nargin <= 1
                epipolarPlaneImageType = lightfield.epi.enums.EpipolarPlaneImageTypes.EPI_IK();
            end
            
            % Convert lightfield data to instance
            lightfieldInstance = utils.enums.Classes.LIGHTFIELD().convert(lightfieldInstance);

            % Define default coordinate pairs for epipolar plane image
            if nargin <= 2
                if epipolarPlaneImageType == lightfield.epi.enums.EpipolarPlaneImageTypes.EPI_IK
                    imageRay_ijkl = [ 1, round(lightfieldInstance.size.numberPixels_j / 2) ...
                                    , 1, round(lightfieldInstance.size.numberMicrolenses_l / 2) ];
                else
                    imageRay_ijkl = [ round(lightfieldInstance.size.numberPixels_i / 2), 1 ...
                                    , round(lightfieldInstance.size.numberMicrolenses_k / 2), 1 ];
                end
            end
            
            % Convert image ray to instance
            imageRay_ijkl = utils.enums.Classes.IMAGE_RAY().convert(imageRay_ijkl);
            
            % Obtain epipolar plane image from lightfield. The epipolar
            % plane image should have the format:
            %       channels x microlens k x pixel i
            %   or  
            %       channels x microlens l x pixel j
            if epipolarPlaneImageType == lightfield.epi.enums.EpipolarPlaneImageTypes.EPI_IK
                epipolarPlaneImageData = lightfieldInstance.data(:,:,imageRay_ijkl.j,:,imageRay_ijkl.l);
            else
                epipolarPlaneImageData = lightfieldInstance.data(:,imageRay_ijkl.i,:,imageRay_ijkl.k,:);
            end                
            epipolarPlaneImageData = permute(epipolarPlaneImageData,epipolarPlaneImageType.encodingOrder);
            
            % Create epipolar plane image instance
            self = lightfield.epi.EpipolarPlaneImage(epipolarPlaneImageData);
        end
        
        function EpipolarPlaneImageMontageFromLightfield( lightfieldInstance ...
                                                        , imageRay_ijkl ...
                                                        , neighborhood_ijkl )
            %
            % Create epipolar plane image instance from lightfield data.
            %
            % INPUTS:
            %   1. lightfieldInstance - lightfield instance.
            %   2. imageRay_ijkl - image ray coordinates to be considered
            %   for retrieving the epipolar plane image data.
            %
            narginchk(1,3);

            % Convert lightfield data to instance
            lightfieldInstance = utils.enums.Classes.LIGHTFIELD().convert(lightfieldInstance);

            % Define default coordinates for epipolar plane image montage
            if nargin <= 1
                imageRay_ijkl = [ round(lightfieldInstance.size.numberPixels_i / 2) ...
                                , round(lightfieldInstance.size.numberPixels_j / 2) ...
                                , round(lightfieldInstance.size.numberMicrolenses_k / 2) ...
                                , round(lightfieldInstance.size.numberMicrolenses_l / 2) ];
            end
            
            if nargin <= 2
                neighborhood_ijkl = [0,0,1,1];
            end
            
            % Convert image ray to instance
            imageRay_ijkl = utils.enums.Classes.IMAGE_RAY().convert(imageRay_ijkl);
            
            % Convert image ray to instance
            neighborhood_ijkl = utils.enums.Classes.IMAGE_RAY().convert(neighborhood_ijkl);

            % Obtain epipolar plane image from lightfield
            % The EPI has the format:   microlens k x pixel i.
            % The montage has the format: 
            %           microlens k x pixel i x channel x frames
            imageSequence = [];
            for viewpoint_j = (imageRay_ijkl.j - neighborhood_ijkl.j):(imageRay_ijkl.j + neighborhood_ijkl.j)
                epipolarPlaneImages = lightfieldInstance.data( :,:,viewpoint_j,: ...
                                                             ,   imageRay_ijkl.l - neighborhood_ijkl.l ...
                                                               : imageRay_ijkl.l + neighborhood_ijkl.l);
                imageSequence = cat(4,imageSequence,permute(epipolarPlaneImages,[2,4,1,5,3]));
            end

            % Display image sequence using montage
            figure(1)
            montage(imageSequence)
            set(gca,'position',[0 0 1 1],'units','normalized')
            axis tight; 
            axis off;

            % The EPI has the format:   microlens l x pixel j.
            % The montage has the format: 
            %           microlens l x pixel j x channel x frames
            imageSequence = [];
            for viewpoint_i = (imageRay_ijkl.i - neighborhood_ijkl.i):(imageRay_ijkl.i + neighborhood_ijkl.i)
                epipolarPlaneImages = lightfieldInstance.data( :,viewpoint_i,: ...
                                                             ,   imageRay_ijkl.k - neighborhood_ijkl.k ...
                                                               : imageRay_ijkl.k + neighborhood_ijkl.k, : );
                imageSequence = cat(4,imageSequence,permute(epipolarPlaneImages,[3,5,1,4,2]));
            end                
            
            % Display image sequence using montage
            figure(2)
            montage(imageSequence);
            set(gca,'position',[0 0 1 1],'units','normalized')
            axis tight; 
            axis off;
        end
    end
end