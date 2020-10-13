classdef ColorCorrectionMethods
    %COLORCORRECTIONMETHODS
    %   Color correction methods.
    
    properties (~Dependent, SetAccess = immutable, GetAccess = public)
        symbol          % Color correction method symbol
        method          % Method to compute color correction method
        description     % Description of color correction method
    end
    
    methods
        function self = ColorCorrectionMethods(symbol,method,description)
            %
            % Color correction method enumeration instance.
            %
            self.symbol = symbol;
            self.method = method;
            self.description = description;
        end
    end
    
    methods
        function correctedImage = correct(self,imageData)
            %
            % Correct color image.
            %
            % INPUTS:
            %   1. imageData - image data.
            %
            narginchk(2,2);
            
            % Correct color image.
            correctedImage = self.method(imageData);
        end
    end
    
    methods (Static)
        function correctedImage = shadesOfGrayNormalization(imageData,pNorm)
            %
            % Correct image colors by normalizing each color channel of the
            % image. For more information please refer to:
            %   Barata, Catarina, M. Emre Celebi, and Jorge S. Marques. 
            %       "Improving dermoscopy image classification using color 
            %       constancy." IEEE journal of biomedical and health 
            %       informatics 19.3 (2015): 1146-1152.            
            %
            % INPUTS:
            %   1. imageData - image data.
            %   2. pNorm - norm that should be considered for estimating
            %   the illuminant. Default is 6.
            %
            narginchk(1,2);
            
            % Define default correction method.
            %       p = 1   -> Gray world
            %       p = inf -> max-RGB
            if nargin <= 1
                pNorm = 6;
            end
            
            % Transform to image instance
            imageData = utils.enums.Classes.IMAGE().convert(imageData);

            % Obtain weight for each color channel. The weight corresponds
            % to the Minkowski norm with the degree given by the p-norm.
            if pNorm == 1           % Gray World
                colorVector = sum( sum( abs(imageData.data),3 ),2 ) ...                         % Integration on image pixels
                           ./ (imageData.size.numberPixels_u * imageData.size.numberPixels_v);  % Number of pixels in channel
            elseif pNorm == Inf     % Max-RGB
                colorVector = max( max( abs(imageData.data), [],3 ), [],2 );
            else
                colorVector = (    sum( sum( abs(imageData.data) .^ pNorm,3 ),2 ) ...                  % Integration on image pixels
                                ./ (imageData.size.numberPixels_u * imageData.size.numberPixels_v) ... % Number of pixels in channel
                              ).^(1/pNorm);
            end
            
            % Normalize the color vector considering an l2-norm.
            colorVector = colorVector ./ sqrt( sum( colorVector.^2 ) );
            
            % Obtain new image with corrected colors using the von Kris
            % diagonal model
            correctedImage = nan(imageData.size.data);

            % Obtain weight to correct channel
            channelCorrectionWeight = (1/sqrt(3)) ./ colorVector;
            for iChannel = 1:imageData.size.numberChannels
                % Obtain normalized color channel
                correctedImage(iChannel,:,:) = channelCorrectionWeight(iChannel) .* imageData.data(iChannel,:,:);
            end
            
            % Transform to image
            correctedImage = image.Image(correctedImage);
        end
        
        function correctedImage = meanVarianceNormalization(imageData)
            %
            % Correct image colors by setting the mean value of the image
            % to the middle value of the range defined for the image type.
            % The variance is considered to be unitary.
            %
            % This method corrects illumination of the images and does not
            % treat image channels differently.
            %
            % INPUTS:
            %   1. imageData - image data.
            %
            narginchk(1,1);
            
            RANGE_MEAN = 0.5;
            RANGE_STANDARD_DEVIATION = 0.5;
            
            % Transform to image instance
            imageData = utils.enums.Classes.IMAGE().convert(imageData);

            % If image is a color image, convert to gray level image
            if imageData.size.numberChannels > 1
                correctedImage = imageData.rgb2gray();
            end
            
            % Compute mean and standard deviation for image
            imageMean = mean(correctedImage.data(:));
            imageStandardDeviation = std(correctedImage.data(:));
            
            % Correct mean and standard deviation of image
            a = RANGE_STANDARD_DEVIATION / imageStandardDeviation;
            b = RANGE_MEAN - (RANGE_STANDARD_DEVIATION / imageStandardDeviation) * imageMean;
            correctedImage = image.Image( a .* imageData.data + b );
        end
    end
    
    enumeration
        GRAY_WORLD ( 'gray_world' ...
                   , @(data) image.enums.ColorCorrectionMethods.shadesOfGrayNormalization(data,1) ...
                   , 'Normalize each color channel considering the gray world algorithm.');
        MAX_RGB    ( 'max_rgb' ...
                   , @(data) image.enums.ColorCorrectionMethods.shadesOfGrayNormalization(data,Inf) ...
                   , 'Normalize each color channel considering the max-RGB algorithm.');
        SHADES_OF_GRAY ( 'gray' ...
                       , @(data) image.enums.ColorCorrectionMethods.shadesOfGrayNormalization(data) ...
                       , 'Normalize each color channel considering the shades of gray algorithm.');
        MEAN_VARIANCE  ( 'mean_var'  ...
                       , @(data) image.enums.ColorCorrectionMethods.meanVarianceNormalization(data) ...
                       , 'Normalize image by having the mean in the middle of the image range and unit variance.');
    end
end