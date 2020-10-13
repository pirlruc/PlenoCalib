classdef (Abstract) TemplateEpipolarPlaneImage < abstract.TemplateImage
    %TEMPLATEEPIPOLARPLANEIMAGE
    %   Epipolar plane image template.
    
    properties (Abstract, Dependent)
        size            % Epipolar plane image size
    end

    methods
        function self = TemplateEpipolarPlaneImage(varargin)
            %
            % Template epipolar plane image instance.
            %
            % INPUTS:
            %   1. epiData - epipolar plane image data.
            %
            narginchk(0,1);
            
            % Create super class instance.
            self = self@abstract.TemplateImage();
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
    end
    
    methods
        function [theta_df_du,coherence] = obtainSlopeUsingStructureTensor( ...
                                              self ...
                                            , gradientFilter_u ...
                                            , tensorSmoothingFilter ...
                                            , coherenceMeasurement ...
                                            , interpolationMethod )
            %
            % Obtain the lines' slopes from analyzing the epipolar plane
            % image using the structure tensor. The structure tensor 
            % returns a coherence measurement. The slope of the lines is 
            % returned as an angle. If you want to retrieve the slope use
            % the tangent function. The angle provided is related with
            % depth. The inverse is related with disparity.
            %
            % INPUTS:
            %   1. gradientFilter_u - gradient filter to be considered to
            %   compute the image gradient. The gradient filter must be
            %   provided with a mask defined in the u-direction. If you
            %   want to perform smoothing, incorporate the smoothing on the
            %   gradient filter mask. The default operator uses a Sobel
            %   mask with a 3 x 3 gaussian mask.
            %   2. tensorSmoothingFilter - filter to be used for smoothing 
            %   the structure tensor. The default filter is a 5 x 5  
            %   gaussian mask.
            %   3. coherenceMeasurement  - coherence measurement to be
            %   considered for the structure tensor. The default is the 
            %   standard coherence measurement.
            %   4. interpolationMethod   - interpolation method. For more
            %   information see InterpolationMethods enum or see the
            %   imresize help. The default interpolation method is bicubic.
            %
            narginchk(1,5);
            
            % Define default gradient filter
            if nargin <= 1
                gaussian = image.filters.Filter(image.filters.Gaussian([3,3],0.8).gaussian);
                gradient = image.filters.Filter(image.gradient.enums.DerivativeMasks.SOBEL().mask_u);
                gradientFilter_u = image.filters.FilterCollection([gaussian,gradient]);
                gradientFilter_u = image.filters.Filter(gradientFilter_u.mask);
            end
            
            % Define default smoothing filter
            if nargin <= 2
                tensorSmoothingFilter = image.filters.Gaussian([5,5],3.2).gaussian;
            end
            
            % Define default coherence measurement
            if nargin <= 3
                coherenceMeasurement = image.gradient.enums.CoherenceMeasurements.CONTINUOUS();
            end
            
            % Define default interpolation method
            if nargin <= 4
                interpolationMethod = image.enums.InterpolationMethods.BICUBIC;
            end
            
            % If image data is empty, throw error
            if isempty(self.data)
                error = MException( 'TemplateEpipolarPlaneImage:obtainDisparityUsingStructureTensor:noData' ...
                                  , 'Epipolar plane image data is not defined...' );
                error.throw();
            end
            
            % Convert input data to instances
            gradientFilter_u      = utils.enums.Classes.FILTER().convert(gradientFilter_u);
            tensorSmoothingFilter = utils.enums.Classes.FILTER().convert(tensorSmoothingFilter);

            % Interpolate image to avoid aliasing. Namely, we are going to
            % duplicate the size of the image
            originalEpiSize = self.size;
            newEpiSize      = self.size.data * 2;
            newEpiSize(image.enums.Formats.INTERNAL_FORMAT().channelDimension) = ...
                              self.size.numberChannels;
            interpolatedEpi = self.interpolate(newEpiSize,interpolationMethod);
            
            % Smooth image and apply gradient
            gradient = image.gradient.Gradient.GradientFromImage( interpolatedEpi ...
                                                                , gradientFilter_u);
            
            % Obtain structure tensor
            structureTensor = image.gradient.StructureTensorField.StructureTensorFieldFromGradient(gradient);
            
            % Smooth and interpolate the structure tensor back to the
            % original image size
            structureTensor = structureTensor.smooth(tensorSmoothingFilter);
            structureTensor = structureTensor.interpolate( originalEpiSize ...
                                                         , interpolationMethod);
            
            % Decompose the structure tensor and obtain coherence
            % measurement.
            [coherence,phase] = structureTensor.decompose(coherenceMeasurement);
            
            % The depth is related with the slopes of the lines in the
            % epipolar plane images. The minimum eigenvalue has an
            % eigenvector that coincides with the lines in the epipolar
            % plane images. The slope is returned as an angle and the
            % tangent relates directly with depth. If you use the inverse
            % of the tangent, the value is related with disparity.
            theta_df_du = phase.minimum;
        end
    end
end
