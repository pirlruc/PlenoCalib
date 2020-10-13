classdef StructureTensorField
    %STRUCTURETENSORFIELD
    %   2D Structure tensor field methods and properties.
    properties
        data = []       % Structure tensor field data. 
            % The structure tensor has the following format: 
            %       channel x pixels_u x pixels_v x dimension x dimension
    end
    
    properties (Dependent)
        j11             % Entry (1,1) of structure tensor
        j12             % Entry (1,2) of structure tensor
        j21             % Entry (2,1) of structure tensor
        j22             % Entry (2,2) of structure tensor
        size            % Number of entries in structure tensor
    end
    
    methods
        function self = StructureTensorField(varargin)
            %
            % Create structure tensor field instance.
            %
            % INPUTS:
            %   1. structureTensorFieldData - structure tensor field data.
            %
            narginchk(0,1);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
        
        function self = set.data(self, newStructureTensorFieldData)
            self.data = newStructureTensorFieldData;
        end
        
        function entry = get.j11(self)
            % If data is not provided, return NaN
            if isempty(self.data)
                entry = nan;
            else
                entry = self.data(:,:,:,1,1);
            end
        end
        
        function entry = get.j12(self)
            % If data is not provided, return NaN
            if isempty(self.data)
                entry = nan;
            else
                entry = self.data(:,:,:,1,2);
            end
        end
        
        function entry = get.j21(self)
            % If data is not provided, return NaN
            if isempty(self.data)
                entry = nan;
            else
                entry = self.data(:,:,:,2,1);
            end
        end
        
        function entry = get.j22(self)
            % If data is not provided, return NaN
            if isempty(self.data)
                entry = nan;
            else
                entry = self.data(:,:,:,2,2);
            end
        end
        
        function tensorSize = get.size(self)
            % If data is not provided, return empty list
            if isempty(self.data)
                tensorSize = image.gradient.StructureTensorSize();
            else
                tensorSize = image.gradient.StructureTensorSize(size(self.data));
            end
        end
    end
    
    methods
        function self = interpolate(self,newSize,interpolationMethod)
            %
            % Interpolate structure tensor information.
            %
            % INPUTS:
            %   1. newSize - new size to be considered for image. The
            %   default size is 0.5 times the image size.
            %   2. interpolationMethod - interpolation method. For more
            %   information see InterpolationMethods enum or see the
            %   imresize help. The default interpolation method is bicubic.
            %
            narginchk(1,3);
            
            % Define default interpolation method
            if nargin <= 2
                interpolationMethod = image.enums.InterpolationMethods.BICUBIC;
            end
            
            % Define default size
            internalFormat = image.enums.Formats.INTERNAL_FORMAT();
            if nargin <= 1
                newSize = self.size.data * 0.5;
                newSize(internalFormat.channelDimension) = self.size.numberChannels;
                newSize = round(newSize);
            end
            newSize = utils.enums.Classes.IMAGE_SIZE().convert(newSize);

            % If structure tensor field data is empty, throw error
            if isempty(self.data)
                error = MException( 'StructureTensorField:interpolate:noData' ...
                                  , 'Structure tensor field data is not defined...' );
                error.throw();
            end
            
            % Encode to matlab format
            matlabFormat = image.enums.Formats.MATLAB_FORMAT();
            matlab_j11   = matlabFormat.encode(self.j11);
            matlab_j12   = matlabFormat.encode(self.j12);
            matlab_j22   = matlabFormat.encode(self.j22);
            
            % Interpolate image. imresize needs matlab format to resize
            % color images. 
            matlab_j11 = imresize( matlab_j11 ...
                                 , [newSize.numberPixels_v, newSize.numberPixels_u] ...
                                 , interpolationMethod.method );
            matlab_j12 = imresize( matlab_j12 ...
                                 , [newSize.numberPixels_v, newSize.numberPixels_u] ...
                                 , interpolationMethod.method );
            matlab_j22 = imresize( matlab_j22 ...
                                 , [newSize.numberPixels_v, newSize.numberPixels_u] ...
                                 , interpolationMethod.method );
                             
            % Decode matlab format to internal format. Create new structure
            % since we are changing the dimensions.
            self = image.gradient.StructureTensorField();
            self.data(:,:,:,1,1) = matlabFormat.decode(matlab_j11);
            self.data(:,:,:,1,2) = matlabFormat.decode(matlab_j12);
            self.data(:,:,:,2,2) = matlabFormat.decode(matlab_j22);
            self.data(:,:,:,2,1) = self.data(:,:,:,1,2);
        end
        
        function self = smooth(self,smoothingFilter)
            %
            % Apply smoothing filter to structure tensor field. 
            % 
            % INPUTS:
            %   1. smoothingFilter - smoothing filter to be applied to
            %   structure tensor field data. The default smoothing filter
            %   is a 5 x 5 gaussian.
            %
            narginchk(1,2);
            
            % If smoothing filter is not applied, define a gaussian filter
            if nargin <= 1
                smoothingFilter = image.filters.Gaussian([5,5],1).gaussian;
            end
            smoothingFilter = utils.enums.Classes.FILTER().convert(smoothingFilter);
            
            % If structure tensor field is empty, throw error
            if isempty(self.data)
                error = MException( 'StructureTensorField:smooth:noData' ...
                                  , 'Structure tensor field data not defined...');
                error.throw();
            end
            
            % Smooth structure tensor field
            self.data(:,:,:,1,1) = smoothingFilter.apply(self.j11).data;
            self.data(:,:,:,1,2) = smoothingFilter.apply(self.j12).data;
            self.data(:,:,:,2,2) = smoothingFilter.apply(self.j22).data;
            self.data(:,:,:,2,1) = self.data(:,:,:,1,2);
        end
        
        function [coherence,phase,eigenvectors,eigenvalues,magnitude] = ...
                            decompose(self,coherenceMeasurement)
            %
            % Perform eigenvalue decomposition for each entry of the
            % structure tensor.
            %
            % INPUTS:
            %   1. coherenceMeasurement - coherence measurement type to be
            %   used. If no coherence measurement is provided, assume the
            %   standard coherence measurement.
            %
            narginchk(1,2);
            
            if nargin <= 1
                coherenceMeasurement = image.gradient.enums.CoherenceMeasurements.CONTINUOUS();
            end
            
            % Throw error if data is not defined
            if isempty(self.data)
                error = MException( 'StructureTensorField:decompose:noData' ...
                                  , 'Structure tensor field data not defined...' );
                error.throw();
            end
            
            % Decompose structure tensor
            % Obtain orientation corresponding to the eigenvector
            % associated with the maximum eigenvalue
            phase.maximum = 0.5 .* angle(self.j11 - self.j22 + 1j .* 2 .* self.j12);
            
            % Obtain orientation corresponding to the eigenvector
            % associated with the minimum eigenvalue
            phase.minimum = phase.maximum + pi/2;
            
            % Restrict values to the interval -\pi to +\pi
            % Correct angle < -\pi
            indices = phase.minimum < -pi;
            phase.minimum(indices) = phase.minimum(indices) + 2 * pi;

            % Correct angle > +\pi
            indices = phase.minimum > +pi;
            phase.minimum(indices) = phase.minimum(indices) - 2 * pi;

            % Obtain coherence
            coherence = coherenceMeasurement.coherence(self);
            
            if nargout >= 3
                % Obtain eigenvector associated with maximum eigenvalue
                eigenvectors.maximum(:,:,:,1) =  cos(phase.maximum);
                eigenvectors.maximum(:,:,:,2) =  sin(phase.maximum);

                % Obtain eigenvector associated with minimum eigenvalue
                eigenvectors.minimum(:,:,:,1) = -sin(phase.maximum);
                eigenvectors.minimum(:,:,:,2) =  cos(phase.maximum);
            end
            
            % Obtain eigenvalues
            if nargout >= 4
                eigenvalues.maximum = 0.5 .* ( self.j11 + self.j22 ...
                                             + sqrt( (self.j11 - self.j22).^2 + 4 .* self.j12.^2 ) );
                eigenvalues.minimum = 0.5 .* ( self.j11 + self.j22 ...
                                             - sqrt( (self.j11 - self.j22).^2 + 4 .* self.j12.^2 ) );
            end
            
            % Obtain magnitude and phase
            if nargout >= 5
                magnitude = eigenvalues.maximum - eigenvalues.minimum;
            end
        end
        
        function show(self)
            %
            % Show tensor field data.
            %
            
            % If structure tensor field is empty, throw error
            if isempty(self.data)
                error = MException( 'StructureTensorField:show:noData' ...
                                  , 'Structure tensor field data not defined...');
                error.throw();
            end
            
            % Transform to internal format
            matlabFormat = image.enums.Formats.MATLAB_FORMAT();
            matlabFormat_j11  = matlabFormat.encode(self.j11);
            matlabFormat_j12  = matlabFormat.encode(self.j12);
            matlabFormat_j22  = matlabFormat.encode(self.j22);
            matlabFormat_data = cat(4,matlabFormat_j11 ,matlabFormat_j12);
            matlabFormat_data = cat(4,matlabFormat_data,matlabFormat_j12);
            matlabFormat_data = cat(4,matlabFormat_data,matlabFormat_j22);
            displayRange = [min(matlabFormat_data(:)),max(matlabFormat_data(:))];
            
            % Display montage of gradients side by side
            montage(matlabFormat_data,'DisplayRange',displayRange,'Size',[2,2]);
            set(gca,'position',[0 0 1 1],'units','normalized')
            axis tight;
            axis off;
        end
        
        function quiver(self,imageData,coherenceMeasurement)
            %
            % Show tensor field data and corresponding eigenvectors.
            %
            % INPUTS:
            %   1. imageData - image data.
            %   2. coherenceMeasurement - coherence measurement type to be
            %   used. If no coherence measurement is provided, assume the
            %   standard coherence measurement.
            %
            narginchk(1,3);
            
            % If image is not provided, consider a zero valued image.
            if nargin <= 1
                imageData = zeros( self.size.numberChannels ...
                                 , self.size.numberPixels_u ...
                                 , self.size.numberPixels_v );
            end
            
            % Use anisotropy coherence measurement
            if nargin <= 2
                coherenceMeasurement = image.gradient.enums.CoherenceMeasurements.CONTINUOUS();
            end
            
            % If structure tensor field is empty, throw error
            if isempty(self.data)
                error = MException( 'StructureTensorField:quiver:noData' ...
                                  , 'Structure tensor field data not defined...');
                error.throw();
            end
            
            % Convert to image instance
            imageInstance = utils.enums.Classes.IMAGE().convert(imageData);
            
            % Obtain eigenvector coordinates and remain in matlab format
            [matlab_u,matlab_v] = meshgrid( 1:self.size.numberPixels_u ...
                                          , 1:self.size.numberPixels_v );
            
            % Obtain eigenvectors and eigenvales
           	[coherence,~,eigenvectors] = self.decompose(coherenceMeasurement);
                        
            % Plot image and eigenvector directions. Quiver uses internal
            % matlab format
            matlabFormat = image.enums.Formats.MATLAB_FORMAT();
            matlab_minimum_x = matlabFormat.encode(coherence .* eigenvectors.minimum(:,:,:,1));
            matlab_minimum_y = matlabFormat.encode(coherence .* eigenvectors.minimum(:,:,:,2));
            matlab_maximum_x = matlabFormat.encode(coherence .* eigenvectors.maximum(:,:,:,1));
            matlab_maximum_y = matlabFormat.encode(coherence .* eigenvectors.maximum(:,:,:,2));
            
            % Display image and structure tensor
            imageInstance.show();
            set(gca,'position',[0 0 1 1],'units','normalized')
            axis tight;
            axis off;
            hold on;
            quiver(matlab_u, matlab_v, matlab_maximum_x, matlab_maximum_y,'b');         
            quiver(matlab_u, matlab_v, matlab_minimum_x, matlab_minimum_y,'r');         
            hold off;
        end
    end
    
    methods (Static)
        function self = StructureTensorFieldFromGradient(gradient)
            %
            % Create structure tensor field instance from gradient.
            % Remember that the structure tensor corresponds to a product
            % of gradients that are computed using filtering. To avoid
            % aliasing do not forget to upsample the image/gradient or 
            % low-pass the image/gradient.
            %
            % INPUTS:
            %   1. gradientInstance - gradient instance.
            %

            % Compute structure tensor. We assume that the image/gradient 
            % already interpolated.
            internalFormat = image.enums.Formats.INTERNAL_FORMAT();
            tensorField    = nan( 1, gradient.size.numberPixels_u,gradient.size.numberPixels_v ...
                                , 2,2);
            tensorField(:,:,:,1,1) = sum(gradient.du.^2, internalFormat.channelDimension);
            tensorField(:,:,:,2,2) = sum(gradient.dv.^2, internalFormat.channelDimension);
            tensorField(:,:,:,1,2) = sum(gradient.du .* gradient.dv, internalFormat.channelDimension);
            tensorField(:,:,:,2,1) = tensorField(:,:,:,1,2);

            % Sum the structure tensor channels
            tensorField =  tensorField ./ gradient.size.numberChannels;
            
            % Set to zero components that have value less that numerical
            % precision
            tensorField(abs(tensorField) < eps) = 0;
            
            % Define structure tensor field
            self = image.gradient.StructureTensorField(tensorField);
        end
    end
end