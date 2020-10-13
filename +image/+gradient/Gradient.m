classdef Gradient
    %GRADIENT
    %   Image gradient methods and properties.
    
    properties
        du = []         % Image gradient on the u-direction (width)
        dv = []         % Image gradient on the v-direction (height)
    end
    
    properties (Dependent)
        data            % Combined information of the gradients. 
                        % The gradients are combined in the 4th dimension.
        magnitude       % Magnitude of the image gradient
        phase           % Orientation of the gradient
        size            % Size of the image gradient
    end
    
    methods
        function self = Gradient(varargin)
            %
            % Create image gradient instance.
            %
            % INPUTS:
            %   1. du - image gradient on the u-direction.
            %   2. dv - image gradient on the v-direction.
            %
            narginchk(0,2)
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.du = varargin{1};
                end
                
                if nargin >= 2
                    self.dv = varargin{2};
                end
            end
        end
        
        function self = set.du(self, new_du)
            self.du = new_du;
        end
        
        function self = set.dv(self, new_dv)
            self.dv = new_dv;
        end
        
        function data = get.data(self)
            data = cat(4,self.du,self.dv);
        end
        
        function magnitude = get.magnitude(self)
            %
            % Obtain the magnitude of the image gradient. The magnitude
            % corresponds to the norm of the image gradient on the u- and
            % v-direction.
            %
            magnitude = sqrt(sum(self.data.^2,4));
        end
        
        function phase = get.phase(self)
            %
            % Obtain the direction of the image gradient. The image
            % gradient is defined by a vector [du,dv]. So the angle is
            % defined as arctan(dv/du). The sign of the phase is obtained
            % with a minus to agree with the convention of matlab. Postive
            % gradients define negative angles.
            %
            if isempty(self.du)         % v-derivative is defined
                phase = math.Angle(-atan2(self.dv,zeros(size(self.dv))));
            elseif isempty(self.dv)     % u-derivative is defined
                phase = math.Angle(-atan2(zeros(size(self.du)),self.du));
            else                        
                phase = math.Angle(-atan2(self.dv,self.du));
            end
        end
        
        function derivativeSize = get.size(self)
            % If derivative data is empty, consider the default size values
            if isempty(self.du) && isempty(self.dv)
                derivativeSize = image.ImageSize();
            elseif isempty(self.du) 	% v-derivative is defined
                derivativeSize = image.ImageSize(size(self.dv));
            else                        % u-derivative is defined
                derivativeSize = image.ImageSize(size(self.du));
            end
        end
    end
    
    methods
        function self = interpolate(self,newSize,interpolationMethod)
            %
            % Interpolate gradient information in order to avoid aliasing.
            % This is useful when filtering results are being multiplied.
            %
            % INPUTS:
            %   1. newSize - new size to be considered for gradient. The
            %   default size is 2 times the gradient size.
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
            if nargin <= 1
                newSize = self.size;
                newSize.data = [ newSize.numberChannels ...
                               , newSize.numberPixels_u * 2 ...
                               , newSize.numberPixels_v * 2 ];
            end
            newSize = utils.enums.Classes.IMAGE_SIZE().convert(newSize);

            % If gradient data is empty, throw error
            if isempty(self.du) && isempty(self.dv)
                error = MException( 'Gradient:interpolate:noData' ...
                                  , 'Gradient data is not defined...' );
                error.throw();
            end
            
            % Only one of the gradients could be filled
            if isempty(self.du)
                dI_du = zeros(size(self.dv));
            else
                dI_du = self.du;
            end
            
            if isempty(self.dv)
                dI_dv = zeros(size(self.du));
            else
                dI_dv = self.dv;
            end
            
            % Encode to matlab format
            matlabFormat = image.enums.Formats.MATLAB_FORMAT();
            matlab_dI_du = matlabFormat.encode(dI_du);
            matlab_dI_dv = matlabFormat.encode(dI_dv);
            
            % Interpolate gradient. imresize needs matlab format to resize
            % color images. Decode matlab format to internal format
            matlab_dI_du = imresize( matlab_dI_du ...
                                   , [newSize.numberPixels_v, newSize.numberPixels_u] ...
                                   , interpolationMethod.method );
            matlab_dI_dv = imresize( matlab_dI_dv ...
                                   , [newSize.numberPixels_v, newSize.numberPixels_u] ...
                                   , interpolationMethod.method );
            self.du = matlabFormat.decode(matlab_dI_du);
            self.dv = matlabFormat.decode(matlab_dI_dv);
        end

        function show(self)
            %
            % Show image gradients side by side
            %
            
            % If there is not data, throw error
            if isempty(self.du) && isempty(self.dv)
                error = MException( 'Gradient:show:noData' ...
                                  , 'Gradient data is not defined...' );
                error.throw();
            end
            
            % Transform to internal format
            matlabFormat = image.enums.Formats.MATLAB_FORMAT();
            matlabFormat_data = matlabFormat.encode(self.data);
            displayRange = [min(matlabFormat_data(:)),max(matlabFormat_data(:))];
            
            % Display montage of gradients side by side
            montage(matlabFormat_data,'DisplayRange',displayRange,'Size',[1,2]);
            set(gca,'position',[0 0 1 1],'units','normalized')
            axis tight;
            axis off;
        end
        
        function quiver(self, imageData)
            %
            % Show gradient directions as a vector on image.
            %
            % INPUTS:
            %   1. imageData - image data.
            %
            narginchk(1,2);
            
            % If image is not provided, consider a zero valued image.
            if nargin <= 1
                imageData = zeros( self.size.numberChannels ...
                                 , self.size.numberPixels_u ...
                                 , self.size.numberPixels_v );
            end
            
            % Convert to image data
            imageInstance = utils.enums.Classes.IMAGE().convert(imageData);
            
            % Obtain gradient coordinates and remain in matlab format
            [matlab_u,matlab_v] = meshgrid( 1:self.size.numberPixels_u ...
                                          , 1:self.size.numberPixels_v );

            % Sum the components of the several channels of the gradient
            internalFormat = image.enums.Formats.INTERNAL_FORMAT();
            matlabFormat   = image.enums.Formats.MATLAB_FORMAT();
            dI_du = sum(self.du,internalFormat.channelDimension);
            dI_dv = sum(self.dv,internalFormat.channelDimension);
            
            % Plot image and gradient directions. Quiver uses internal
            % matlab format
            matlab_dI_du = matlabFormat.encode(dI_du);
            matlab_dI_dv = matlabFormat.encode(dI_dv);
            
            imageInstance.show();
            set(gca,'position',[0 0 1 1],'units','normalized')
            axis tight;
            axis off;
            hold on;
            quiver(matlab_u, matlab_v, matlab_dI_du, matlab_dI_dv);         
            hold off;
        end
        
        function write(self,outputFilepath)
            %
            % Write gradient data to files.
            %
            % INPUTS:
            %   1. outputFilepath - local filepath to write gradient data.
            %
            narginchk(1,2);
            
            % If output filepath is not provided, consider default
            % filename.
            if nargin <= 1
                outputFilepath = 'gradient_output.png';
            end
            
            % Obtain fileparts
            [filepath,filename,extension] = fileparts(outputFilepath);
            % If extension is not provided consider png extension.
            if isempty(extension)
                extension = '.png';
            end
            
            % If gradient data is empty, throw error
            if isempty(self.du) && isempty(self.dv)
                error = MException( 'Gradient:write:noData' ...
                                  , 'Gradient data is not defined...' );
                error.throw();
            end
            
            % The data for imwrite, must be formatted in matlab format.
            matlabFormat = image.enums.Formats.MATLAB_FORMAT();
            rootFilepath = [filepath, filesep, filename];
            imwrite(matlabFormat.encode(self.phase.anglesInRadians), [rootFilepath, '_phase' ,extension]);
            imwrite(matlabFormat.encode(self.magnitude), [rootFilepath,'_magnitude',extension]);
            
            if ~isempty(self.du)
                imwrite(matlabFormat.encode(self.du), [rootFilepath,'_du',extension]);
            end
            
            if ~isempty(self.dv)
                imwrite(matlabFormat.encode(self.dv), [rootFilepath,'_dv',extension]);
            end
        end
    end
    
    methods (Static)
        function self = GradientFromImage(imageInstance,gradientFilter_u)
            %
            % Create gradient instance from image data. 
            %
            % INPUTS:
            %   1. imageInstance  - image data. The image data should be
            %   provided in internal format.
            %   2. gradientFilter_u - gradient filter to be considered to
            %   compute the image gradient. The gradient filter must be
            %   provided with a mask defined in the u-direction. If you
            %   want to perform smoothing, incorporate the smoothing on the
            %   gradient filter mask. The default operator uses a Sobel
            %   mask.
            %
            narginchk(1,2)

            % If gradient filter is not provided, consider the Sobel
            % operator mask using convolution and symmetric boundary
            % conditions
            if nargin <= 1
                gradientFilter_u = image.gradient.enums.DerivativeMasks.SOBEL().mask_u;
            end
            
            % Convert input data to instances and obtain u- and v-direction
            % gradient filters
            imageInstance    = utils.enums.Classes.IMAGE().convert(imageInstance);
            gradientFilter_u = utils.enums.Classes.FILTER().convert(gradientFilter_u);
            gradientFilter_v = gradientFilter_u;
            
            % The mask in v-direction is the transpose if we have the mask
            % in matlab format.
            matlabFormat  = image.enums.Formats.MATLAB_FORMAT();
            matlab_mask_v = matlabFormat.encode(gradientFilter_v.mask)';
            gradientFilter_v.mask = matlabFormat.decode(matlab_mask_v);
            
            % Obtain image gradients in u- and v-direction
            dI_du = gradientFilter_u.apply(imageInstance);
            dI_dv = gradientFilter_v.apply(imageInstance);
            
            % Create gradient instance
            self = image.gradient.Gradient( dI_du.data, dI_dv.data );
        end
    end
end
