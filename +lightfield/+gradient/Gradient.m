classdef Gradient
    %GRADIENT
    %   Lightfield gradient methods and properties.
    
    properties
        di = []         % Lightfield gradient on the i-direction
        dj = []         % Lightfield gradient on the j-direction
        dk = []         % Lightfield gradient on the k-direction
        dl = []         % Lightfield gradient on the v-direction
    end
    
    properties (Dependent)
        data            % Combined information of the gradients. 
                        % The gradients are combined in the 6th dimension.
        magnitude       % Magnitude of the lightfield gradient
        magnitude_ik    % Magnitude of the lightfield gradient in (i,k)-direction
        magnitude_jl    % Magnitude of the lightfield gradient in (j,l)-direction
        size            % Size of the lightfield gradient
    end
    
    methods
        function self = Gradient(varargin)
            %
            % Create lightfield gradient instance.
            %
            % INPUTS:
            %   1. di - lightfield gradient on the i-direction.
            %   2. dj - lightfield gradient on the j-direction.
            %   3. dk - lightfield gradient on the k-direction.
            %   4. dl - lightfield gradient on the l-direction.
            %
            narginchk(0,4)
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.di = varargin{1};
                end
                
                if nargin >= 2
                    self.dj = varargin{2};
                end
                
                if nargin >= 3
                    self.dk = varargin{3};
                end
                
                if nargin >= 4
                    self.dl = varargin{4};
                end
            end
        end
        
        function self = set.di(self, new_di)
            self.di = new_di;
        end
        
        function self = set.dj(self, new_dj)
            self.dj = new_dj;
        end
        
        function self = set.dk(self, new_dk)
            self.dk = new_dk;
        end
        
        function self = set.dl(self, new_dl)
            self.dl = new_dl;
        end
        
        function data = get.data(self)
            if ~isempty(self.di)
                data = self.di;
                data = cat(6,data,self.dj);
                data = cat(6,data,self.dk);
                data = cat(6,data,self.dl);
            elseif ~isempty(self.dj)
                data = self.dj;
                data = cat(6,data,self.dk);
                data = cat(6,data,self.dl);
            elseif ~isempty(self.dk)
                data = self.dk;
                data = cat(6,data,self.dl);
            else
                data = self.dl;
            end
        end
        
        function magnitude = get.magnitude(self)
            %
            % Obtain the magnitude of the lightfield gradient. The 
            % magnitude corresponds to the norm of the lightfield gradient 
            % on the i-, j-, k- and l-direction.
            %
            magnitude = sqrt(sum(self.data.^2,6));
        end
        
        function magnitude = get.magnitude_ik(self)
            %
            % Obtain the magnitude of the lightfield gradient considering
            % the epipolar plane image (i,k). The magnitude corresponds to 
            % the norm of the lightfield gradients in directions (i,k).
            %
            % Concatenate gradients
            gradient_ik = cat(6,self.di,self.dk);
            
            % Obtain magnitudes for gradients on epipolar plane images
            % using the l2-norm for the gradient vector
            magnitude = sqrt(sum(gradient_ik.^2,6));
        end
        
        function magnitude = get.magnitude_jl(self)
            %
            % Obtain the magnitude of the lightfield gradient considering
            % the epipolar plane image (j,l). The magnitude corresponds to 
            % the norm of the lightfield gradients in directions (j,l).
            %
            % Concatenate gradients
            gradient_jl = cat(6,self.dj,self.dl);
            
            % Obtain magnitudes for gradients on epipolar plane images
            % using the l2-norm for the gradient vector
            magnitude = sqrt(sum(gradient_jl.^2,6));
        end
        
        function derivativeSize = get.size(self)
            % If derivative data is empty, consider the default size values
            if ~isempty(self.di)        % i-derivative is defined
                derivativeSize = lightfield.LightfieldSize(size(self.di));
            elseif ~isempty(self.dj) 	% j-derivative is defined
                derivativeSize = lightfield.LightfieldSize(size(self.dj));
            elseif ~isempty(self.dk) 	% k-derivative is defined
                derivativeSize = lightfield.LightfieldSize(size(self.dk));
            elseif ~isempty(self.dl) 	% l-derivative is defined
                derivativeSize = lightfield.LightfieldSize(size(self.dl));
            else                        % No derivative is defined
                derivativeSize = lightfield.LightfieldSize();
            end
        end
    end
    
    methods (Static)
        function self = GradientFromLightfieldUsingEpipolarPlaneImages( ...
                                                lightfieldInstance ...
                                              , gradientFilter_k )
            %
            % Create gradient instance from lightfield data by computing 2D
            % gradients on epipolar plane images.
            %
            % INPUTS:
            %   1. lightfieldInstance  - lightfield data. The lightfield 
            %   data should be provided in internal format.
            %   2. gradientFilter_k - gradient filter to be considered to
            %   compute the epipolar plane image gradient. The gradient 
            %   filter must be provided with a mask defined in the 
            %   k-direction. The default operator uses a Sobel mask.
            %
            narginchk(1,2)

            % If gradient filter is not provided, consider the Sobel
            % operator mask using convolution and symmetric boundary
            % conditions
            if nargin <= 1
                gradientFilter_k = image.gradient.enums.DerivativeMasks.SOBEL().mask_u;
            end
            
            % Convert input data to instances and obtain gradient filters
            % for each direction
            lightfieldInstance = utils.enums.Classes.LIGHTFIELD().convert(lightfieldInstance);
            gradientFilter_k = utils.enums.Classes.FILTER().convert(gradientFilter_k);
            gradientFilter_i = gradientFilter_k;
            
            % The mask in i-direction is the transpose if we have the mask
            % in matlab format.
            matlabFormat  = image.enums.Formats.MATLAB_FORMAT();
            matlab_mask_i = matlabFormat.encode(gradientFilter_i.mask)';
            gradientFilter_i.mask = matlabFormat.decode(matlab_mask_i);
            
            % Obtain epipolar plane image formats from lightfield
            epipolarPlaneImages_ik = lightfield.epi.enums.EpipolarPlaneImageTypes.EPI_IK().encode(lightfieldInstance.data);
            epipolarPlaneImages_jl = lightfield.epi.enums.EpipolarPlaneImageTypes.EPI_JL().encode(lightfieldInstance.data);
                                 
            % Obtain epipolar plane image gradients 
            dL_di = gradientFilter_i.apply(epipolarPlaneImages_ik);
            dL_dj = gradientFilter_i.apply(epipolarPlaneImages_jl);
            dL_dk = gradientFilter_k.apply(epipolarPlaneImages_ik);
            dL_dl = gradientFilter_k.apply(epipolarPlaneImages_jl);
            
            % Transform gradients to lightfield format
            dL_di = lightfield.epi.enums.EpipolarPlaneImageTypes.EPI_IK().decode(dL_di.data);
            dL_dj = lightfield.epi.enums.EpipolarPlaneImageTypes.EPI_JL().decode(dL_dj.data);
            dL_dk = lightfield.epi.enums.EpipolarPlaneImageTypes.EPI_IK().decode(dL_dk.data);
            dL_dl = lightfield.epi.enums.EpipolarPlaneImageTypes.EPI_JL().decode(dL_dl.data);
            
            % Create gradient instance
            self = lightfield.gradient.Gradient(dL_di,dL_dj,dL_dk,dL_dl);
        end
        
        function self = GradientFromLightfieldUsingViewpointAndMicrolensImages( ...
                                                lightfieldInstance ...
                                              , viewpointGradient_k ...
                                              , microlensGradient_i )
            %
            % Create gradient instance from lightfield data by computing 2D
            % gradients on viewpoint and microlens images.
            %
            % INPUTS:
            %   1. lightfieldInstance  - lightfield data. The lightfield 
            %   data should be provided in internal format.
            %   2. viewpointGradient_k - gradient filter to be considered 
            %   compute the viewpoint image gradient. The gradient filter 
            %   must be provided with a mask defined in the k-direction.
            %   The default operator uses a Sobel mask.
            %   3. microlensGradient_i - gradient filter to be considered 
            %   compute the microlens image gradient. The gradient filter 
            %   must be provided with a mask defined in the i-direction.
            %   The default operator uses a Sobel mask.
            %
            narginchk(1,3)

            % If gradient filter is not provided, consider the Sobel
            % operator mask using convolution and symmetric boundary
            % conditions
            if nargin <= 1
                viewpointGradient_k = image.gradient.enums.DerivativeMasks.SOBEL().mask_u;
            end
            
            if nargin <= 2
                microlensGradient_i = image.gradient.enums.DerivativeMasks.SOBEL().mask_u;
            end
            
            % Obtain gradient from microlens images
            microlensGradient = lightfield.gradient.Gradient.GradientFromLightfieldUsingMicrolensImages( lightfieldInstance ...
                                                                                                       , microlensGradient_i );

            % Obtain gradient from viewpoint images
            viewpointGradient = lightfield.gradient.Gradient.GradientFromLightfieldUsingViewpointImages( lightfieldInstance ...
                                                                                                       , viewpointGradient_k );
            
            % Create gradient instance
            self = lightfield.gradient.Gradient( microlensGradient.di,microlensGradient.dj ...
                                               , viewpointGradient.dk,viewpointGradient.dl );
        end
        
        function self = GradientFromLightfieldUsingMicrolensImages( lightfieldInstance ...
                                                                  , microlensGradient_i )
            %
            % Create gradient instance from lightfield data by computing 2D
            % gradients on microlens images.
            %
            % INPUTS:
            %   1. lightfieldInstance  - lightfield data. The lightfield 
            %   data should be provided in internal format.
            %   2. microlensGradient_i - gradient filter to be considered 
            %   compute the microlens image gradient. The gradient filter 
            %   must be provided with a mask defined in the i-direction.
            %   The default operator uses a Sobel mask.
            %
            narginchk(1,2)

            % If gradient filter is not provided, consider the Sobel
            % operator mask using convolution and symmetric boundary
            % conditions
            if nargin <= 1
                microlensGradient_i = image.gradient.enums.DerivativeMasks.SOBEL().mask_u;
            end
            
            % Convert input data to instances and obtain gradient filters
            % for each direction
            lightfieldInstance = utils.enums.Classes.LIGHTFIELD().convert(lightfieldInstance);
            microlensGradient_i = utils.enums.Classes.FILTER().convert(microlensGradient_i);
            microlensGradient_j = microlensGradient_i;
            
            % The mask in j-direction is the transpose if we have the mask
            % in matlab format.
            matlabFormat  = image.enums.Formats.MATLAB_FORMAT();
            matlab_mask_j = matlabFormat.encode(microlensGradient_j.mask)';
            microlensGradient_j.mask = matlabFormat.decode(matlab_mask_j);
            
            % Obtain microlens image formats from lightfield
            microlensImages = permute( lightfieldInstance.data ...
                                     , lightfield.image.MicrolensImage.ENCODING_ORDER );
                                 
            % Obtain viewpoint and microlens image gradients 
            dL_di = microlensGradient_i.apply(microlensImages);
            dL_dj = microlensGradient_j.apply(microlensImages);
            
            % Transform gradients to lightfield format
            dL_di = permute(dL_di.data,lightfield.image.MicrolensImage.DECODING_ORDER);
            dL_dj = permute(dL_dj.data,lightfield.image.MicrolensImage.DECODING_ORDER);
            
            % Create gradient instance
            self = lightfield.gradient.Gradient(dL_di,dL_dj,[],[]);
        end
        
        function self = GradientFromLightfieldUsingViewpointImages( lightfieldInstance ...
                                                                  , viewpointGradient_k )
            %
            % Create gradient instance from lightfield data by computing 2D
            % gradients on viewpoint images.
            %
            % INPUTS:
            %   1. lightfieldInstance  - lightfield data. The lightfield 
            %   data should be provided in internal format.
            %   2. viewpointGradient_k - gradient filter to be considered 
            %   compute the viewpoint image gradient. The gradient filter 
            %   must be provided with a mask defined in the k-direction.
            %   The default operator uses a Sobel mask.
            %
            narginchk(1,2)

            % If gradient filter is not provided, consider the Sobel
            % operator mask using convolution and symmetric boundary
            % conditions
            if nargin <= 1
                viewpointGradient_k = image.gradient.enums.DerivativeMasks.SOBEL().mask_u;
            end
            
            % Convert input data to instances and obtain gradient filters
            % for each direction
            lightfieldInstance = utils.enums.Classes.LIGHTFIELD().convert(lightfieldInstance);
            viewpointGradient_k = utils.enums.Classes.FILTER().convert(viewpointGradient_k);
            viewpointGradient_l = viewpointGradient_k;
            
            % The mask in l--direction is the transpose if we have the mask
            % in matlab format.
            matlabFormat  = image.enums.Formats.MATLAB_FORMAT();
            matlab_mask_l = matlabFormat.encode(viewpointGradient_l.mask)';
            viewpointGradient_l.mask = matlabFormat.decode(matlab_mask_l);
            
            % Obtain viewpoint image formats from lightfield
            viewpointImages = permute( lightfieldInstance.data ...
                                     , lightfield.image.ViewpointImage.ENCODING_ORDER );
                                 
            % Obtain viewpoint and microlens image gradients 
            dL_dk = viewpointGradient_k.apply(viewpointImages);
            dL_dl = viewpointGradient_l.apply(viewpointImages);
            
            % Transform gradients to lightfield format
            dL_dk = permute(dL_dk.data,lightfield.image.ViewpointImage.DECODING_ORDER);
            dL_dl = permute(dL_dl.data,lightfield.image.ViewpointImage.DECODING_ORDER);
            
            % Create gradient instance
            self = lightfield.gradient.Gradient([],[],dL_dk,dL_dl);
        end
        
        function self = GradientFromLightfield(lightfieldInstance)
            %
            % Create gradient instance from lightfield data. The gradient
            % computes the central difference in each combination of
            % dimensions of the lightfield.
            %
            % INPUTS:
            %   1. lightfieldInstance  - lightfield data. The lightfield 
            %   data should be provided in internal format.
            %
            narginchk(1,1)

            % Convert input data to instances and obtain u- and v-direction
            % gradient filters
            lightfieldInstance = utils.enums.Classes.LIGHTFIELD().convert(lightfieldInstance);

            % Obtain central difference for each dimension of the
            % lightfield
            [~,di,dj,dk,dl] = gradient(lightfieldInstance.data);
            
            % Create gradient instance
            self = lightfield.gradient.Gradient(di,dj,dk,dl);
        end
    end
end
