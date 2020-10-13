classdef StructureTensor
    %STRUCTURETENSOR
    %   Lightfield structure tensor methods and properties.
    %
    %   The structure tensor is a collection of 2D-EPI analysis.
    %   The structure tensor has the following format: 
    %       channel x pixel_i x pixel_j x microlens_k x microlens_l x dimension x dimension
    
    properties
        data_ik = []        % Structure tensor data from (i,k)-analysis
        data_jl = []        % Structure tensor data from (j,l)-analysis
    end
    
    properties (Dependent)
        j11_ik
        j22_ik
        j12_ik
        j21_ik
        j11_jl
        j22_jl
        j12_jl
        j21_jl
        size            % Number of entries in structure tensor
    end
    
    methods
        function self = StructureTensor(varargin)
            %
            % Lightfield structure tensor instance.
            %
            % INPUTS:
            %   1. data_ik - structure tensor data from (i,k)-analysis.
            %   2. data_jl - structure tensor data from (j,l)-analysis.
            %
            narginchk(0,2);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.data_ik = varargin{1};
                end
                
                if nargin >= 2
                    self.data_jl = varargin{2};
                end
            end
        end
        
        function self = set.data_ik(self, newStructureTensorData_ik)
            self.data_ik = newStructureTensorData_ik;
        end
        
        function self = set.data_jl(self, newStructureTensorData_jl)
            self.data_jl = newStructureTensorData_jl;
        end
        
        function tensorSize = get.size(self)
            % If data is not provided, return empty list
            if isempty(self.data_ik) && isempty(self.data_jl)
                tensorSize = lightfield.gradient.StructureTensorSize();
            elseif ~isempty(self.data_ik)
                tensorSize = lightfield.gradient.StructureTensorSize(size(self.data_ik));
            else
                tensorSize = lightfield.gradient.StructureTensorSize(size(self.data_jl));
            end
        end
        
        function entry = get.j11_ik(self)
            % If data is not provided, return NaN
            if isempty(self.data_ik)
                entry = nan;
            else
                entry = self.data_ik(:,:,:,:,:,1);
            end
        end
        
        function entry = get.j12_ik(self)
            % If data is not provided, return NaN
            if isempty(self.data_ik)
                entry = nan;
            else
                entry = self.data_ik(:,:,:,:,:,2);
            end
        end
        
        function entry = get.j21_ik(self)
            entry = self.j12_ik;
        end
        
        function entry = get.j22_ik(self)
            % If data is not provided, return NaN
            if isempty(self.data_ik)
                entry = nan;
            else
                entry = self.data_ik(:,:,:,:,:,3);
            end
        end
        
        function entry = get.j11_jl(self)
            % If data is not provided, return NaN
            if isempty(self.data_jl)
                entry = nan;
            else
                entry = self.data_jl(:,:,:,:,:,1);
            end
        end
        
        function entry = get.j12_jl(self)
            % If data is not provided, return NaN
            if isempty(self.data_jl)
                entry = nan;
            else
                entry = self.data_jl(:,:,:,:,:,2);
            end
        end
        
        function entry = get.j21_jl(self)
            entry = self.j12_jl;
        end
        
        function entry = get.j22_jl(self)
            % If data is not provided, return NaN
            if isempty(self.data_jl)
                entry = nan;
            else
                entry = self.data_jl(:,:,:,:,:,3);
            end
        end
    end
    
    methods
        function self = interpolate(self,newSize,interpolationMethod)
            %
            % Interpolate structure tensor information.
            %
            % INPUTS:
            %   1. newSize - new size to be considered for structure
            %   tensor. The default size is 0.5 times the tensor size.
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
                newSize = self.size.data * 0.5;
                newSize(lightfield.enums.Formats.INTERNAL_FORMAT().channelDimension) = ...
                          self.size.numberChannels;
                newSize.data = round(newSize.data);
            end
            newSize = utils.enums.Classes.LIGHTFIELD_SIZE().convert(newSize);

            % If structure tensor field data is empty, throw error
            if isempty(self.data_ik) && isempty(self.data_jl)
                error = MException( 'StructureTensor:interpolate:noData' ...
                                  , 'Structure tensor data is not defined...' );
                error.throw();
            end
            
            % Create lightfield instances with the structure tensor data
            tensor_ik = lightfield.Lightfield(self.data_ik);
            tensor_jl = lightfield.Lightfield(self.data_jl);
            
            % Interpolate structure tensor
            tensor_ik = tensor_ik.interpolate( newSize, interpolationMethod ...
                                             , lightfield.epi.EpipolarPlaneImage );
            tensor_jl = tensor_jl.interpolate( newSize,interpolationMethod ...
                                             , lightfield.epi.EpipolarPlaneImage );
            
            % Create new structure tensor with new size
            self = lightfield.gradient.StructureTensor();
            self.data_ik = tensor_ik.data;
            self.data_jl = tensor_jl.data;
        end
        
        function self = smooth(self,smoothingFilter)
            %
            % Apply smoothing filter to structure tensor. 
            % 
            % INPUTS:
            %   1. smoothingFilter - smoothing filter to be applied to
            %   structure tensor field data. The default smoothing filter
            %   is a 5 x 5 gaussian.
            %
            narginchk(1,2);
            
            % If smoothing filter is not applied, define a gaussian filter
            if nargin <= 1
                smoothingFilter = image.filters.Gaussian([5,5],3.2).gaussian;
            end
            smoothingFilter = utils.enums.Classes.FILTER().convert(smoothingFilter);
            
            % If structure tensor field is empty, throw error
            if isempty(self.data_ik) && isempty(self.data_jl)
                error = MException( 'StructureTensor:smooth:noData' ...
                                  , 'Structure tensor data not defined...');
                error.throw();
            end
            
            % Convert structure tensor to (i,k) epipolar plane images
            tensorEPI_ik = lightfield.epi.enums.EpipolarPlaneImageTypes.EPI_IK().encode(self.data_ik);
            tensorEPI_ik = smoothingFilter.apply(tensorEPI_ik).data;
            self.data_ik = lightfield.epi.enums.EpipolarPlaneImageTypes.EPI_IK().decode(tensorEPI_ik);
            
            % Convert structure tensor to (j,l) epipolar plane images
            tensorEPI_jl = lightfield.epi.enums.EpipolarPlaneImageTypes.EPI_JL().encode(self.data_jl);
            tensorEPI_jl = smoothingFilter.apply(tensorEPI_jl).data;
            self.data_jl = lightfield.epi.enums.EpipolarPlaneImageTypes.EPI_JL().decode(tensorEPI_jl);
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
            internalFormat = lightfield.enums.Formats.INTERNAL_FORMAT();
                            
            % Obtain (i,k)-structure tensor field
            tensorField_ik = nan( 1, gradient.size.numberPixels_i,gradient.size.numberPixels_j ...
                                   , gradient.size.numberMicrolenses_k,gradient.size.numberMicrolenses_l ...
                                , 3);
            tensorField_ik(:,:,:,:,:,1) = sum(gradient.di.^2, internalFormat.channelDimension);
            tensorField_ik(:,:,:,:,:,2) = sum(gradient.di .* gradient.dk, internalFormat.channelDimension);
            tensorField_ik(:,:,:,:,:,3) = sum(gradient.dk.^2, internalFormat.channelDimension);
            
            % Obtain (j,l)-structure tensor field
            tensorField_jl = nan( 1, gradient.size.numberPixels_i,gradient.size.numberPixels_j ...
                                   , gradient.size.numberMicrolenses_k,gradient.size.numberMicrolenses_l ...
                                , 3);
            tensorField_jl(:,:,:,:,:,1) = sum(gradient.dj.^2, internalFormat.channelDimension);
            tensorField_jl(:,:,:,:,:,2) = sum(gradient.dj .* gradient.dl, internalFormat.channelDimension);
            tensorField_jl(:,:,:,:,:,3) = sum(gradient.dl.^2, internalFormat.channelDimension);
            
            % Sum the structure tensor channels
            tensorField_ik =  tensorField_ik ./ gradient.size.numberChannels;
            tensorField_jl =  tensorField_jl ./ gradient.size.numberChannels;
            
            % Set to zero components that have value less that numerical
            % precision
            tensorField_ik(abs(tensorField_ik) < eps) = 0;
            tensorField_jl(abs(tensorField_jl) < eps) = 0;
            
            % Define structure tensor field
            self = lightfield.gradient.StructureTensor();
            self.data_ik = tensorField_ik;
            self.data_jl = tensorField_jl;
        end
    end
end

