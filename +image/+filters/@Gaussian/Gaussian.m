classdef Gaussian
    %GAUSSIAN
    %   2D Gaussian filter mask.
    
    properties (Constant)
        MATLAB_FORMAT = image.enums.Formats.MATLAB_FORMAT()
            % Matlab format encoder and decoder
    end
    
    properties
        maskSize          = image.Pixel([3,3])  % Filter mask size
        standardDeviation = 0.5                 
            % Standard deviation of the Gaussian distribution
    end
    
    properties (Dependent)
        gaussian                    % Gaussian distribution values
        derivative_u                % Gaussian u-direction derivative values
        derivative_v                % Gaussian v-direction derivative values
    end
    
    methods
        function self = Gaussian(varargin)
            %
            % Gaussian distribution instance.
            %
            % INPUTS:
            %   1. maskSize - filter mask size.
            %   2. standardDeviation - standard deviation of Gaussian
            %   distribution.
            %
            narginchk(0,2);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.maskSize = varargin{1};
                end
                
                if nargin >= 2
                    self.standardDeviation = varargin{2};
                end
            end
        end
           
        function self = set.maskSize(self,newMaskSize)
            self.maskSize = utils.enums.Classes.PIXEL().convert(newMaskSize);
        end
        
        function self = set.standardDeviation(self,newStandardDeviation)
            self.standardDeviation = newStandardDeviation;
        end
        
        function gaussian = get.gaussian(self)
            % Obtain gaussian mask. The mask is already normalized with
            % L1-norm. Switch columns to get mask in matlab format.
            matlab_gaussian = fspecial( 'gaussian' ...
                                      , [self.maskSize.v,self.maskSize.u] ...
                                      , self.standardDeviation);
            
            % Convert to internal format
            gaussian = self.MATLAB_FORMAT.decode(matlab_gaussian);
        end
        
        function derivative = get.derivative_u(self)
            % Compute the convolution not considering the first and last
            % columns of the full convolution. Convolution needs matlab
            % format.
            matlab_gaussian = self.MATLAB_FORMAT.encode(self.gaussian);
            matlab_mask     = self.MATLAB_FORMAT.encode(image.gradient.enums.DerivativeMasks.CENTRAL_DIFFERENCE().mask_u);
            derivative = conv2(matlab_gaussian,matlab_mask,'same');
            
            % Normalize the derivative mask using the L1-norm.
            derivative = derivative ./ sum(abs(derivative(:)));
            
            % Convert back to internal structure
            derivative = self.MATLAB_FORMAT.decode(derivative);
        end
        
        function derivative = get.derivative_v(self)
            % The v-direction derivative mask will be the transposed of the
            % u-direction derivative mask. The mask is normalized using the
            % L1-norm.
            
            % The transpose is only true if the mask is provided in matlab
            % format.
            derivative = self.MATLAB_FORMAT.encode(self.derivative_u)';
            derivative = self.MATLAB_FORMAT.decode(derivative);
        end
    end
end