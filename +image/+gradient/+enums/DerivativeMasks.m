classdef DerivativeMasks
    %DERIVATIVEMASKS
    %   Masks to be used to compute image gradients. The masks are provided
    %   in internal format.
    
    properties (~Dependent, SetAccess = immutable, GetAccess = public)
        mask_u          % 2D masks to obtain image derivatives in the u-direction
        mask_v          % 2D masks to obtain image derivatives in the v-direction
        description     % Description of the masks
    end
    
    methods
        function self = DerivativeMasks(mask, description)
            %
            % Derivative masks enumeration instance. Notice that the masks
            % are normalized using L1-norm.
            %
            self.mask_u      = mask ./ sum(abs(mask(:)));
            self.mask_v      = self.mask_u';
            self.description = description;
            
            % Convert masks to internal format
            self.mask_u = image.enums.Formats.MATLAB_FORMAT().decode(self.mask_u);
            self.mask_v = image.enums.Formats.MATLAB_FORMAT().decode(self.mask_v);
        end
    end
    
    enumeration
        SOBEL   ([+1 0 -1; +2  0 -2 ; +1 0 -1], 'Sobel edge mask')
        PREWITT ([+1 0 -1; +1  0 -1 ; +1 0 -1], 'Prewitt edge mask')
        SCHARR  ([+3 0 -3; +10 0 -10; +3 0 -3], 'Scharr edge mask')
        CENTRAL_DIFFERENCE ([+1 0 -1], 'Central difference mask')
        FORWARD_DIFFERENCE ([+1 -1], 'Forward difference mask')
        BACKWARD_DIFFERENCE([-1 +1], 'Backward difference mask')
        LAPLACIAN_4N([ 0 +1  0; +1 -4 +1;  0 +1  0],'Laplacian edge mask with 4-neighborhoods')
        LAPLACIAN_8N([+1 +1 +1; +1 -8 +1; +1 +1 +1],'Laplacian edge mask with 8-neighborhoods')
    end
end