classdef Pixel < math.Vector
    %Pixel
    %   Utility to represent pixel data.
    
    properties (Constant)
        NUMBER_COMPONENTS = 2       % Number of components in pixel
    end
    
    properties (Dependent)
        u       % u-component for pixel
        v       % v-component for pixel
    end
    
    methods
        function self = Pixel(varargin)
            %
            % Pixel instance.
            %
            % INPUTS:
            %   1. pixelData - pixel entries. To represent several pixels,
            %      use a matrix with each pixel being represented as a new
            %      column.
            %   2. homogeneous - homogeneous flag to indicate if pixel is
            %      in homogeneous coordinates or not.
            %
            narginchk(0,2);
            
            % Create super class instance
            self = self@math.Vector();
            
            if ~isempty(varargin)
                if nargin >= 2
                    self.homogeneous = varargin{2};
                end
                
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
        
        function data = get.u(self)
            % If data is empty, return empty list
            if isempty(self.data)
                data = [];
            else
                data = self.data(1,:);
            end
        end
        
        function data = get.v(self)
            % If data is empty, return empty list
            if isempty(self.data)
                data = [];
            else
                data = self.data(2,:);
            end
        end
    end
    
    methods (Access = protected)
        function data = correctInputData(self,newData)
            %
            % Correct input data to have points defined in each column.
            %
            % INPUT:
            %   1. newData - input data to point.
            %
            
            % Different points should be provided in different columns
            newData = math.Vector(newData);
            if     self.homogeneous == false ...
                && newData.numberComponents ~= self.NUMBER_COMPONENTS
                data = newData.data';
            elseif self.homogeneous == true ...
                && newData.numberComponents ~= (self.NUMBER_COMPONENTS + 1)
                data = newData.data';
            else
                data = newData.data;
            end
        end
    end
    
    methods
        function pixels = obtainVectors(self,pixelPositions)
            %
            % Obtain pixels from the list of pixels in current object.
            %
            % INPUTS:
            %   1. pixelPositions - positions of the pixels in the list of
            %      pixels.
            %
            
            % Obtain vector
            vectors = obtainVectors@math.Vector(self,pixelPositions);
            
            % Convert to pixel structure
            pixels      = image.Pixel();
            pixels.homogeneous = vectors.homogeneous;
            pixels.data        = vectors.data;
        end
        
        function self = normalize(self)
            %
            % Normalize pixels in order to have the pixels centered at the
            % origin with an average distance of sqrt(2).
            %
            
            % Set homogeneous coordinates
            self = self.setHomogeneousCoordinates();

            % Normalize pixels
            self.data = self.obtainNormalizationMatrix * self.data;
            
            % Remove homogeneous coordinates
            self = self.removeHomogeneousCoordinates();
        end
        
        function normalizationMatrix = obtainNormalizationMatrix(self)
            %
            % Obtain normalization matrix to have the pixels centered at 
            % the origin with an average distance of sqrt(2).
            %
            % Obtain mean for pixels (u,v) 
            shift_u = mean(self.u);
            shift_v = mean(self.v);

            % Scale pixels to have root mean squared distance equal to
            % sqrt(2).
            scalePixels = sqrt(2) / sqrt(mean( (self.u - shift_u).^2 ...
                                             + (self.v - shift_v).^2 ));

            % Obtain normalization matrix
            normalizationMatrix  = [ scalePixels 0 -scalePixels * shift_u ...
                                   ; 0 scalePixels -scalePixels * shift_v ...
                                   ; 0 0 1 ];
        end
        
        function line = pixels2line(self)
            % 
            % Convert pixels to line parameters.
            %
            
            % Obtain line parameters
            lineFitting = math.LineFitting();
            lineFitting.observationInputs  = self.u;
            lineFitting.observationOutputs = self.v;

            % Obtain line parameters
            lineParameters = lineFitting.fit();
            lineParameters = [lineParameters.a;lineParameters.m;lineParameters.b];

            % Convert to vector class
            line = math.Vector(lineParameters,false);
        end
    end
end
