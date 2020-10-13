classdef PolynomialFitting
    %POLYNOMIALFITTING
    %   Utility to fit a polynomial using least-squares to the data 
    %   provided.
    
    properties
        observationInputs  = []     % Observed inputs for the polynomial  (x)
        observationOutputs = []     % Observed outputs for the polynomial (y)
        polynomialDegree   = 0      % Polynomial degree
    end
    
    properties (Dependent)
        numberObservations          % Number of points observed for the polynomial
    end
    
    methods
        function self = PolynomialFitting(varargin)
            %
            % Create polynomial fitting instance.
            %
            % INPUTS:
            %   1. observationInputs  - input points observed for the
            %   polynomial.
            %   2. observationOutputs - output points observed for the 
            %   polynomial.
            %   3. polynomialDegree   - polynomial degree.
            %
            narginchk(0,3);
            
            if ~isempty(varargin)
                if nargin >= 3
                    self.polynomialDegree = varargin{3};
                end
                
                if nargin >= 2
                    self.observationOutputs = varargin{2};
                end
                
                if nargin >= 1
                    self.observationInputs = varargin{1};
                end
            end
        end
        
        function self = set.polynomialDegree(self,newPolynomialDegree)
            self.polynomialDegree = newPolynomialDegree;
        end
        
        function self = set.observationInputs(self,newObservationInputs)
            self.observationInputs = newObservationInputs;
        end
        
        function self = set.observationOutputs(self,newObservationOutputs)
            self.observationOutputs = newObservationOutputs;
        end
        
        function number = get.numberObservations(self)
            number = length(self.observationInputs);
        end
    end
    
    methods
        function [polynomialParameters,rmse,sse,residuals] = fit(self)
            %
            % Obtain parameters for polynomial considering a least-squares
            % approximation.
            % This methodology allows to obtain the parameters that 
            % describe the polynomial:
            %   aN_hat * x^N + ... + a2_hat * x^2 + a1_hat * x + a0_hat = y
            %
            
            % If no data is provided, do not fit data
            if isempty(self.observationInputs) || isempty(self.observationOutputs)
                error = MException( 'PolynomialFitting:fit:noData' ...
                                  , 'No input or output data provided...' );
                error.throw();
            end
            
            % Construct observation matrix and observation output
            observationOutput = self.observationOutputs(:);
            observationMatrix = [];
            for degree = self.polynomialDegree:-1:0
                observationMatrix = cat(2, observationMatrix, self.observationInputs(:) .^ degree);
            end
            
            % This problem has the form of Ax = b and therefore can be
            % solved using x = ((A'A)^-1)A'b.
            parameters = mldivide(observationMatrix,observationOutput);
            
            % Set to zero parameters that are smaller than precision
            parameters(abs(parameters) <= eps) = 0;
            
            % Obtain residuals and errors for least squares problem
            if nargout > 1
                residuals = abs(observationMatrix * parameters - observationOutput);
                sse  = sum(residuals.^2);               % Sum of squared errors
                rmse = sqrt(sse / length(residuals));   % Root mean squared error
            end    
            
            % Transform polynomial parameters to structure
            for iDegree = 1:self.polynomialDegree + 1
                polynomialParameters.(['a' num2str(self.polynomialDegree + 1 - iDegree)]) = parameters(iDegree);
            end
        end
    end
end
