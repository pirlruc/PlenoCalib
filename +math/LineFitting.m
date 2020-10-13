classdef LineFitting
    %LINEFITTING
    %   Utility to fit a line using least-squares to the data provided.
    
    properties
        observationInputs  = []     % Observed inputs for the line  (x)
        observationOutputs = []     % Observed outputs for the line (y)
    end
    
    properties (Dependent)
        numberObservations          % Number of points observed for the line
    end
    
    methods
        function self = LineFitting(varargin)
            %
            % Create line fitting instance.
            %
            % INPUTS:
            %   1. observationInputs  - input points observed for the line.
            %   2. observationOutputs - output points observed for the 
            %      line.
            %
            narginchk(0,2);
            
            if ~isempty(varargin)
                if nargin >= 2
                    self.observationOutputs = varargin{2};
                end
                
                if nargin >= 1
                    self.observationInputs = varargin{1};
                end
            end
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
        function [lineParameters,rmse,sse,residuals] = fit(self)
            %
            % Obtain parameters for line considering an homogeneous
            % representation for the line using a least-squares
            % methodology.
            % This methodology allows to obtain the parameters a_hat, m_hat
            % and b_hat that describes the line:
            %       a_hat * y + m_hat * x + b_hat = 0
            %
            
            % If no data is provided, do not fit data
            if isempty(self.observationInputs) || isempty(self.observationOutputs)
                error = MException( 'lineFitting:fit:noData' ...
                                  , 'No input or output data provided...' );
                error.throw();
            end
            
            % If only 2 observation are given use cross product
            if self.numberObservations == 2
                point_A        = [self.observationInputs(1),self.observationOutputs(1),1];
                point_B        = [self.observationInputs(2),self.observationOutputs(2),1];
                lineParameters = cross(point_A,point_B);
                lineParameters = lineParameters ./ sqrt(sum(lineParameters.^2));
                
                if nargout > 1
                    residuals = 0;
                    sse  = 0;
                    rmse = 0;
                end
                
            % If more than 2 pixels are given use line fitting
            else
                % Construct observation matrix
                observationMatrix = [ self.observationOutputs(:) ...
                                    , self.observationInputs(:) ...
                                    , ones(self.numberObservations,1) ];

                % The solution to this problem is given by the null space of
                % the observation matrix. The null space is obtained by
                % applying a single-value decomposition (SVD) and selecting
                % the "eigenvector" corresponding to the null "eigenvalue".
                [~,~,eigenvectors] = svd(observationMatrix);
                lineParameters     = eigenvectors(:, end);

                % Set to zero parameters that are smaller than precision
                lineParameters(abs(lineParameters) <= eps) = 0;

                % Obtain residuals and errors for least squares problem
                if nargout > 1
                    residuals = abs(observationMatrix * lineParameters);
                    sse  = sum(residuals.^2);               % Sum of squared errors
                    rmse = sqrt(sse / length(residuals));   % Root mean squared error
                end    
            end
            % Transform line parameters to structure
            lineParameters = struct( 'a', lineParameters(1) ...
                                   , 'm', lineParameters(2) ...
                                   , 'b', lineParameters(3) );
        end
    end
    
    methods (Static)
        function newLineParameters = convertLineParameters_y(lineParameters)
            %
            % Convert line parameters to have:
            %       y = -(1/a_hat) * (m_hat * x + b_hat).
            %
            newLineParameters = struct( 'a', 1 ...
                                      , 'm', - lineParameters.m /  lineParameters.a ...
                                      , 'b', - lineParameters.b /  lineParameters.a );
        end
        
        function newLineParameters = convertLineParameters_x(lineParameters)
            %
            % Convert line parameters to have:
            %       x = -(1/m_hat) * (a_hat * x + b_hat).
            %
            newLineParameters = struct( 'a', - lineParameters.a /  lineParameters.m ...
                                      , 'm', 1 ...
                                      , 'b', - lineParameters.b /  lineParameters.m );
        end
    end
end
