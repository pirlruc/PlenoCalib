classdef PlaneFitting
    %PLANEFITTING
    %   Utility to fit a plane using least-squares to the data provided.
    
    properties
        observations = math.Point()     % Observations for the plane in 3D space (x,y,z)
    end
    
    properties (Dependent)
        numberObservations              % Number of points observed for the plane
    end
    
    methods
        function self = PlaneFitting(varargin)
            %
            % Create plane fitting instance.
            %
            % INPUTS:
            %   1. observations - points observed for the plane in the 3D
            %   space (x,y,z). Each observation should be provided in
            %   different columns.
            %
            narginchk(0,1);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.observations = varargin{1};
                end
            end
        end
        
        function self = set.observations(self,newObservations)
            self.observations = utils.enums.Classes.POINT().convert(newObservations);
        end
        
        function number = get.numberObservations(self)
            number = self.observations.numberVectors;
        end
    end
    
    methods
        function [planeParameters,rmse,sse,residuals] = fit(self)
            %
            % Obtain parameters for plane considering an homogeneous
            % representation for the plane using a least-squares
            % methodology.
            % This methodology allows to obtain the parameters a_hat,
            % b_hat, c_hat and d_hat that describes the plane:
            %       a_hat * x + b_hat * y + c_hat * z + d_hat = 0
            %
            
            % If no data is provided, do not fit data
            if isempty(self.observations)
                error = MException( 'planeFitting:fit:noData' ...
                                  , 'No data provided...' );
                error.throw();
            end
            
            % Construct observation matrix
            observationMatrix = [ self.observations.x(:) ...
                                , self.observations.y(:) ...
                                , self.observations.z(:) ...
                                , ones(self.numberObservations,1) ];

            % The solution to this problem is given by the null space of
            % the observation matrix. The null space is obtained by
            % applying a single-value decomposition (SVD) and selecting
            % the "eigenvector" corresponding to the null "eigenvalue".
            [~,~,eigenvectors] = svd(observationMatrix);
            planeParameters    = eigenvectors(:, end);
            
            % Set to zero parameters that are smaller than precision
            planeParameters(abs(planeParameters) <= eps) = 0;
            
            % Obtain residuals and errors for least squares problem
            if nargout > 1
                residuals = abs(observationMatrix * planeParameters);
                sse  = sum(residuals.^2);               % Sum of squared errors
                rmse = sqrt(sse / length(residuals));   % Root mean squared error
            end    
            
            % Transform line parameters to structure
            planeParameters = struct( 'a', planeParameters(1) ...
                                    , 'b', planeParameters(2) ...
                                    , 'c', planeParameters(3) ...
                                    , 'd', planeParameters(4) );
        end
    end
end
