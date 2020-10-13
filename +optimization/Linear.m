classdef Linear
    %LINEAR 
    %   Linear optimization solver. This object solves equation systems of
    %   the form A x = b or A x = 0.
    
    properties
        observationMatrix = []      % Observation matrix A
        observationVector = []      % Observation vector b
    end
    
    properties (Dependent)
        numberVariables             % Number of variables
        numberEquations             % Number of equations
        systemRank                  % Rank of the system of equations
    end
    
    methods
        function self = Linear(varargin)
            %
            % Create linear optimization.
            %
            narginchk(0,2);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.observationMatrix = varargin{1};
                end
                
                if nargin >= 2
                    self.observationVector = varargin{2};
                end
            end
        end
        
        function number = get.numberVariables(self)
            number = size(self.observationMatrix,2);
        end
        
        function number = get.numberEquations(self)
            number = size(self.observationMatrix,1);
        end
        
        function value = get.systemRank(self)
            value = rank(self.observationMatrix);
        end
    end
    
    methods
        function [solution,errorFlag,rmse,residuals] = solve(self)
            %
            % Solve system of equations of the form A x = b or A x = 0.
            %
            % The solution is good if errorFlag = 0. If the rank of the
            % matrix is not sufficient to get a solution, the errorFlag =
            % -1.

            % Validate if matrix has sufficient rank to obtain a solution
            if self.systemRank < size(self.observationMatrix,2) - 1
                errorFlag = -1;
            else
                errorFlag = 0;
            end

            % Determine which system of equations should be solved
            if sum(self.observationVector(:)) == 0
                if nargout > 2
                    [solution,rmse,residuals] = self.solveAx_0();
                else
                    solution = self.solveAx_0();
                end
            else
                if nargout > 2
                    [solution,rmse,residuals] = self.solveAx_b();
                else
                    solution = self.solveAx_b();
                end
            end
        end
        
        function [solution,rmse,residuals] = solveAx_0(self)
            %
            % Solve system of equations of the form A x = 0.
            %
            % This solution is based on the code of Qi Zhang.
            %
            
            % Obtain solution using singular value decomposition
            [~,singularValues,singularVectors] = svd(self.observationMatrix);
            singularValues = diag(singularValues);
            
            % If number of equations is greater than the number of 
            % variables, obtain the solution with the lowest singular value
            if self.numberEquations > self.numberVariables 
                [~,minimumIndex] = min(abs(singularValues));
                
                % If singular value is negative, the solution is the
                % negative of the one computed.
                if singularValues(minimumIndex) >= 0
                    solutionSVD = singularVectors(:,minimumIndex);
                else
                    solutionSVD = -singularVectors(:,minimumIndex);
                end
            else
                solutionSVD = singularVectors(:,end);
            end
            
            % Obtain error associated with the singular value decomposition
            errorSVD = norm(self.observationMatrix * solutionSVD);

            % Obtain solution using eigenvalue decomposition
            [eigenvectors,eigenvalues] = eig(self.observationMatrix' * self.observationMatrix);
            eigenvalues      = diag(eigenvalues);
            [~,minimumIndex] = min(abs(eigenvalues));

            % If eigenvalue is negative, the solution is the negative of 
            % the one computed.
            if eigenvalues(minimumIndex) >= 0
                solutionEIG = eigenvectors(:,minimumIndex);
            else
                solutionEIG = -eigenvectors(:,minimumIndex);
            end
            
            errorEIG = norm(self.observationMatrix * solutionEIG);
      
            % Obtain the solution with the smaller error
            if errorSVD < errorEIG
                solution = solutionSVD;
            else
                solution = solutionEIG;
            end
            
            % Obtain errors
            if nargout > 1
                residuals = abs(self.observationMatrix * solution);
                sse  = sum(residuals.^2);               % Sum of squared errors
                rmse = sqrt(sse / length(residuals));   % Root mean squared error
            end
        end
        
        function [solution,rmse,residuals] = solveAx_b(self)
            %
            % Solve system of equations of the form A x = b.
            %

            % Obtain solution
            solution = mldivide(self.observationMatrix,self.observationVector);
            
            % Obtain errors
            if nargout > 1
                residuals = abs(self.observationMatrix * solution - self.observationVector);
                sse  = sum(residuals.^2);               % Sum of squared errors
                rmse = sqrt(sse / length(residuals));   % Root mean squared error
            end
        end
    end
end
