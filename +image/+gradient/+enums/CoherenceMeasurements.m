classdef CoherenceMeasurements
    %COHERENCEMEASUREMENTS
    %   Coherence measurements available for structure tensor.
    
    properties (~Dependent, SetAccess = immutable, GetAccess = public)
        symbol          % Coherence measurement
        description     % Description of the masks
        method          % Method to compute the coherence measurement
    end
    
    methods
        function self = CoherenceMeasurements(symbol, description, method)
            %
            % Coherence measurements enumeration instance.
            %
            self.symbol = symbol;
            self.method = method;
            self.description = description;
        end
    end
    
    methods 
        function value = coherence(self,data)
            %
            % Obtain coherence measurement.
            %
            % INPUTS:
            %   1. data - data to compute coherence measurement.
            %
            narginchk(1,2);
            
            % Convert to structure tensor field
            if ~ismatrix(data) == 2
                tensor = utils.enums.Classes.STRUCTURE_TENSOR_FIELD().convert(data);
            
            % Convert to structure tensor if it is not a structure tensor
            % field
            elseif ~isa(data,utils.enums.Classes.STRUCTURE_TENSOR_FIELD().class)
                tensor = utils.enums.Classes.STRUCTURE_TENSOR().convert(data);
            
            % Data is already a tensor
            else 
                tensor = data;
            end
            
            % If data is empty, return NaN
            if isempty(tensor.data)
                value = nan;
            else
                % Obtain coherence
                value = self.method(tensor); 
            end
        end
    end
    
    methods (Static)
        function value = obtainAnisotropy(tensor)
            %
            % Anisotropy coherence measurement.
            %
            % INPUTS:
            %   1. data - data to compute coherence measurement.
            %
            narginchk(1,1);
            
            % Obtain magnitude and eigenvalues
            magnitude         = sqrt((tensor.j22 - tensor.j11).^2 + 4 .* tensor.j12.^2);
            minimumEigenvalue = 0.5 .* (tensor.j11 + tensor.j22 - magnitude);
            maximumEigenvalue = 0.5 .* (tensor.j11 + tensor.j22 + magnitude);
            
            % Obtain coherence measurement
            value = ( magnitude  ...
                  ./ (minimumEigenvalue + maximumEigenvalue) ).^2;

            % The coherence should be zero for all indices were the sum of 
            % the eigenvalues is zero.
            indices = minimumEigenvalue + maximumEigenvalue <= eps;
            value(indices) = 0;
        end
        
        function value = obtainContinuous(tensor)
            %
            % Standard continuous coherence measurement.
            %
            % INPUTS:
            %   1. data - data to compute coherence measurement.
            %
            narginchk(1,1);
            
            % Obtain continuous coherence measurement
            value = sqrt((tensor.j22 - tensor.j11).^2 + 4 .* tensor.j12.^2);
        end

        function value = obtainForstner(tensor)
            %
            % Forstner coherence measurement.
            %
            % INPUTS:
            %   1. data - data to compute coherence measurement.
            %
            narginchk(1,1);
            
            % Obtain Forstner measurement
            value  =  (tensor.j11 .* tensor.j22 - tensor.j12.^2) ...
                   ./ (tensor.j11 + tensor.j22);
        end

        function value = obtainCornerResponseFunction(tensor)
            %
            % Corner response function coherence measurement.
            %
            % INPUTS:
            %   1. data - data to compute coherence measurement.
            %
            narginchk(1,1);
            
            % Obtain corner response function
            value  = (tensor.j11 .* tensor.j22 - tensor.j12.^2) ...
                   - 0.04 .* (tensor.j11 + tensor.j22).^2;
        end
        
        function value = obtainConditionNumber(tensor)
            %
            % Condition number coherence measurement.
            %
            % INPUTS:
            %   1. data - data to compute coherence measurement.
            %
            narginchk(1,1);
            
            % Obtain eigenvalues
            magnitude         = sqrt((tensor.j22 - tensor.j11).^2 + 4 .* tensor.j12.^2);
            minimumEigenvalue = 0.5 .* (tensor.j11 + tensor.j22 - magnitude);
            maximumEigenvalue = 0.5 .* (tensor.j11 + tensor.j22 + magnitude);
            
            % Obtain condition number
            value = maximumEigenvalue ./ minimumEigenvalue;
        end
        
        function value = obtainMininmumEigenvalue(tensor)
            %
            % Obtain minimum eigenvalue.
            %
            % INPUTS:
            %   1. data - data to compute coherence measurement.
            %
            narginchk(1,1);
            
            % Obtain minimum eigenvalues
            magnitude         = sqrt((tensor.j22 - tensor.j11).^2 + 4 .* tensor.j12.^2);
            minimumEigenvalue = 0.5 .* (tensor.j11 + tensor.j22 - magnitude);
            value = minimumEigenvalue;
        end            
    end
        
    enumeration
        ANISOTROPY ('anysotropy'     , [ 'sqrt((J22 - J11)^2 + 4 * J12^2) / (J22 + J11)^2 = ' ...
                                       , '((eig_max - eig_min) / (eig_max - eig_min))^2' ] ...
                   , @(data) image.gradient.enums.CoherenceMeasurements.obtainAnisotropy(data) );
        CONTINUOUS ('continuous'     , [ 'sqrt((J22 - J11)^2 + 4 * J12^2) = ' ...
                                       , 'eig_max - eig_min' ] ...
                   , @(data) image.gradient.enums.CoherenceMeasurements.obtainContinuous(data) );
        FORSTNER   ('forstner'       , '(J11 * J22 - J12^2) / (J11 + J22)' ...
                   , @(data) image.gradient.enums.CoherenceMeasurements.obtainForstner(data) );
        SHI_TOMASI ('shiTomasi' , [ '0.5 * (J11 + J22 - sqrt((J22 - J11)^2 + 4 * J12^2)) =' ...
                                  , 'min(eig_max,eig_min)' ] ...
                   , @(data) image.gradient.enums.CoherenceMeasurements.obtainMininmumEigenvalue(data) );
        CORNER_RESPONSE   ('cornerResponse' , '(J11 * J22 - J12^2) - 0.04 * (J11 + J22)' ...
                          , @(data) image.gradient.enums.CoherenceMeasurements.obtainCornerResponseFunction(data) );
        CONDITION_NUMBER  ('conditionNumber', [ '(J11 * J22 - J12^2) - 0.04 * (J11 + J22) = ' ...
                                              , 'eig_max / eig_min' ] ...
                          , @(data) image.gradient.enums.CoherenceMeasurements.obtainConditionNumber(data) );
    end
end