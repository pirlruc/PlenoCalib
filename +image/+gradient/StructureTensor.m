classdef StructureTensor
    %STRUCTURETENSOR 
    %   2D Structure tensor properties.
    
    properties
        data  = []      % Structure tensor data
    end
    
    properties (Dependent)
        j11             % Entry (1,1) of structure tensor
        j12             % Entry (1,2) of structure tensor
        j21             % Entry (2,1) of structure tensor
        j22             % Entry (2,2) of structure tensor
    end
    
    methods
        function self = StructureTensor(varargin)
            %
            % Create structure tensor instance.
            %
            % INPUTS:
            %   1. data  - structure tensor data.
            %
            narginchk(0,1);
            
            if ~isempty(varargin)
                if nargin >= 1
                    self.data = varargin{1};
                end
            end
        end
        
        function self = set.data(self, newStructureTensorData)
            self.data = newStructureTensorData;
        end
        
        function entry = get.j11(self)
            % If data is not provided, return NaN
            if isempty(self.data)
                entry = nan;
            else
                entry = self.data(1,1);
            end
        end
        
        function entry = get.j12(self)
            % If data is not provided, return NaN
            if isempty(self.data)
                entry = nan;
            else
                entry = self.data(1,2);
            end
        end
        
        function entry = get.j21(self)
            % If data is not provided, return NaN
            if isempty(self.data)
                entry = nan;
            else
                entry = self.data(2,1);
            end
        end
        
        function entry = get.j22(self)
            % If data is not provided, return NaN
            if isempty(self.data)
                entry = nan;
            else
                entry = self.data(2,2);
            end
        end
    end
    
    methods 
        function [coherence,phase,eigenvectors,eigenvalues,magnitude] = ...
                            decompose(self,coherenceMeasurement)
            %
            % Perform eigenvalue decomposition for the structure tensor.
            %
            % INPUTS:
            %   1. coherenceMeasurement - coherence measurement type to be
            %   used. If no coherence measurement is provided, assume the
            %   standard coherence measurement.
            %
            narginchk(1,2);
            
            if nargin <= 1
                coherenceMeasurement = image.gradient.enums.CoherenceMeasurements.CONTINUOUS();
            end
            
            % Throw error if data is not defined
            if isempty(self.data)
                error = MException( 'StructureTensor:decompose:noData' ...
                                  , 'Structure tensor data not defined...' );
                error.throw();
            end
            
            % Decompose structure tensor
            [rightEigenvectorsMatrix,eigenvaluesMatrix] = eig(self.data);

            % Select indices to obtain eigenvalues from eigenvalues matrix
            eigenvaluesIndices = sub2ind( size(eigenvaluesMatrix) ...
                                        , 1:size(eigenvaluesMatrix,1) ...
                                        , 1:size(eigenvaluesMatrix,2) );
            eigenvalues        = eigenvaluesMatrix(eigenvaluesIndices);
            
            % Obtain sorted indices and eigenvalues. This way we have the
            % eigenvalues sorted from minimum to maximum.
            [eigenvalues,indices] = sort(eigenvalues);

            % Apply sorting to eigenvectors
            eigenvectors = math.Vector(rightEigenvectorsMatrix(:,indices));
            
            % Obtain phase of structure tensor
            % Obtain components of eigenvectors
            y = eigenvectors.obtainComponents(2).data;
            x = eigenvectors.obtainComponents(1).data;
                
            % Obtain phase associated with each eigenvector
            phase = atan2(y,x)';

            % Give the eigenvalues, eigenvectors and phase in a structured 
            % form
            eigenvalues  = struct( 'minimum', eigenvalues(1) ...
                                 , 'maximum', eigenvalues(2) );
            eigenvectors = struct( 'minimum', eigenvectors.obtainVectors(1) ...
                                 , 'maximum', eigenvectors.obtainVectors(2) );
            phase = struct( 'minimum', phase(1) ...
                          , 'maximum', phase(2) );

            % Obtain magnitude of structure tensor
            %     magnitude = lambda_max - lambda_min
            magnitude    = eigenvalues.maximum - eigenvalues.minimum;
            
            % Obtain coherence measurement
            coherence = coherenceMeasurement.coherence(self);
        end
    end
    
    methods (Static)
        function self = StructureTensorFromStructureTensorField( structureTensor ...
                                                               , pixelInStructureTensor )
            %
            % Create structure tensor instance from structure tensor field.
            %
            % INPUTS:
            %   1. tensor - structure tensor instance.
            %   2. pixel  - pixel position.
            %
            narginchk(2,2);
            
            % Convert pixel data to pixel object
            pixelInStructureTensor = utils.enums.Classes.PIXEL().convert(pixelInStructureTensor);
            
            % Create structure tensor instance. The structure tensor field
            % has structure:
            %       channels x pixels_u x pixels_v x dimension x dimension 
            % This structure is modified to have a structure tensor with:
            %       dimension x dimension 
            self = image.gradient.StructureTensor( ...
                            permute( structureTensor.data( :, pixelInStructureTensor.u ...
                                                            , pixelInStructureTensor.v,:,:) ...
                                   , [4,5,2,3,1] ));
        end
    end
end
