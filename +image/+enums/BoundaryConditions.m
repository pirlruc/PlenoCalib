classdef BoundaryConditions
    %BOUNDARYCONDITIONS
    %   Filter boundary conditions.
    
    properties (~Dependent, SetAccess = immutable, GetAccess = public)
        method          % Boundary condition keyword
        description     % Description of boundary condition
    end
    
    methods
        function self = BoundaryConditions(method, description)
            %
            % Boundary conditions enumeration instance.
            %
            self.method      = method;
            self.description = description;
        end
    end
    
    enumeration
        SYMMETRIC ('symmetric', [ 'Input array values outside the bounds of the array are computed by ' ...
                                  'mirror-reflecting the array across the array border.' ]);
        REPLICATE ('replicate', [ 'Input array values outside the bounds of the array are assumed to ' ...
                                  'equal the nearest array border value.' ]);
        PERIODIC  ('circular' , [ 'Input array values outside the bounds of the array are computed by ' ...
                                  'implicitly assuming the input array is periodic.' ]);
    end
end