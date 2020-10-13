classdef DistortionModels
    %DISTORTIONMODELS
    %   Enumerate distortion models that can be used for the lytro camera.
    
    properties
        keyword           % Distortion model keyword
        description       % Distortion model description
        firstParameters   % Parameters to optimize in first estimation step
        secondParameters  % Parameters to optimize in second estimation step
        refineParameters  % Parameters to optimize in final estimation step
    end
    
    methods
        function self = DistortionModels( keyword, description ...
                                        , firstParameters, secondParameters, refineParameters )
            %
            % Distortion model instance.
            %
            self.keyword          = keyword;
            self.description      = description;
            self.firstParameters  = firstParameters;
            self.secondParameters = secondParameters;
            self.refineParameters = refineParameters;
        end
    end
    
    methods
        function parameters = initializeParameters(self)
            %
            % Initialize distortion parameters.
            %
            parameters = zeros(1,length(self.refineParameters));
        end
        
        function parameters = encode(self,distortion)
            %
            % Encode distortion model into a vector of parameters.
            %
            % INPUTS: 
            %   1. distortion - distortion model data.
            %
            narginchk(2,2);
            
            % Obtain parameters from distortion model
            if self == camera.models.plenoptic.enums.DistortionModels.DANSEREAU_ORIGINAL()
                parameters = [distortion.radial(1:3);distortion.center.data];
            elseif self == camera.models.plenoptic.enums.DistortionModels.DANSEREAU_SIMPLIFIED()
                parameters = distortion.radial(1:3);
            elseif self == camera.models.plenoptic.enums.DistortionModels.BOK_ORIGINAL()
                parameters = distortion.radial(1:3);
            elseif self == camera.models.plenoptic.enums.DistortionModels.BOK_DISTORTION_CENTER()
                parameters = [distortion.radial(1:3);distortion.center.data];
            elseif self == camera.models.plenoptic.enums.DistortionModels.JOHANNSEN_ORIGINAL()
                parameters = [distortion.radial(1:3);distortion.tangential;distortion.center.data];
            elseif self == camera.models.plenoptic.enums.DistortionModels.ZELLER_ORIGINAL()
                parameters = [distortion.radial(1:3);distortion.tangential];
            elseif self == camera.models.plenoptic.enums.DistortionModels.ZHANG_ORIGINAL()
                parameters = [distortion.radial(1:3);distortion.scaling.data(:,1)];
            elseif self == camera.models.plenoptic.enums.DistortionModels.ZHANG_DISTORTION_CENTER()
                parameters = [distortion.radial(1:3);distortion.scaling.data(:,1);distortion.center.data];
            elseif self == camera.models.plenoptic.enums.DistortionModels.MONTEIRO_ORIGINAL()
                parameters = [distortion.radial(1:3);distortion.center.data;distortion.scaling.data(:,1);distortion.scaling.data(:,2)];
            end
            parameters = parameters';
        end
        
        function distortion = decode(self,parameters)
            %
            % Decode vector of parameters into a distortion model.
            %
            % INPUTS: 
            %   1. parameters - parameters of distortion model.
            %
            narginchk(2,2);

            % Obtain distortion model from parameters
            distortion = camera.models.plenoptic.Distortion();
            distortion.distortionModel = self;
            if ~isempty(parameters)
                if self == camera.models.plenoptic.enums.DistortionModels.DANSEREAU_ORIGINAL()
                    distortion.radial = [parameters(1:3),0]';
                    distortion.center = parameters(4:5)';
                elseif self == camera.models.plenoptic.enums.DistortionModels.DANSEREAU_SIMPLIFIED()
                    distortion.radial = [parameters(1:3),0]';
                elseif self == camera.models.plenoptic.enums.DistortionModels.BOK_ORIGINAL()
                    distortion.radial = [parameters(1:3),0]';
                elseif self == camera.models.plenoptic.enums.DistortionModels.BOK_DISTORTION_CENTER()
                    distortion.radial = [parameters(1:3),0]';
                    distortion.center = parameters(4:5)';
                elseif self == camera.models.plenoptic.enums.DistortionModels.JOHANNSEN_ORIGINAL()
                    distortion.radial     = [parameters(1:3),0]';
                    distortion.tangential = parameters(4:6)';
                    distortion.center     = parameters(7:8)';
                elseif self == camera.models.plenoptic.enums.DistortionModels.ZELLER_ORIGINAL()
                    distortion.radial     = [parameters(1:3),0]';
                    distortion.tangential = parameters(4:6)';
                elseif self == camera.models.plenoptic.enums.DistortionModels.ZHANG_ORIGINAL()
                    distortion.radial  = [parameters(1:3),0]';
                    distortion.scaling = [parameters(4),0;parameters(5),0];
                elseif self == camera.models.plenoptic.enums.DistortionModels.ZHANG_DISTORTION_CENTER()
                    distortion.radial  = [parameters(1:3),0]';
                    distortion.scaling = [parameters(4),0;parameters(5),0];
                    distortion.center  = parameters(6:7)';
                elseif self == camera.models.plenoptic.enums.DistortionModels.MONTEIRO_ORIGINAL()
                    distortion.radial  = [parameters(1:3),0]';
                    distortion.center  = parameters(4:5)';
                    distortion.scaling = [parameters(6:7)',parameters(8:9)'];
                end
            end
        end
    end
    
    enumeration
        DANSEREAU_ORIGINAL     ( 'directions_radialdistort_center'   ...
                               , [ 'Radial distortion with optimization of distortion center. ' ...
                                 , 'The radial distortion is applied to directions (u,v).' ] ...
                               , 1:5, [], 1:5 );
        DANSEREAU_SIMPLIFIED   ( 'directions_radialdistort_nocenter' ...
                               , [ 'Radial distortion with distortion center at (0,0). ' ...
                                 , 'The radial distortion is applied to directions (u,v).' ] ...
                               , 1:3, [], 1:3 );
        BOK_ORIGINAL           ( 'positions_directions_radialdistort_nocenter' ...
                               , [ 'Radial distortion with distortion center at (0,0). ' ...
                                 , 'The radial distortion is applied to directions (u,v) and positions (s,t).' ] ...
                               , 1:3, [], 1:3 );
        BOK_DISTORTION_CENTER  ( 'positions_directions_radialdistort_center' ...
                               , [ 'Radial distortion with optimization of distortion centers. ' ...
                                 , 'The radial distortion is applied to directions (u,v) and positions (s,t).' ] ...
                               , 1:5, [], 1:5 );
        JOHANNSEN_ORIGINAL     ( 'directions_distort_center' ...
                               , [ 'Radial and tangential distortion with optimization of distortion center. ' ...
                                 , 'The distortion is applied to directions (u,v).' ] ...
                               , [1:3,7,8], 4:6, 1:8 );
        ZELLER_ORIGINAL        ( 'directions_distort_nocenter' ...
                               , [ 'Radial and tangential distortion with distortion center at (0,0). ' ...
                                 , 'The distortion is applied to directions (u,v).' ] ...
                               , 1:6, [], 1:6 );
        ZHANG_ORIGINAL         ( 'positions_directions_mixture_radialdistort_nocenter' ...
                               , [ 'Radial distortion with distortion center at (0,0). ' ...
                                 , 'The radial distortion is applied to directions (u,v) ' ...
                                 , 'considering an offset proportional to the position (s,t).' ] ...
                               , 1:5, [], 1:5 );
        ZHANG_DISTORTION_CENTER( 'positions_directions_mixture_radialdistort_center' ...
                               , [ 'Radial distortion with optimization of distortion center. ' ...
                                 , 'The radial distortion is applied to directions (u,v) ' ...
                                 , 'considering an offset proportional to the position (s,t).' ] ...
                               , 1:7, [], 1:7 );
        MONTEIRO_ORIGINAL      ( 'positions_directions_mixture_radialdistort_center_dansereau_tip' ...
                               , [ 'Radial distortion with optimization of distortion center. ' ...
                                 , 'The radial distortion is applied to directions (u,v) ' ...
                                 , 'considering a dependency on higher order terms of the position (s,t).' ] ...
                               , 1:9, [], 1:9 );
    end
end

