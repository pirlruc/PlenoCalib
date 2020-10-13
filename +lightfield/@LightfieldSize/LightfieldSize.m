classdef LightfieldSize < abstract.TemplateSize
    %LIGHTFIELDSIZE
    %   Utility to define the lightfield size structure.
    
    properties (Dependent)
        numberPixels_i              % Number of horizontal pixels in microlens
        numberPixels_j              % Number of vertical pixels in microlens
        numberMicrolenses_k         % Number of horizontal microlenses
        numberMicrolenses_l         % Number of vertical microlenses
        numberChannels              % Lightfield number of channels
    end
    
    methods
        function self = LightfieldSize(varargin)
            %
            % Lightfield size instance.
            %
            % INPUTS:
            %   1. ligthfieldSize - lightfield size array.
            %
            narginchk(0,1)
            
            self@abstract.TemplateSize(varargin{:});
        end
        
        function number = get.numberPixels_i(self)
            number = self.numberItemsInDimension(2);
        end
        
        function number = get.numberPixels_j(self)
            number = self.numberItemsInDimension(3);
        end
        
        function number = get.numberMicrolenses_k(self)
            number = self.numberItemsInDimension(4);
        end
        
        function number = get.numberMicrolenses_l(self)
            number = self.numberItemsInDimension(5);
        end
        
        function number = get.numberChannels(self)
            number = self.numberItemsInDimension(1);
        end
    end
end
