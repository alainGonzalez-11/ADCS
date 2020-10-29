classdef DCM
    %DCM Director Cosines Matrix Class
    %   A class with the basic methods corresponding to the DCM.
    
    properties
    end
    
    methods (Static)
        function rates = getRates(DCM, w)
            %getRates Construct the rates matric from the DCM and the angular
            %rates vector
            %   Detailed explanation goes here
            w = CrossMatrix(w);
            rates = -w * DCM;
        end
        
        function addition = method1(DCM1,DCM2)
            %addition Addition of DCMs
            %   AC = AB * BC
            addition = DCM1 * DCM2;
        end
    end
end

