classdef MRP
    %MRP Modified Rodrigues Parameters Class
    %   A class with the basic methods corresponding to the MRPs.
    
    properties
    end
    
    methods (Static)
        function rates = getRates(mrp, w)
            %getRates Construct the rates matric from the DCM and the angular
            %rates vector
            %   Detailed explanation goes here
            rates = (1/4) * ((1 - mrp' * mrp)* eye(3) + 2 * CrossMatrix(mrp) + 2 * (mrp * mrp')) * w;
        end
        
        function addition = method1(DCM1,DCM2)
            %addition Addition of DCMs
            %   AC = AB * BC
            addition = DCM1 * DCM2;
        end
    end
end

